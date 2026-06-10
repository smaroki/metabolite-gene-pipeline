import pandas as pd
import gc
import numpy as np
from pathlib import Path
from matchms.importing import load_from_mgf
from matchms.filtering import (
    default_filters,
    normalize_intensities,
    select_by_mz,
    add_losses,
    require_minimum_number_of_peaks
)
from ms2query import MS2Library
from ms2query.utils import SettingsRunMS2Query

# ============================================================================
# CONFIGURATION
# ============================================================================
INPUT_MGF   = Path(r"C:\Users\palokaich\Desktop\FF2.mgf")
OUTPUT_CSV  = Path(r"D:\ms2query_results\Analogues.csv")
LIBRARY_DIR = Path(r"D:\downloads\ms2query_lib-20260610T102616Z-3-001\ms2query_lib")

# SEARCH PARAMETERS
MIN_PEAKS      = 3
MZ_MIN         = 0
MZ_MAX         = 2000
TOP_N_RESULTS  = 20

print("=" * 80)
print("MS2Query: ANALOG SEARCH")
print("=" * 80)

# ============================================================================
# STEP 1: Validate Library Directory
# ============================================================================
print("\n[1/6] Validating library directory...")

required_files = {
    "sqlite_file_name"           : "ms2query_library.sqlite",
    "s2v_model_file_name"        : "spec2vec_model.model",
    "ms2ds_model_file_name"      : "ms2deepscore_model.pt",
    "ms2query_model_file_name"   : "ms2query_model.onnx",
    "s2v_embeddings_file_name"   : "s2v_embeddings.parquet",
    "ms2ds_embeddings_file_name" : "ms2ds_embeddings.parquet",
}

if not LIBRARY_DIR.exists():
    raise FileNotFoundError(
        f"Library directory not found:\n  {LIBRARY_DIR}\n"
        f"Please check the path and ensure the drive is mounted."
    )

missing = []
for arg, fname in required_files.items():
    fpath = LIBRARY_DIR / fname
    if not fpath.exists():
        missing.append(fname)
    else:
        size_mb = fpath.stat().st_size / (1024 ** 2)
        print(f"   ✓ {fname}  ({size_mb:.1f} MB)")

if missing:
    raise FileNotFoundError(
        f"The following required library files are missing:\n"
        + "\n".join(f"   - {f}" for f in missing)
    )

print("   All library files present.")

# ============================================================================
# STEP 2: Load and Preprocess Spectra
# ============================================================================
print(f"\n[2/6] Loading spectra from:\n   {INPUT_MGF}")

if not INPUT_MGF.exists():
    raise FileNotFoundError(f"MGF input file not found:\n  {INPUT_MGF}")

spectra = list(load_from_mgf(str(INPUT_MGF)))
print(f"   Loaded {len(spectra)} spectra")

if len(spectra) > 0:
    first = spectra[0]
    print(f"\n   First spectrum preview:")
    print(f"      Precursor m/z : {first.metadata.get('precursor_mz', 'N/A')}")
    print(f"      Number of peaks: {len(first.peaks.mz)}")

print(f"\n[3/6] Preprocessing spectra...")
print(f"   Parameters: min_peaks={MIN_PEAKS}, mz_range=[{MZ_MIN}, {MZ_MAX}]")

processed = []
skipped   = 0

for i, spec in enumerate(spectra):
    try:
        s = default_filters(spec)
        if s is None:
            skipped += 1
            continue
        s = normalize_intensities(s)
        s = select_by_mz(s, mz_from=MZ_MIN, mz_to=MZ_MAX)
        if s is None:
            skipped += 1
            continue
        s = require_minimum_number_of_peaks(s, n_required=MIN_PEAKS)
        if s is None:
            skipped += 1
            continue
        s = add_losses(s)
        processed.append(s)
    except Exception as e:
        print(f"   Warning: Spectrum {i} failed preprocessing: {e}")
        skipped += 1

print(f"   Passed : {len(processed)}/{len(spectra)} spectra")
print(f"   Skipped: {skipped} spectra (too few peaks or filtering failure)")

if len(processed) == 0:
    raise ValueError(
        "No spectra survived preprocessing. "
        "Try lowering MIN_PEAKS or widening the MZ_MIN/MZ_MAX range."
    )

# Free raw spectra from memory
del spectra
gc.collect()

# ============================================================================
# STEP 3: Load MS2Query Library
# ============================================================================
print(f"\n[4/6] Loading MS2Query library...")

ms2library = MS2Library(
    sqlite_file_name           = str(LIBRARY_DIR / "ms2query_library.sqlite"),
    s2v_model_file_name        = str(LIBRARY_DIR / "spec2vec_model.model"),
    ms2ds_model_file_name      = str(LIBRARY_DIR / "ms2deepscore_model.pt"),
    ms2query_model_file_name   = str(LIBRARY_DIR / "ms2query_model.onnx"),
    s2v_embeddings_file_name   = str(LIBRARY_DIR / "s2v_embeddings.parquet"),
    ms2ds_embeddings_file_name = str(LIBRARY_DIR / "ms2ds_embeddings.parquet"),
)

print("   Library loaded successfully.")

# ============================================================================
# STEP 4: Configure Search Settings
# ============================================================================
print("\n[5/6] Configuring search settings...")

settings = SettingsRunMS2Query(
    nr_of_top_analogs_to_save   = TOP_N_RESULTS,  # 20, from your config
    minimal_ms2query_metascore  = 0,              # retain all matches
    preselection_cut_off        = MZ_MAX,         # 2000, from your config
)

print(f"   Top analogs to save     : {TOP_N_RESULTS}")
print(f"   Minimal metascore       : 0 (all matches retained)")
print(f"   Preselection m/z cutoff : {MZ_MAX}")


# ============================================================================
# STEP 5: Run Analog Search
# ============================================================================
print(f"\n[6/6] Running analog search on {len(processed)} spectra...")
print("   This may take several minutes depending on dataset size.\n")

all_results = []

for i, df_result in enumerate(ms2library.analog_search_yield_df(
    query_spectra=processed,
    settings=settings,
    progress_bar=True
)):
    all_results.append(df_result)
    if (i + 1) % 10 == 0:
        print(f"   Progress: {i + 1}/{len(processed)} spectra processed")

print(f"\n   Search complete. Results collected for {len(all_results)} spectra.")

# ============================================================================
# STEP 6: Combine, Format, and Save Results
# ============================================================================
print("\n[Post-processing] Formatting and saving results...")

if len(all_results) == 0:
    print("WARNING: No results returned. Check your MGF file and library compatibility.")
else:
    df_combined = pd.concat(all_results, ignore_index=True)

    print(f"   Total raw matches : {len(df_combined)}")
    print(f"   Columns available : {list(df_combined.columns)}")

    # ── Inspect for SMILES / InChIKey columns ────────────────────────────
    smiles_candidates  = [c for c in df_combined.columns if 'smiles'  in c.lower()]
    inchikey_candidates= [c for c in df_combined.columns if 'inchi'   in c.lower()]

    if smiles_candidates:
        print(f"   SMILES columns found   : {smiles_candidates}")
    else:
        print("   WARNING: No SMILES columns detected in results.")

    if inchikey_candidates:
        print(f"   InChIKey columns found : {inchikey_candidates}")
    else:
        print("   WARNING: No InChIKey columns detected in results.")

    # ── Sample row for debugging ──────────────────────────────────────────
    print(f"\n   Sample row (first match):")
    first_row = df_combined.iloc[0]
    for col in df_combined.columns:
        val = first_row[col]
        if isinstance(val, str) and len(val) > 100:
            val = val[:100] + "..."
        print(f"      {col}: {val}")

    # ── Build output rows ─────────────────────────────────────────────────
    output_data = []

    for _, row in df_combined.iterrows():

        spectrum_id   = row.get('query_spectrum_nr',
                         row.get('spectrum_id', 'unknown'))
        precursor_mz  = row.get('precursor_mz_query',
                         row.get('precursor_mz', 'N/A'))
        score         = row.get('ms2query_model_prediction',
                         row.get('analog_score',
                         row.get('score', 0.0)))
        name          = row.get('analog_compound_name',
                         row.get('compound_name',
                         row.get('name', 'Unknown')))

        # SMILES — try multiple possible column names
        smiles = ''
        for col in ['smiles', 'analog_smiles', 'canonical_smiles',
                    'isomeric_smiles', 'SMILES']:
            if col in row.index and pd.notna(row[col]):
                smiles = str(row[col])
                break

        # InChIKey — try multiple possible column names
        inchikey = ''
        for col in ['inchikey', 'analog_inchikey', 'inchikey_14',
                    'InChIKey', 'inchi_key']:
            if col in row.index and pd.notna(row[col]):
                inchikey = str(row[col])
                break

        formula            = row.get('formula', row.get('molecular_formula', ''))
        precursor_mz_match = row.get('precursor_mz_library',
                              row.get('library_precursor_mz', ''))

        output_data.append({
            'spectrum_id'       : spectrum_id,
            'query_precursor_mz': precursor_mz,
            'analog_score'      : round(float(score), 4),
            'compound_name'     : name,
            'smiles'            : smiles,
            'inchikey'          : inchikey,
            'formula'           : formula,
            'match_precursor_mz': precursor_mz_match,
        })

    df_output = pd.DataFrame(output_data)

    # ── SMILES coverage report ────────────────────────────────────────────
    smiles_count = (df_output['smiles'] != '').sum()
    print(f"\n   Rows with SMILES   : {smiles_count}/{len(df_output)}")

    # ── Sort and rank ─────────────────────────────────────────────────────
    df_output = df_output.sort_values(
        ['spectrum_id', 'analog_score'],
        ascending=[True, False]
    )
    df_output['rank'] = df_output.groupby('spectrum_id').cumcount() + 1
    df_output = df_output[df_output['rank'] <= TOP_N_RESULTS]

    # ── Save ──────────────────────────────────────────────────────────────
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df_output.to_csv(OUTPUT_CSV, index=False)

    # ── Summary statistics ────────────────────────────────────────────────
    print(f"\n{'=' * 80}")
    print(f"ANALYSIS COMPLETE")
    print(f"{'=' * 80}")
    print(f"\n   Output file      : {OUTPUT_CSV}")
    print(f"   Spectra analyzed : {len(processed)}")
    print(f"   Spectra matched  : {df_output['spectrum_id'].nunique()}")
    print(f"   Total rows saved : {len(df_output)}")

    print(f"\n   Score statistics:")
    print(f"      Mean   : {df_output['analog_score'].mean():.4f}")
    print(f"      Median : {df_output['analog_score'].median():.4f}")
    print(f"      Max    : {df_output['analog_score'].max():.4f}")
    print(f"      Min    : {df_output['analog_score'].min():.4f}")

    print(f"\n   Top 5 matches by score:")
    top5 = df_output.nlargest(5, 'analog_score')
    for i, (_, row) in enumerate(top5.iterrows(), 1):
        print(f"      {i}. {str(row['compound_name'])[:55]}")
        print(f"         Score: {row['analog_score']:.4f}  |  "
              f"Spectrum: {row['spectrum_id']}")

    print(f"\n   Sample output (first 10 rows):")
    display_cols   = ['spectrum_id', 'compound_name', 'analog_score', 'smiles', 'inchikey']
    available_cols = [c for c in display_cols if c in df_output.columns]
    print(df_output[available_cols].head(10).to_string(index=False, max_colwidth=60))

    print(f"\n{'=' * 80}")
    print(f"Done! Results saved to: {OUTPUT_CSV}")
    print(f"{'=' * 80}")

# ============================================================================
# ALTERNATIVE: Direct CSV Export (simpler, no post-processing control)
# Uncomment below and comment out Steps 5-6 above to use this instead.
# ============================================================================
# print("\n[Alternative] Direct CSV export...")
# ms2library.analog_search_store_in_csv(
#     query_spectra=processed,
#     results_csv_file_location=str(OUTPUT_CSV),
#     settings=settings
# )
# print(f"Results saved directly to: {OUTPUT_CSV}")
