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
from ms2query import create_library_object_from_one_dir
from ms2query.utils import SettingsRunMS2Query

# ============================================================================
# CONFIGURATION
# ============================================================================
INPUT_MGF = Path(r"C:\Users\dell\OneDrive\Desktop\Raw files\FM2B5short.mgf")
OUTPUT_CSV = Path(r"C:\Users\dell\Desktop\Final Manuscript Data\Original pipeline output\FM2B5_predictions.csv")
LIBRARY_DIR = Path(r"C:\Users\dell\OneDrive\Desktop\Raw files\ms2query_lib")

# RELAXED PARAMETERS
MIN_PEAKS = 3
MZ_MIN = 0
MZ_MAX = 2000
TOP_N_RESULTS = 20

print("=" * 80)
print("MS2Query: ANALOG SEARCH (Fixed API)")
print("=" * 80)

# ============================================================================
# STEP 1: Load Spectra
# ============================================================================
print("\n[1/5] Loading spectra...")
spectra = list(load_from_mgf(str(INPUT_MGF)))
print(f"Loaded {len(spectra)} spectra")

if len(spectra) > 0:
    first = spectra[0]
    print(f"\nFirst spectrum:")
    print(f"   Precursor m/z: {first.metadata.get('precursor_mz', 'N/A')}")
    print(f"   Peaks: {len(first.peaks.mz)}")

# ============================================================================
# STEP 2: Preprocess Spectra
# ============================================================================
print(f"\n[2/5] Preprocessing...")

processed = []
for i, spec in enumerate(spectra):
    try:
        s = default_filters(spec)
        if s is None:
            continue
        s = normalize_intensities(s)
        s = select_by_mz(s, mz_from=MZ_MIN, mz_to=MZ_MAX)
        if s is None:
            continue
        s = require_minimum_number_of_peaks(s, n_required=MIN_PEAKS)
        if s is None:
            continue
        s = add_losses(s)
        processed.append(s)
    except Exception as e:
        print(f"   Warning: Spectrum {i} failed: {e}")

print(f"{len(processed)}/{len(spectra)} spectra passed preprocessing")

if len(processed) == 0:
    raise ValueError("No spectra survived preprocessing!")

del spectra
gc.collect()

# ============================================================================
# STEP 3: Load Library
# ============================================================================
print("\n[3/5] Loading MS2Query library...")
ms2library = create_library_object_from_one_dir(str(LIBRARY_DIR))
print(f"Library loaded: {ms2library.ionization_mode} mode")

# ============================================================================
# STEP 4: Configure Search Settings (OPTIONAL)
# ============================================================================
print("\n[4/5] Configuring search...")

# Use default settings or customize if needed
# Common parameters: ms2query_score_cutoff, ms2ds_score_cutoff
settings = None  # Use defaults for now

# If you want to filter by score later, uncomment:
# settings = SettingsRunMS2Query(
#     ms2query_score_cutoff=0.5,  # Minimum score threshold
# )

print(f"   Using default settings (will get all matches)")

# ============================================================================
# STEP 5: Run Analog Search - THE CORRECT WAY!
# ============================================================================
print(f"\n[5/5] Running analog search on {len(processed)} spectra...")
print("   (This may take a few minutes...)\n")

# Method 1: Use analog_search_yield_df (recommended for processing)
all_results = []

for i, df_result in enumerate(ms2library.analog_search_yield_df(
    query_spectra=processed,
    settings=settings,
    progress_bar=True
)):
    # df_result is a DataFrame with results for one spectrum
    all_results.append(df_result)
    
    if (i + 1) % 10 == 0:
        print(f"   Processed {i+1}/{len(processed)} spectra")

print(f"\nSearch complete!")

# ============================================================================
# STEP 6: Combine and Format Results
# ============================================================================
print("\n[6/6] Formatting results...")

if len(all_results) == 0:
    print("WARNING: No results found!")
else:
    # Combine all DataFrames
    df_combined = pd.concat(all_results, ignore_index=True)
    
    print(f"   Total matches: {len(df_combined)}")
    print(f"   Columns: {list(df_combined.columns)}")
    
    # Show what columns are available
    print(f"\nAvailable result columns:")
    for col in df_combined.columns:
        print(f"   - {col}")
    
    # First, let's see what columns we actually have
    print(f"\nInspecting columns for SMILES/InChIKey...")
    print(f"   All columns: {list(df_combined.columns)}")
    
    # Look for SMILES in various possible column names
    smiles_candidates = [col for col in df_combined.columns if 'smiles' in col.lower()]
    inchikey_candidates = [col for col in df_combined.columns if 'inchi' in col.lower()]
    
    if smiles_candidates:
        print(f"   Found SMILES columns: {smiles_candidates}")
    if inchikey_candidates:
        print(f"   Found InChIKey columns: {inchikey_candidates}")
    
    # Show sample row to debug
    if len(df_combined) > 0:
        print(f"\nSample row (first match):")
        first_row = df_combined.iloc[0]
        for col in df_combined.columns:
            val = first_row[col]
            if isinstance(val, str) and len(val) > 100:
                val = val[:100] + "..."
            print(f"   {col}: {val}")
    
    # Format output - limit to TOP_N_RESULTS per spectrum
    output_data = []
    
    for _, row in df_combined.iterrows():
        # Get query spectrum info
        spectrum_id = row.get('query_spectrum_nr', row.get('spectrum_id', 'unknown'))
        precursor_mz = row.get('precursor_mz_query', row.get('precursor_mz', 'N/A'))
        
        # Get match info - try multiple column names
        score = row.get('ms2query_model_prediction', row.get('analog_score', row.get('score', 0.0)))
        name = row.get('analog_compound_name', row.get('compound_name', row.get('name', 'Unknown')))
        
        # Try to get SMILES from different possible columns
        smiles = ''
        for smiles_col in ['smiles', 'analog_smiles', 'canonical_smiles', 'isomeric_smiles', 'SMILES']:
            if smiles_col in row.index and pd.notna(row[smiles_col]):
                smiles = str(row[smiles_col])
                break
        
        # Try to get InChIKey from different possible columns  
        inchikey = ''
        for inchi_col in ['inchikey', 'analog_inchikey', 'inchikey_14', 'InChIKey', 'inchi_key']:
            if inchi_col in row.index and pd.notna(row[inchi_col]):
                inchikey = str(row[inchi_col])
                break
        
        # Get additional metadata if available
        formula = row.get('formula', row.get('molecular_formula', ''))
        precursor_mz_match = row.get('precursor_mz_library', row.get('library_precursor_mz', ''))
        
        output_data.append({
            'spectrum_id': spectrum_id,
            'query_precursor_mz': precursor_mz,
            'analog_score': round(float(score), 4),
            'compound_name': name,
            'smiles': smiles,
            'inchikey': inchikey,
            'formula': formula,
            'match_precursor_mz': precursor_mz_match,
        })
    
    df_output = pd.DataFrame(output_data)
    
    # Check if we got any SMILES
    smiles_count = df_output['smiles'].notna().sum()
    print(f"\n   Rows with SMILES: {smiles_count}/{len(df_output)}")
    
    # Sort by spectrum and score
    df_output = df_output.sort_values(
        ['spectrum_id', 'analog_score'], 
        ascending=[True, False]
    )
    
    # Add rank and keep only TOP_N_RESULTS per spectrum
    df_output['rank'] = df_output.groupby('spectrum_id').cumcount() + 1
    df_output = df_output[df_output['rank'] <= TOP_N_RESULTS]
    
    # Save
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df_output.to_csv(OUTPUT_CSV, index=False)
    
    # ========================================================================
    # STATISTICS
    # ========================================================================
    print(f"\n{'=' * 80}")
    print(f"ANALYSIS COMPLETE")
    print(f"{'=' * 80}")
    
    print(f"\nOutput: {OUTPUT_CSV}")
    print(f"Total matches: {len(df_output)}")
    print(f"Spectra analyzed: {len(processed)}")
    print(f"Spectra with matches: {df_output['spectrum_id'].nunique()}")
    
    print(f"\nScore Statistics:")
    print(f"   Mean:   {df_output['analog_score'].mean():.4f}")
    print(f"   Median: {df_output['analog_score'].median():.4f}")
    print(f"   Max:    {df_output['analog_score'].max():.4f}")
    print(f"   Min:    {df_output['analog_score'].min():.4f}")
    
    print(f"\nTOP 5 MATCHES:")
    top5 = df_output.nlargest(5, 'analog_score')
    for i, (_, row) in enumerate(top5.iterrows(), 1):
        print(f"   {i}. {row['compound_name'][:50]}")
        print(f"      Score: {row['analog_score']:.4f} | Spectrum: {row['spectrum_id']}")
    
    print(f"\nSample Results (with SMILES):")
    display_cols = ['spectrum_id', 'compound_name', 'analog_score', 'smiles', 'inchikey']
    available_cols = [col for col in display_cols if col in df_output.columns]
    print(df_output[available_cols].head(10).to_string(index=False, max_colwidth=60))
    
    print(f"\nDone!")

# ============================================================================
# ALTERNATIVE: Direct CSV Export (even simpler!)
# ============================================================================
# Uncomment this to use the direct CSV export method:
#
# print("\n[Alternative] Direct CSV export...")
# ms2library.analog_search_store_in_csv(
#     query_spectra=processed,
#     results_csv_file_location=str(OUTPUT_CSV),
#     settings=settings
# )
# print(f"Results saved directly to: {OUTPUT_CSV}")
