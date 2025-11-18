import os
import pandas as pd
from rdkit import Chem
from pathlib import Path

# ============================================================================
# CONFIGURATION - Match your MS2Query output
# ============================================================================
INPUT_FILE = Path(r"C:\Users\dell\Desktop\Final Manuscript Data\Original pipeline output\FM2B5_predictions.csv")
OUTPUT_FILE = Path(r"C:\Users\dell\Desktop\Final Manuscript Data\Original pipeline output\Validated_SMILES.csv")

print("=" * 80)
print("RDKit SMILES Validation")
print("=" * 80)

def validate_smiles(smiles: str) -> str:
    """
    Check if SMILES is valid using RDKit. 
    Return canonical SMILES if valid, else None.
    """
    if pd.isna(smiles) or smiles == '' or not isinstance(smiles, str):
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)  # standardize form
    except Exception:
        return None

def get_molecular_properties(smiles: str) -> dict:
    """
    Calculate molecular properties from SMILES.
    Returns dict with molecular weight, formula, etc.
    """
    if pd.isna(smiles) or smiles == '':
        return {
            'molecular_weight': None,
            'num_atoms': None,
            'num_heavy_atoms': None,
            'num_rings': None
        }
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                'molecular_weight': None,
                'num_atoms': None,
                'num_heavy_atoms': None,
                'num_rings': None
            }
        
        from rdkit.Chem import Descriptors
        return {
            'molecular_weight': round(Descriptors.MolWt(mol), 2),
            'num_atoms': mol.GetNumAtoms(),
            'num_heavy_atoms': mol.GetNumHeavyAtoms(),
            'num_rings': Descriptors.RingCount(mol)
        }
    except Exception:
        return {
            'molecular_weight': None,
            'num_atoms': None,
            'num_heavy_atoms': None,
            'num_rings': None
        }

def main():
    # Create output file
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"\nInput:  {INPUT_FILE}")
    print(f"Output: {OUTPUT_FILE}")
    
    # Load input file
    print(f"\n[1/4] Loading predictions...")
    df = pd.read_csv(INPUT_FILE)
    print(f"Loaded {len(df)} rows")
    
    # Check for SMILES column (try different names)
    smiles_col = None
    for col in ['smiles', 'SMILES', 'Smiles', 'analog_smiles']:
        if col in df.columns:
            smiles_col = col
            break
    
    if smiles_col is None:
        raise ValueError(f"No SMILES column found! Available columns: {list(df.columns)}")
    
    print(f"Found SMILES column: '{smiles_col}'")
    
    # Count non-empty SMILES
    non_empty = df[smiles_col].notna() & (df[smiles_col] != '')
    print(f"   Non-empty SMILES: {non_empty.sum()}/{len(df)}")
    
    # Validate and standardize
    print(f"\n[2/4] Validating SMILES with RDKit...")
    df['validated_smiles'] = df[smiles_col].apply(validate_smiles)
    df['is_valid'] = df['validated_smiles'].notnull()
    
    valid_count = df['is_valid'].sum()
    print(f"Valid SMILES: {valid_count}/{len(df)} ({100*valid_count/len(df):.1f}%)")
    
    # Calculate molecular properties for valid SMILES
    print(f"\n[3/4] Calculating molecular properties...")
    props = df['validated_smiles'].apply(get_molecular_properties)
    df['molecular_weight'] = props.apply(lambda x: x['molecular_weight'])
    df['num_atoms'] = props.apply(lambda x: x['num_atoms'])
    df['num_heavy_atoms'] = props.apply(lambda x: x['num_heavy_atoms'])
    df['num_rings'] = props.apply(lambda x: x['num_rings'])
    
    print(f"Properties calculated")
    
    # Reorder columns for better readability
    base_cols = ['spectrum_id', 'compound_name', 'analog_score', 'rank']
    smiles_cols = [smiles_col, 'validated_smiles', 'is_valid']
    prop_cols = ['molecular_weight', 'num_atoms', 'num_heavy_atoms', 'num_rings']
    other_cols = [col for col in df.columns if col not in base_cols + smiles_cols + prop_cols]
    
    # Reorder
    available_base = [col for col in base_cols if col in df.columns]
    df = df[available_base + smiles_cols + prop_cols + other_cols]
    
    # Save results
    print(f"\n[4/4] Saving results...")
    output_path = OUTPUT_FILE
    df.to_csv(output_path, index=False)
    
    # ========================================================================
    # STATISTICS
    # ========================================================================
    print(f"\n{'=' * 80}")
    print(f"VALIDATION COMPLETE")
    print(f"{'=' * 80}")
    
    print(f"\nResults saved to: {output_path}")
    print(f"Total rows: {len(df)}")
    print(f"Valid SMILES: {valid_count} ({100*valid_count/len(df):.1f}%)")
    print(f"Invalid SMILES: {len(df) - valid_count} ({100*(len(df)-valid_count)/len(df):.1f}%)")
    
    # Show invalid cases
    invalid_df = df[~df['is_valid'] & df[smiles_col].notna()]
    if len(invalid_df) > 0:
        print(f"\nWARNING: Invalid SMILES examples (first 5):")
        for i, (_, row) in enumerate(invalid_df.head(5).iterrows(), 1):
            print(f"   {i}. {row['compound_name']}: '{row[smiles_col][:60]}'")
    
    # Molecular weight statistics for valid compounds
    valid_df = df[df['is_valid']]
    if len(valid_df) > 0:
        print(f"\nMolecular Weight Statistics (valid compounds):")
        print(f"   Mean:   {valid_df['molecular_weight'].mean():.2f}")
        print(f"   Median: {valid_df['molecular_weight'].median():.2f}")
        print(f"   Min:    {valid_df['molecular_weight'].min():.2f}")
        print(f"   Max:    {valid_df['molecular_weight'].max():.2f}")
        
        # Show top 5 valid matches
        print(f"\nTOP 5 VALID MATCHES (by score):")
        top_valid = valid_df.nlargest(5, 'analog_score')
        for i, (_, row) in enumerate(top_valid.iterrows(), 1):
            print(f"   {i}. {row['compound_name'][:40]}")
            print(f"      Score: {row['analog_score']:.4f} | MW: {row['molecular_weight']:.1f}")
            print(f"      SMILES: {row['validated_smiles'][:60]}...")
    
    # Sample output
    print(f"\nSample Validated Results (first 10 rows):")
    display_cols = ['spectrum_id', 'compound_name', 'analog_score', 'is_valid', 'molecular_weight']
    available_display = [col for col in display_cols if col in df.columns]
    print(df[available_display].head(10).to_string(index=False))
    
    print(f"\nDone!")

if __name__ == "__main__":
    main()
