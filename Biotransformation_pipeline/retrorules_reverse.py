import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

# ============================================================================
# Configuration
# ============================================================================
Validated_SMILES_path = r"D:\OneDrive\Desktop\SMAROKI\ValidatedSMILES.csv"
retrorules_path = r"E:\Downloads\retrorules_rr02_rp2_hs (1)\retrorules_rr02_rp2_hs\retrorules_rr02_rp2_flat_all.csv"

# Output location
output_folder = r"D:\OneDrive\Desktop\SMAROKI"
output_file = os.path.join(output_folder, "Predicted_Precursors.csv")

# Confidence score threshold
CONFIDENCE_THRESHOLD = 0.5

# ============================================================================
# Setup
# ============================================================================
os.makedirs(output_folder, exist_ok=True)

# ============================================================================
# Load and Filter RetroRules
# ============================================================================
print("Loading RetroRules...")
retrorules_df = pd.read_csv(retrorules_path)
print(f"Total rules loaded: {len(retrorules_df)}")

# Filter for reverse reactions
retrorules_df = retrorules_df[retrorules_df['Reaction direction'] == -1]
print(f"Reverse-direction rules: {len(retrorules_df)}")

# Check required RetroRules columns
required_retrorules_cols = ['Rule ID', 'Rule', 'EC number', 'Score normalized']
for col in required_retrorules_cols:
    if col not in retrorules_df.columns:
        raise ValueError(f"Required column '{col}' not found in RetroRules file.")

# Keep needed columns
retrorules_df = retrorules_df[['Rule ID', 'Rule', 'EC number', 'Score normalized']].dropna()
retrorules_df.columns = ['rule_id', 'reaction_smarts', 'ec_number', 'score']

# Apply confidence filter
retrorules_df = retrorules_df[retrorules_df['score'] >= CONFIDENCE_THRESHOLD]
print(f"Rules after score filtering (>= {CONFIDENCE_THRESHOLD}): {len(retrorules_df)}")

if len(retrorules_df) == 0:
    raise ValueError(f"No rules remain after applying threshold {CONFIDENCE_THRESHOLD}.")

print("\nScore distribution of selected rules:")
print(retrorules_df['score'].describe())

# ============================================================================
# Load Input Metabolites
# ============================================================================
print("\nLoading validated SMILES...")
smiles_df = pd.read_csv(Validated_SMILES_path)
print(f"Total metabolites loaded: {len(smiles_df)}")
print("Columns found:")
print(smiles_df.columns.tolist())

# Required columns in your metabolite file
required_input_cols = ['validated_smiles', 'Metabolite Name']
for col in required_input_cols:
    if col not in smiles_df.columns:
        raise ValueError(f"Required column '{col}' not found in input file.")

# ============================================================================
# Predict Precursors
# ============================================================================
print("\nPredicting biochemical precursors...")
results = []
valid_metabolites = 0
invalid_metabolites = 0

for idx, row in tqdm(smiles_df.iterrows(), total=len(smiles_df), desc="Processing metabolites"):
    validated_smiles = row['validated_smiles']
    metabolite_name = row['Metabolite Name']
    hmdb_id = row.get('HMDB IDs', None)

    # Skip missing validated SMILES
    if pd.isna(validated_smiles) or str(validated_smiles).strip() == '':
        invalid_metabolites += 1
        continue

    # If metabolite name is missing, give fallback name
    if pd.isna(metabolite_name) or str(metabolite_name).strip() == '':
        metabolite_name = f"Metabolite_{idx}"

    # Parse validated_smiles
    try:
        product_mol = Chem.MolFromSmiles(str(validated_smiles))
        if product_mol is None:
            invalid_metabolites += 1
            continue
        valid_metabolites += 1
    except Exception:
        invalid_metabolites += 1
        continue

    # Apply RetroRules
    for _, rule in retrorules_df.iterrows():
        smarts = rule['reaction_smarts']
        rule_id = rule['rule_id']
        ec_number = rule['ec_number']
        score = rule['score']

        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            if rxn is None:
                continue

            products = rxn.RunReactants((product_mol,))

            for prod in products:
                if prod and len(prod) > 0:
                    try:
                        precursor_smiles = Chem.MolToSmiles(prod[0])
                        results.append({
                            'Metabolite Name': metabolite_name,
                            'HMDB IDs': hmdb_id,
                            'Validated Metabolite SMILES': validated_smiles,
                            'Predicted Precursor SMILES': precursor_smiles,
                            'Rule ID': rule_id,
                            'EC Number': ec_number,
                            'Reaction SMARTS': smarts,
                            'Rule Confidence Score': score
                        })
                    except Exception:
                        continue

        except Exception:
            continue

# ============================================================================
# Save Results
# ============================================================================
results_df = pd.DataFrame(results)

# Optional: reorder columns
priority_cols = [
    'Metabolite Name',
    'HMDB IDs',
    'Validated Metabolite SMILES',
    'Predicted Precursor SMILES',
    'Rule ID',
    'EC Number',
    'Reaction SMARTS',
    'Rule Confidence Score'
]
existing_priority = [c for c in priority_cols if c in results_df.columns]
other_cols = [c for c in results_df.columns if c not in priority_cols]
results_df = results_df[existing_priority + other_cols]

results_df.to_csv(output_file, index=False)

print(f"\n{'='*70}")
print("REVERSE BIOTRANSFORMATION COMPLETE")
print(f"{'='*70}")
print(f"Confidence threshold used : {CONFIDENCE_THRESHOLD}")
print(f"Rules applied             : {len(retrorules_df)}")
print(f"Valid metabolites         : {valid_metabolites}")
print(f"Invalid/skipped           : {invalid_metabolites}")
print(f"Total precursors predicted: {len(results_df)}")
print("\nOutput saved to:")
print(f"  {output_file}")
print(f"{'='*70}")

