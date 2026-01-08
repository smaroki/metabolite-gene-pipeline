
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

# ============================================================================
# Configuration
# ============================================================================
Validated_SMILES_path = r"C:\Users\dell\Desktop\Final Manuscript Data\Original pipeline output\Validated_SMILES.csv"
retrorules_path = r"C:\Users\dell\Desktop\Raw files\results\retrorules_rr02_rp2_hs (1)\retrorules_rr02_rp2_hs\retrorules_rr02_rp2_flat_all.csv"

output_folder = r"C:\Users\dell\Desktop\Final Manuscript Data\Original pipeline output\precursor_output3"

# Confidence score threshold (based on sensitivity analysis - see Methods)
CONFIDENCE_THRESHOLD = 0.5

# ============================================================================
# Setup
# ============================================================================
os.makedirs(output_folder, exist_ok=True)
output_file = os.path.join(output_folder, "Predicted_Precursors.csv")

# ============================================================================
# Load and Filter RetroRules
# ============================================================================
print("Loading RetroRules...")
retrorules_df = pd.read_csv(retrorules_path)
print(f"Total rules loaded: {len(retrorules_df)}")

# Filter for reverse reactions (Reaction direction = -1)
retrorules_df = retrorules_df[retrorules_df['Reaction direction'] == -1]
print(f"Reverse-direction rules: {len(retrorules_df)}")

# Select and rename columns
if 'Score normalized' not in retrorules_df.columns:
    raise ValueError("'Score normalized' column not found in RetroRules file!")

retrorules_df = retrorules_df[['Rule ID', 'Rule', 'EC number', 'Score normalized']].dropna()
retrorules_df.columns = ['rule_id', 'reaction_smarts', 'ec_number', 'score']

# Apply confidence score threshold
retrorules_df = retrorules_df[retrorules_df['score'] >= CONFIDENCE_THRESHOLD]
print(f"Rules after score filtering (>= {CONFIDENCE_THRESHOLD}): {len(retrorules_df)}")

if len(retrorules_df) == 0:
    raise ValueError(f"No rules remain after applying threshold {CONFIDENCE_THRESHOLD}. Check your data.")

print(f"\nScore distribution of selected rules:")
print(retrorules_df['score'].describe())

# ============================================================================
# Load Validated SMILES
# ============================================================================
print("\nLoading validated SMILES...")
smiles_df = pd.read_csv(Validated_SMILES_path)
print(f"Total metabolites loaded: {len(smiles_df)}")

if 'validated_smiles' not in smiles_df.columns:
    raise ValueError("Input file must contain 'validated_smiles' column.")

# ============================================================================
# Predict Precursors
# ============================================================================
print("\nPredicting biochemical precursors...")
results = []
valid_metabolites = 0
invalid_metabolites = 0

for idx, row in tqdm(smiles_df.iterrows(), total=len(smiles_df), desc="Processing metabolites"):
    validated_smiles = row['validated_smiles']
    
    # Skip invalid SMILES
    if pd.isna(validated_smiles) or validated_smiles == '':
        invalid_metabolites += 1
        continue
    
    metabolite_id = f"Metabolite_{idx}"
    
    # Parse SMILES
    try:
        product_mol = Chem.MolFromSmiles(validated_smiles)
        if product_mol is None:
            invalid_metabolites += 1
            continue
        valid_metabolites += 1
    except Exception:
        invalid_metabolites += 1
        continue
    
    # Apply all rules
    for _, rule in retrorules_df.iterrows():
        smarts = rule['reaction_smarts']
        rule_id = rule['rule_id']
        ec_number = rule['ec_number']
        score = rule['score']
        
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((product_mol,))
            
            for prod in products:
                if prod and len(prod) > 0:
                    precursor_smiles = Chem.MolToSmiles(prod[0])
                    results.append({
                        'Original Metabolite ID': metabolite_id,
                        'Original Metabolite SMILES': validated_smiles,
                        'Predicted Precursor SMILES': precursor_smiles,
                        'Rule ID': rule_id,
                        'EC Number': ec_number,
                        'Reaction SMARTS': smarts,
                        'Rule Confidence Score': score
                    })
        except Exception:
            continue

# ============================================================================
# Save Results
# ============================================================================
results_df = pd.DataFrame(results)
results_df.to_csv(output_file, index=False)

print(f"\n{'='*70}")
print("REVERSE BIOTRANSFORMATION COMPLETE")
print(f"{'='*70}")
print(f"Confidence threshold used: {CONFIDENCE_THRESHOLD}")
print(f"Rules applied: {len(retrorules_df)}")
print(f"Valid metabolites processed: {valid_metabolites}")
print(f"Invalid/skipped metabolites: {invalid_metabolites}")
print(f"Total precursors predicted: {len(results_df)}")
print(f"\nOutput saved to:")
print(f"  {output_file}")
print(f"{'='*70}")
