import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

Validated_SMILES_path = r"C:\Users\dell\Desktop\Final Manuscript Data\Original pipeline output\Validated_SMILES.csv"
retrorules_path = r"C:\Users\dell\OneDrive\Desktop\Raw files\results\retrorules_rr02_rp2_hs (1)\retrorules_rr02_rp2_hs\retrorules_rr02_rp2_flat_all.csv"

output_folder = r"C:\Users\dell\Desktop\Final Manuscript Data\Original pipeline output\precursor_output"
os.makedirs(output_folder, exist_ok=True)
output_file = os.path.join(output_folder, "Predicted_Precursors_with_RxnType.csv")

# Load retrorules
print("Loading retrorules...")
retrorules_df = pd.read_csv(retrorules_path)
print(f"Total retrorules loaded: {len(retrorules_df)}")

retrorules_df = retrorules_df[retrorules_df['Reaction direction'] == -1]
print(f"Retrorules with direction -1: {len(retrorules_df)}")

retrorules_df = retrorules_df[['Rule ID', 'Rule', 'EC number']].dropna()
retrorules_df.columns = ['rule_id', 'reaction_smarts', 'ec_number']
print(f"Valid retrorules after cleanup: {len(retrorules_df)}")

# Load SMILES
print("\nLoading validated SMILES...")
smiles_df = pd.read_csv(Validated_SMILES_path)
print(f"Total SMILES loaded: {len(smiles_df)}")
print(f"Columns in file: {smiles_df.columns.tolist()}")

if 'validated_smiles' not in smiles_df.columns:
    raise ValueError("The file must contain 'validated_smiles' column.")

results = []
valid_metabolites = 0
invalid_metabolites = 0

print("\nProcessing metabolites...")
for idx, row in tqdm(smiles_df.iterrows(), total=len(smiles_df)):
    validated_smiles = row['validated_smiles']
    
    # Skip if SMILES is NaN or empty
    if pd.isna(validated_smiles) or validated_smiles == '':
        invalid_metabolites += 1
        continue
    
    # Use index as identifier
    metabolite_id = f"Metabolite_{idx}"
    
    try:
        # FIXED: Use validated_smiles instead of metabolite_smiles
        product_mol = Chem.MolFromSmiles(validated_smiles)
        if product_mol is None:
            invalid_metabolites += 1
            continue
        valid_metabolites += 1
    except:
        invalid_metabolites += 1
        continue
    
    for _, rule in retrorules_df.iterrows():
        smarts = rule['reaction_smarts']
        rule_id = rule['rule_id']
        ec_number = rule['ec_number']
        
        try:
            rxn = AllChem.ReactionFromSmarts(smarts)
            products = rxn.RunReactants((product_mol,))
            
            for prod in products:
                if prod and len(prod) > 0:
                    precursor_smiles = Chem.MolToSmiles(prod[0])
                    results.append({
                        'Original Metabolite ID': metabolite_id,
                        'Original Metabolite SMILES': validated_smiles,  # FIXED
                        'Predicted Precursor Name': '',
                        'Predicted Precursor SMILES': precursor_smiles,
                        'Reaction Type': 'reverse',
                        'Rule ID': rule_id,
                        'EC Number': ec_number,
                        'Reaction SMARTS': smarts
                    })
        except Exception:
            continue

# Save results
results_df = pd.DataFrame(results)
results_df.to_csv(output_file, index=False)

print(f"\n{'='*60}")
print(f"✅ Prediction complete!")
print(f"{'='*60}")
print(f"Valid metabolites processed: {valid_metabolites}")
print(f"Invalid/skipped metabolites: {invalid_metabolites}")
print(f"Total precursors predicted: {len(results_df)}")
print(f"Output saved to:\n{output_file}")
print(f"{'='*60}")
