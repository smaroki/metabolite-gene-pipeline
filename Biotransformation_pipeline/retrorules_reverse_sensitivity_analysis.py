import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

Validated_SMILES_path = r"C:\Users\dell\Desktop\Final Manuscript Data\Original pipeline output\Validated_SMILES.csv"
retrorules_path = r"C:\Users\dell\OneDrive\Desktop\Raw files\results\retrorules_rr02_rp2_hs (1)\retrorules_rr02_rp2_hs\retrorules_rr02_rp2_flat_all.csv"

output_folder = r"C:\Users\dell\Desktop\Final Manuscript Data\Original pipeline output\precursor_output"
os.makedirs(output_folder, exist_ok=True)

# Define score thresholds for sensitivity analysis
SCORE_THRESHOLDS = [0.0, 0.3, 0.5, 0.7, 0.9]

# Load retrorules
print("Loading retrorules...")
retrorules_df = pd.read_csv(retrorules_path)
print(f"Total retrorules loaded: {len(retrorules_df)}")

# Filter by reaction direction
retrorules_df = retrorules_df[retrorules_df['Reaction direction'] == -1]
print(f"Retrorules with direction -1: {len(retrorules_df)}")

# USE 'Score normalized' - it's already 0-1 range
if 'Score normalized' not in retrorules_df.columns:
    raise ValueError("'Score normalized' column not found in retrorules file!")

print(f"\nUsing 'Score normalized' as confidence metric")

# Keep essential columns and rename
retrorules_df = retrorules_df[['Rule ID', 'Rule', 'EC number', 'Score normalized']].dropna()
retrorules_df.columns = ['rule_id', 'reaction_smarts', 'ec_number', 'score']

print(f"Valid rules after cleanup: {len(retrorules_df)}")

print(f"\nScore distribution:")
print(retrorules_df['score'].describe())
print(f"\nRules per threshold:")
for threshold in SCORE_THRESHOLDS:
    count = len(retrorules_df[retrorules_df['score'] >= threshold])
    print(f"  Score >= {threshold}: {count} rules")

# Load SMILES
print("\nLoading validated SMILES...")
smiles_df = pd.read_csv(Validated_SMILES_path)
print(f"Total SMILES loaded: {len(smiles_df)}")

if 'validated_smiles' not in smiles_df.columns:
    raise ValueError("The file must contain 'validated_smiles' column.")

# Function to run predictions at a given threshold
def predict_with_threshold(rules_df, smiles_df, threshold):
    """Run retrosynthesis with rules above the given score threshold"""
    filtered_rules = rules_df[rules_df['score'] >= threshold].copy()
    
    results = []
    valid_metabolites = 0
    invalid_metabolites = 0
    
    for idx, row in tqdm(smiles_df.iterrows(), total=len(smiles_df), 
                         desc=f"Threshold {threshold}"):
        validated_smiles = row['validated_smiles']
        
        if pd.isna(validated_smiles) or validated_smiles == '':
            invalid_metabolites += 1
            continue
        
        metabolite_id = f"Metabolite_{idx}"
        
        try:
            product_mol = Chem.MolFromSmiles(validated_smiles)
            if product_mol is None:
                invalid_metabolites += 1
                continue
            valid_metabolites += 1
        except:
            invalid_metabolites += 1
            continue
        
        for _, rule in filtered_rules.iterrows():
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
                            'Predicted Precursor Name': '',
                            'Predicted Precursor SMILES': precursor_smiles,
                            'Reaction Type': 'reverse',
                            'Rule ID': rule_id,
                            'EC Number': ec_number,
                            'Reaction SMARTS': smarts,
                            'Rule Score': score,
                            'Score Threshold': threshold
                        })
            except Exception:
                continue
    
    return results, valid_metabolites, invalid_metabolites

# Run sensitivity analysis
print("\n" + "="*60)
print("SENSITIVITY ANALYSIS: Testing multiple score thresholds")
print("="*60)

all_results = []
summary_stats = []

for threshold in SCORE_THRESHOLDS:
    print(f"\n{'='*60}")
    print(f"Running with score threshold >= {threshold}")
    print(f"{'='*60}")
    
    results, valid, invalid = predict_with_threshold(
        retrorules_df, smiles_df, threshold
    )
    
    all_results.extend(results)
    
    stats = {
        'Threshold': threshold,
        'Rules Used': len(retrorules_df[retrorules_df['score'] >= threshold]),
        'Valid Metabolites': valid,
        'Invalid Metabolites': invalid,
        'Total Precursors': len(results)
    }
    summary_stats.append(stats)
    
    print(f"Rules used: {stats['Rules Used']}")
    print(f"Valid metabolites: {stats['Valid Metabolites']}")
    print(f"Precursors predicted: {stats['Total Precursors']}")
    
    # Save individual threshold results
    threshold_file = os.path.join(
        output_folder, 
        f"Precursors_threshold_{threshold:.1f}.csv"
    )
    pd.DataFrame(results).to_csv(threshold_file, index=False)
    print(f"Saved: {threshold_file}")

# Save combined results
combined_file = os.path.join(output_folder, "Precursors_all_thresholds.csv")
pd.DataFrame(all_results).to_csv(combined_file, index=False)

# Save summary statistics
summary_file = os.path.join(output_folder, "Sensitivity_Analysis_Summary.csv")
summary_df = pd.DataFrame(summary_stats)
summary_df.to_csv(summary_file, index=False)

# Print final summary
print(f"\n{'='*60}")
print("SENSITIVITY ANALYSIS COMPLETE")
print(f"{'='*60}")
print("\nSummary by threshold:")
print(summary_df.to_string(index=False))
print(f"\nAll results saved to: {output_folder}")
print(f"{'='*60}")
