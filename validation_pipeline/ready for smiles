import pandas as pd
import os
import re

# ==============================
# File paths
# ==============================
INPUT_FILE = r"C:\Users\dell\OneDrive\Desktop\Raw files\results\Literature_metabolites_cleaned.csv"
OUTPUT_FILE = os.path.join(os.path.dirname(INPUT_FILE), "Literature_metabolites_ready_for_SMILES.csv")

# ==============================
# Load the cleaned CSV
# ==============================
df = pd.read_csv(INPUT_FILE)
df.columns = df.columns.str.strip()  # Remove hidden spaces

if "Cleaned Name" not in df.columns:
    raise ValueError("❌ The CSV must contain a column named 'Cleaned Name'")

# ==============================
# Function to split multiple metabolites in one cell
# ==============================
def split_metabolites(name):
    if pd.isna(name):
        return []
    # Split by comma or multiple spaces
    parts = re.split(r'[,\s]{2,}|,', str(name))
    parts = [p.strip() for p in parts if p.strip()]
    return parts

# ==============================
# Function to filter vague entries
# ==============================
def is_valid_metabolite(name):
    vague_keywords = [
        "peptide", "metabolite", "metabolites", "fatty acid", "hydroxylation",
        "glycerophosphocholines", "fibrinogen"
    ]
    name_lower = name.lower()
    for kw in vague_keywords:
        if kw in name_lower:
            return False
    return True

# ==============================
# Generate list of individual metabolites
# ==============================
all_metabolites = []
for entry in df["Cleaned Name"]:
    split_names = split_metabolites(entry)
    valid_names = [n for n in split_names if is_valid_metabolite(n)]
    all_metabolites.extend(valid_names)

# Remove duplicates and sort
all_metabolites = sorted(list(set(all_metabolites)))

# ==============================
# Save as CSV
# ==============================
df_ready = pd.DataFrame({"Metabolite Name": all_metabolites})
df_ready.to_csv(OUTPUT_FILE, index=False, encoding="utf-8")

print(f"✅ Preprocessed metabolites saved to:\n{OUTPUT_FILE}")
print(f"🔹 Total individual metabolites ready for PubChem: {len(all_metabolites)}")
