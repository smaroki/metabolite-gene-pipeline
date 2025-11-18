import pandas as pd
import os

# ==============================
# File paths
# ==============================
INPUT_FILE = r"C:\Users\dell\OneDrive\Desktop\Raw files\results\Literature metabolites.xlsx"
OUTPUT_FILE = os.path.join(os.path.dirname(INPUT_FILE), "Literature_metabolites_cleaned.csv")

# ==============================
# Load the Excel file
# ==============================
df = pd.read_excel(INPUT_FILE)

# ==============================
# Check if required column exists
# ==============================
if "Metabolite Name" not in df.columns:
    raise ValueError("❌ The Excel file must contain a column named 'Metabolite Name'")

# ==============================
# Clean the metabolite names
# ==============================
df["Cleaned Name"] = (
    df["Metabolite Name"]
    .astype(str)
    .str.strip()                               # remove leading/trailing spaces
    .str.replace(r"\(.*?\)", "", regex=True)   # remove text inside parentheses
    .str.replace(r"\[.*?\]", "", regex=True)   # remove text inside brackets
    .str.replace("-", " ")                     # replace dashes with spaces
    .str.replace("_", " ")                     # replace underscores with spaces
    .str.replace(r"[^a-zA-Z0-9\s]", "", regex=True)  # remove special characters
    .str.replace(r"\s+", " ", regex=True)      # collapse multiple spaces into one
    .str.strip()                               # trim again
)

# ==============================
# Save the cleaned version
# ==============================
df.to_csv(OUTPUT_FILE, index=False, encoding="utf-8")

print(f"✅ Cleaned file saved as CSV:\n{OUTPUT_FILE}")
print("🔹 Added new column: 'Cleaned Name'")
