import pandas as pd
import requests

# Load your cleaned metabolites CSV
input_path = r"C:\Users\dell\OneDrive\Desktop\Raw files\results\Literature_metabolites_cleaned.csv"
df = pd.read_csv(input_path)

# Ensure 'Cleaned' column is clean
df["Cleaned Name"] = df["Cleaned Name"].astype(str).str.strip()

# Function to fetch SMILES from PubChem
def fetch_smiles(name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/CanonicalSMILES,ConnectivitySMILES/JSON"
    print(f"🔍 Searching: {name}")
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            props = data.get("PropertyTable", {}).get("Properties", [])
            if props:
                smiles = props[0].get("CanonicalSMILES") or props[0].get("ConnectivitySMILES")
                if smiles:
                    print(f"✅ Found: {smiles}")
                    return smiles
                else:
                    print(f"⚠️ No SMILES in response for: {name}")
                    return None
            else:
                print(f"⚠️ No Properties returned for: {name}")
                return None
        else:
            print(f"❌ HTTP error {response.status_code} for: {name}")
            return None
    except Exception as e:
        print(f"❌ Error for {name}: {e}")
        return None

# Apply to all cleaned names
df["SMILES"] = df["Cleaned Name"].apply(fetch_smiles)

# Save output
output_path = r"C:\Users\dell\OneDrive\Desktop\Raw files\results\Metabolites_with_SMILES.csv"
df.to_csv(output_path, index=False)

print(f"\n✅ All done! SMILES saved to: {output_path}")

