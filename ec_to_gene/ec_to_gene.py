import pandas as pd
import requests
import json
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import tempfile

input_file = r"C:\Users\dell\Desktop\Final Manuscript Data\Original pipeline output\precursor_output\Predicted_Precursors_with_RxnType.csv"
output_file = r"C:\Users\dell\Desktop\Final Manuscript Data\Original pipeline output\precursor_output\EC_to_Genes_Uniprot_Annotated.csv"
cache_file = os.path.join(tempfile.gettempdir(), "uniprot_ec_cache.json")

# Load data
df = pd.read_csv(input_file)
print(f"Loaded input file with {len(df)} rows")
print(f"Columns: {df.columns.tolist()}")

if 'EC Number' not in df.columns:
    raise KeyError("Ensure 'EC Number' column exists.")

BASE_URL = "https://rest.uniprot.org/uniprotkb/search"

# Load cache
if os.path.exists(cache_file):
    with open(cache_file, "r") as f:
        cache = json.load(f)
else:
    cache = {}

def generate_ec_queries(ec):
    ec_list = [e.strip() for e in str(ec).split(';')]
    all_queries = []
    for e in ec_list:
        parts = e.split('.')
        if len(parts) == 4:
            all_queries.append(e)
            all_queries.append('.'.join(parts[:3]) + ".-")
            all_queries.append('.'.join(parts[:2]) + ".-.-")
            all_queries.append(parts[0] + ".-.-.-")
        else:
            all_queries.append(e)
    return list(dict.fromkeys(all_queries))

def query_uniprot(ec_query, reviewed=True):
    q = f"ec:{ec_query}"
    if reviewed:
        q += " AND reviewed:true"
    params = {
        "query": q,
        "format": "tsv",
        "fields": "gene_names,protein_name",
        "size": 1
    }
    try:
        r = requests.get(BASE_URL, params=params, timeout=10)
        if r.status_code == 200:
            lines = r.text.strip().split('\n')
            if len(lines) > 1:
                parts = lines[1].split('\t')
                if len(parts) >= 2:
                    return parts[0], parts[1]
    except Exception as e:
        print(f"WARNING: Error querying {ec_query}: {e}")
    return None, None

def process_ec_number(ec_number):
    """Process a single EC number with fallback logic"""
    if ec_number in cache:
        return ec_number, cache[ec_number]
    
    gene = protein = None
    match_type = "No match"
    queries = generate_ec_queries(ec_number)
    
    # Try reviewed first
    for q in queries:
        gene, protein = query_uniprot(q, reviewed=True)
        if gene or protein:
            match_type = "Perfect match"
            break
    
    # Fallback to unreviewed
    if not gene and not protein:
        for q in queries:
            gene, protein = query_uniprot(q, reviewed=False)
            if gene or protein:
                match_type = "Close match"
                break
    
    if not gene and not protein:
        gene = protein = ""
    
    result = (gene, protein, match_type)
    cache[ec_number] = result
    return ec_number, result

# Get unique EC numbers to process
unique_ecs = df['EC Number'].astype(str).str.strip().unique()
ecs_to_process = [ec for ec in unique_ecs if ec not in cache]

print(f"\nTotal EC numbers: {len(unique_ecs)}")
print(f"Already cached: {len(unique_ecs) - len(ecs_to_process)}")
print(f"Need to query: {len(ecs_to_process)}")

# Process in parallel with progress tracking
if ecs_to_process:
    print("\nStarting parallel queries...")
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = {executor.submit(process_ec_number, ec): ec for ec in ecs_to_process}
        
        completed = 0
        for future in as_completed(futures):
            ec_number, result = future.result()
            completed += 1
            if completed % 10 == 0 or completed == len(ecs_to_process):
                print(f"   Progress: {completed}/{len(ecs_to_process)} ({100*completed//len(ecs_to_process)}%)")
    
    # Save cache after processing
    with open(cache_file, "w") as f:
        json.dump(cache, f, indent=2)
    print(f"Cache saved to: {cache_file}")

# Add gene/protein columns to the original dataframe
print("\nAdding gene and protein annotations...")
gene_names = []
protein_names = []
match_types = []

for _, row in df.iterrows():
    ec_number = str(row['EC Number']).strip()
    gene, protein, match_type = cache.get(ec_number, ("", "", "No match"))
    gene_names.append(gene)
    protein_names.append(protein)
    match_types.append(match_type)

# Add new columns to the dataframe
df['Gene Name'] = gene_names
df['Protein Name'] = protein_names
df['Match Type'] = match_types

# Save the result
df.to_csv(output_file, index=False)

print(f"\nDone! Output saved to:\n{output_file}")
print(f"Total results: {len(df)}")
print(f"   Perfect matches: {sum(1 for m in match_types if m == 'Perfect match')}")
print(f"   Close matches: {sum(1 for m in match_types if m == 'Close match')}")
print(f"   No matches: {sum(1 for m in match_types if m == 'No match')}")
