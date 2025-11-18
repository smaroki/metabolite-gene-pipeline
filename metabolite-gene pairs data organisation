import pandas as pd
import requests
import time
from collections import defaultdict
import os
import re

# File paths
GENE_FILE = r"C:\Users\dell\OneDrive\Desktop\Untargeted metabolomics\edges - edges.csv"
METABOLITE_FILE = r"C:\Users\dell\OneDrive\Desktop\Untargeted metabolomics\nodes - nodes.csv"
OUTPUT_DIR = r"C:\Users\dell\OneDrive\Desktop\Untargeted metabolomics\kegg_validation_output"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

def clean_metabolite_name(metabolite_name):
    """Clean and simplify metabolite name for KEGG search"""
    # Convert to lowercase
    cleaned = metabolite_name.lower().strip()
    
    # Remove common suffixes and annotations
    patterns_to_remove = [
        r'\([^)]*\)',  # Remove anything in parentheses
        r'\[[^\]]*\]',  # Remove anything in brackets
        r'\s*\d+[:\-]\d+.*'  # Remove lipid notation like 18:1, 16:0
    ]
    
    for pattern in patterns_to_remove:
        cleaned = re.sub(pattern, '', cleaned)
    
    # Remove special characters except hyphens and spaces
    cleaned = re.sub(r'[^a-z0-9\s\-]', ' ', cleaned)
    
    # Remove extra spaces
    cleaned = ' '.join(cleaned.split())
    
    return cleaned.strip()

def query_kegg_compound(metabolite_name):
    """Search for a metabolite in KEGG and return compound IDs"""
    try:
        # Try original name first
        searches = [metabolite_name]
        
        # Add cleaned version
        cleaned = clean_metabolite_name(metabolite_name)
        if cleaned != metabolite_name.lower():
            searches.append(cleaned)
        
        # Try first word only (often the base compound name)
        first_word = cleaned.split()[0] if cleaned.split() else cleaned
        if len(first_word) > 3 and first_word not in searches:
            searches.append(first_word)
        
        compound_ids = []
        
        for search_term in searches:
            if not search_term or len(search_term) < 3:
                continue
                
            url = f"http://rest.kegg.jp/find/compound/{search_term}"
            response = requests.get(url, timeout=10)
            time.sleep(0.3)
            
            if response.status_code == 200 and response.text.strip():
                for line in response.text.strip().split('\n'):
                    compound_id = line.split('\t')[0]
                    compound_ids.append(compound_id)
                
                if compound_ids:
                    print(f"  Search term used: '{search_term}'")
                    return compound_ids
        
        return []
    except Exception as e:
        print(f"Error querying metabolite {metabolite_name}: {e}")
        return []

def get_genes_for_compound(compound_id):
    """Get all genes associated with a compound through reactions and enzymes"""
    try:
        genes_found = set()
        
        # Get reactions for this compound
        url = f"http://rest.kegg.jp/link/reaction/{compound_id}"
        response = requests.get(url, timeout=10)
        time.sleep(0.3)
        
        if response.status_code == 200 and response.text.strip():
            reactions = [line.split('\t')[1] for line in response.text.strip().split('\n')]
            
            # For each reaction, get enzymes
            for reaction in reactions[:10]:  # Limit to avoid too many queries
                try:
                    url = f"http://rest.kegg.jp/link/enzyme/{reaction}"
                    response = requests.get(url, timeout=10)
                    time.sleep(0.3)
                    
                    if response.status_code == 200 and response.text.strip():
                        enzymes = [line.split('\t')[1] for line in response.text.strip().split('\n')]
                        
                        # For each enzyme, get human genes
                        for enzyme in enzymes:
                            try:
                                url = f"http://rest.kegg.jp/link/hsa/{enzyme}"
                                response = requests.get(url, timeout=10)
                                time.sleep(0.3)
                                
                                if response.status_code == 200 and response.text.strip():
                                    for line in response.text.strip().split('\n'):
                                        gene = line.split('\t')[1]
                                        genes_found.add((gene, enzyme, reaction))
                            except:
                                continue
                except:
                    continue
        
        return genes_found
    except Exception as e:
        print(f"Error getting genes for {compound_id}: {e}")
        return set()

def get_gene_info(gene_id):
    """Get gene information from KEGG"""
    try:
        url = f"http://rest.kegg.jp/get/{gene_id}"
        response = requests.get(url, timeout=10)
        time.sleep(0.3)
        
        if response.status_code == 200 and response.text.strip():
            lines = response.text.strip().split('\n')
            gene_name = None
            for line in lines:
                if line.startswith('SYMBOL'):
                    gene_name = line.split()[1] if len(line.split()) > 1 else None
                    break
            return gene_name
        return None
    except:
        return None

# Main processing
print("="*80)
print("KEGG METABOLITE-GENE PATHWAY ANALYSIS")
print("="*80)

print("\nLoading data files...")
genes_df = pd.read_csv(GENE_FILE)
metabolites_df = pd.read_csv(METABOLITE_FILE)

unique_metabolites = metabolites_df['Metabolite Name'].unique()
unique_genes = genes_df['Gene Name'].unique()

print(f"Loaded {len(unique_metabolites)} unique metabolites")
print(f"Loaded {len(unique_genes)} unique genes")

# Dictionary to store metabolite -> genes mapping
metabolite_to_genes = {}
gene_to_metabolites = defaultdict(list)
all_connections = []

print("\n" + "="*80)
print("STEP 1: Finding KEGG genes for each metabolite")
print("="*80)

for idx, metabolite in enumerate(unique_metabolites, 1):
    print(f"\n[{idx}/{len(unique_metabolites)}] Processing: {metabolite}")
    
    # Get KEGG compound IDs
    compound_ids = query_kegg_compound(metabolite)
    
    if not compound_ids:
        print(f"  ✗ No KEGG compound found")
        continue
    
    print(f"  Found {len(compound_ids)} KEGG compound ID(s): {', '.join(compound_ids[:3])}")
    
    all_genes = set()
    # Get genes for each compound ID (limit to first 2 to save time)
    for compound_id in compound_ids[:2]:
        genes = get_genes_for_compound(compound_id)
        all_genes.update(genes)
    
    if all_genes:
        print(f"  ✓ Found {len(all_genes)} gene-enzyme connections")
        metabolite_to_genes[metabolite] = all_genes
        
        # Store connections
        for gene_id, enzyme, reaction in all_genes:
            gene_name = get_gene_info(gene_id)
            
            connection = {
                'Metabolite': metabolite,
                'KEGG_Compound': compound_ids[0],
                'KEGG_Gene_ID': gene_id,
                'Gene_Symbol': gene_name if gene_name else gene_id,
                'Enzyme': enzyme,
                'Reaction': reaction
            }
            all_connections.append(connection)
            gene_to_metabolites[gene_id].append(metabolite)
    else:
        print(f"  ✗ No genes found")

print("\n" + "="*80)
print("STEP 2: Analyzing results")
print("="*80)

# Save all metabolite-gene connections
connections_df = pd.DataFrame(all_connections)
if not connections_df.empty:
    connections_output = os.path.join(OUTPUT_DIR, "metabolite_gene_connections.csv")
    connections_df.to_csv(connections_output, index=False)
    print(f"\n✓ Saved {len(connections_df)} metabolite-gene connections")
    print(f"  File: {connections_output}")

# Find genes that connect multiple metabolites
multi_metabolite_genes = []
for gene_id, metabolites in gene_to_metabolites.items():
    if len(metabolites) > 1:
        gene_name = get_gene_info(gene_id)
        multi_metabolite_genes.append({
            'KEGG_Gene_ID': gene_id,
            'Gene_Symbol': gene_name if gene_name else gene_id,
            'Number_of_Metabolites': len(metabolites),
            'Connected_Metabolites': '; '.join(set(metabolites))
        })

if multi_metabolite_genes:
    multi_genes_df = pd.DataFrame(multi_metabolite_genes)
    multi_genes_df = multi_genes_df.sort_values('Number_of_Metabolites', ascending=False)
    multi_genes_output = os.path.join(OUTPUT_DIR, "genes_connecting_multiple_metabolites.csv")
    multi_genes_df.to_csv(multi_genes_output, index=False)
    print(f"\n✓ Saved {len(multi_genes_df)} genes connecting multiple metabolites")
    print(f"  File: {multi_genes_output}")
    
    print("\n  Top 5 genes connecting most metabolites:")
    for idx, row in multi_genes_df.head().iterrows():
        print(f"    • {row['Gene_Symbol']} ({row['KEGG_Gene_ID']}): {row['Number_of_Metabolites']} metabolites")

# Find metabolites that share common genes (metabolites in same pathways)
metabolite_pairs_sharing_genes = []
metabolite_list = list(metabolite_to_genes.keys())

print("\n" + "="*80)
print("STEP 3: Finding metabolites connected through common genes/pathways")
print("="*80)

for i, met1 in enumerate(metabolite_list):
    for met2 in metabolite_list[i+1:]:
        genes1 = set([g[0] for g in metabolite_to_genes.get(met1, [])])
        genes2 = set([g[0] for g in metabolite_to_genes.get(met2, [])])
        
        common_genes = genes1 & genes2
        
        if common_genes:
            for gene_id in common_genes:
                gene_name = get_gene_info(gene_id)
                metabolite_pairs_sharing_genes.append({
                    'Metabolite_1': met1,
                    'Metabolite_2': met2,
                    'Shared_Gene_ID': gene_id,
                    'Shared_Gene_Symbol': gene_name if gene_name else gene_id,
                    'Total_Common_Genes': len(common_genes)
                })

if metabolite_pairs_sharing_genes:
    pairs_df = pd.DataFrame(metabolite_pairs_sharing_genes)
    pairs_output = os.path.join(OUTPUT_DIR, "metabolite_pairs_sharing_pathways.csv")
    pairs_df.to_csv(pairs_output, index=False)
    print(f"\n✓ Saved {len(pairs_df)} metabolite pairs sharing common genes")
    print(f"  File: {pairs_output}")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"Total metabolites with KEGG data: {len(metabolite_to_genes)}")
print(f"Total metabolite-gene connections: {len(all_connections)}")
print(f"Genes connecting multiple metabolites: {len(multi_metabolite_genes)}")
print(f"Metabolite pairs sharing pathways: {len(metabolite_pairs_sharing_genes)}")
print(f"\nAll output files saved in: {OUTPUT_DIR}")
print("="*80)
