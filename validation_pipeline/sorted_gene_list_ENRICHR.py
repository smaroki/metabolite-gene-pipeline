import pandas as pd
import os
import re

# --- INPUT AND OUTPUT PATHS ---
input_file = r"C:\Users\dell\OneDrive\Desktop\Raw files\results\precursor_output\EC_to_Genes_Uniprot_Annotated.csv"
output_folder = r"C:\Users\dell\OneDrive\Desktop\Raw files\results\sorted_gene_list"
output_file = os.path.join(output_folder, "unique_sorted_genes.csv")

# --- CREATE OUTPUT FOLDER IF IT DOESN'T EXIST ---
os.makedirs(output_folder, exist_ok=True)

# --- READ CSV FILE ---
df = pd.read_csv(input_file)

# --- EXTRACT AND CLEAN GENE NAMES ---
gene_series = df['Gene Name'].dropna().astype(str)
all_genes = []

for entry in gene_series:
    # Split on comma, semicolon, or any whitespace
    genes = re.split(r'[,\s;]+', entry)
    genes = [gene.strip() for gene in genes if gene.strip() != '']
    all_genes.extend(genes)

# --- REMOVE DUPLICATES AND SORT ---
unique_genes = sorted(set(all_genes))

# --- SAVE TO CSV ---
df_out = pd.DataFrame({'Gene Name': unique_genes})
df_out.to_csv(output_file, index=False)

print(f"✅ Cleaned and sorted unique gene list saved to CSV:\n{output_file}")
