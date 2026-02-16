# Metabolite-Gene Association Pipeline

This repository contains all Python implementations used to build a complete metabolite-to-gene association workflow. It includes MS/MS analogue prediction, SMILES validation, reverse biotransformation (RetroRules + BioTransformer), EC→gene mapping (UniProt API), and optional network table generation for visualization.

Both the **original pipeline** and the **validation pipeline** (for literature-derived metabolites) are provided exactly as used in the manuscript.

---

## Table of Contents

- [1. System Requirements](#1-system-requirements)
- [2. Creating the Required Environments](#2-creating-the-required-environments)
- [3. Repository Structure](#3-repository-structure)
- [4. Running the Pipelines](#4-running-the-pipelines)
- [5. Input File Requirements](#5-input-file-requirements)
- [6. Troubleshooting](#6-troubleshooting)
- [7. Citation](#7-citation)
- 
---

## 1. System Requirements

Multiple Python environments are required due to strict version dependencies.

| Environment | Python Version | Purpose |
|------------|----------------|---------|
| `ms2query-env` | **3.10.8–3.10.18** | MS2Query inference, MatchMS |
| `rdkit_env` | **3.10.x** | RDKit validation + RetroRules reverse biotransformation |
| `py313_env` | **3.13.3** | Prerequisites, validation pipeline|
| `py312_env` | **3.12.7** | Cytoscape-style metabolite–gene table formatting, UniProt mapping |

### Software Dependencies

- **Miniconda/Anaconda** (for environment management)
- **Python 3.10.x, 3.12.x, 3.13.x** (as specified above)
- **Internet connection** (for UniProt API queries)

---

## 2. Creating the Required Environments

First install **Miniconda**:  
[https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

### MS2Query Environment
```bash
conda create -n ms2query-env python=3.10.18
conda activate ms2query-env
pip install ms2query==1.5.4 matchms==0.26.4
```

### RDKit Environment
```bash
conda create -n rdkit_env python=3.10.18
conda activate rdkit_env
conda install -c conda-forge rdkit
```

### Python 3.13.3 Environment

(Prerequisites + validation pipeline)
```bash
conda create -n py313_env python=3.13.3
conda activate py313_env
pip install requests pandas numpy lxml simplejson
```

### Python 3.12.7 Environment

(Cytoscape network table generation + UniProt mapping)
```bash
conda create -n py312_env python=3.12.7
conda activate py312_env
pip install pandas numpy
```

---

## 3. Repository Structure
```
metabolite-gene-pipeline/
│
├── prerequisites/
│   ├── numpy_dependencies.py
│   ├── gensim_installation.py
│   └── cleanup_if_not_running.py
│
├── original_pipeline/
│   ├── ms2query_prediction.py
│   ├── rdkit_validation.py
│   ├── retrorules_reverse.py
│   └── ec_to_gene.py
│
├── validation_pipeline/
│   ├── Literature_metabolites_cleaned/
│   ├── pandas_requests/
│   ├── ready_for_smiles/
│   └── smiles_prediction.py
│
├── metabolite-gene_pairs_data_organisation.py
│
└── README.md
```

---

## 4. Running the Pipelines

### A. Prerequisites (Run First)

**Environment:** `py313_env` (Python 3.13.3)
```bash
conda activate py313_env
Then run:- 
python pre-requisites/numpy_dependencies.py
python pre-requisites/gensim_installation.py
python pre-requisites/cleanup_if_not_running.py   # optional
```

**These scripts ensure:**
- NumPy/Gensim compatibility
- MS2Query and MatchMS model compatibility
- Cleanup if MS2Query cannot import correctly

---

### B. Original Pipeline

#### Step 1 - MS/MS Analogue Prediction (MS2Query)

**Environment:** `ms2query-env`  
1. Before running ms2query_prediction.py, users must download the official MS2Query library files.
Download the MS2 Library (Version V8)

All required files are available at the following Zenodo DOI:

https://doi.org/10.5281/zenodo.13348638
 (Version V8)

Download all 8 files provided in this Zenodo record, as MS2Query requires the complete set for correct operation.

2. After downloading:
Create a directory named: ms2_library/
Place all 8 downloaded files from Zenodo inside this folder.

3. Before running the script, update the path of MS2Lib inside `ms2query_prediction.py`:
Important: Update the Library Path in the Script
Inside ms2query_prediction.py, update the directory path pointing to your MS2 library:
LIBRARY_DIR = "path/to/ms2_library/"
Make sure this path correctly matches the location where you saved the 8 files.


Then run:-

```bash
conda activate ms2query-env
python ms2query_prediction/ms2query_predict.py
```

**Input:**
- `.mgf` file containing MS/MS spectra (path must be specified inside the script)

**Output:**  
CSV containing:
- `spectrum_id`
- `analogue_score`
- Predicted SMILES
- InChIKey
- Compound name
- Rank (top 20 analogues)

---

#### Step 2 - RDKit SMILES Validation

**Environment:** `rdkit_env`
```bash
conda activate rdkit_env
python Biotransformation_pipeline/rdkit_validation.py
```

**Input:**  
CSV from Step 1 (predicted SMILES)

**Output:**  
CSV with:
- Validated SMILES
- InChIKey
- Removed invalid structures

---

#### Step 3 - Reverse Biotransformation (RetroRules + BioTransformer)

1. RetroRules File Required

The workflow requires the RetroRules pre-parsed reaction rule file
retrorules_rr02_rp2_flat_all.csv.
This file must be downloaded manually from the RetroRules archive.

Download Procedure:

a. Visit: https://retrorules.org/dl
b. Scroll to the Archives section.
c. Open the rr02 entry.
d. Under MNX-based Dataset (For RetroPath 2.0, Radii 1–8, Explicit Hydrogens), click Download.
   This will redirect to a Zenodo archive containing a .tar.gz file.
e. Download and extract the archive.
f. Inside the extracted folder, locate:
     retrorules_rr02_rp2_flat_all.csv
g. Place this file into a directory of your choice, for example:
    retrorules/
h. Open your RetroRules-based step script:
    original_pipeline/retrorules_reverse.py
i. Update the RetroRules file path:
    retrorules_path= "resources/retrorules/retrorules_rr02_rp2_flat_all.csv"

Only after this configuration should users run the RetroRules-based reverse biotransformation step.

**Still in `rdkit_env`:**
```bash
python Biotransformation_pipeline/retrorules_reverse.py
```

**Input:**  
CSV with SMILES from RDKit validation

**Output:**
CSV with
- Predicted precursor molecules
- EC numbers
- Reaction rules
- Specificity and confidence categories

### Sensitivity Analysis (Optional)
To reproduce the threshold selection analysis:
```bash
python Biotransformation_pipeline/retrorules_reverse_sensitivity_analysis.py
```
Tests thresholds [0.0, 0.3, 0.5, 0.7, 0.9] and generates comparison statistics.

---

#### Step 4 - EC → Gene Mapping (UniProt API)

**Environment:** `py312_env`
```bash
conda activate py312_env
python ec_to_gene/ec_to_gene.py
```

**What the script does:**
- Hierarchical EC matching
- Swiss-Prot reviewed entries = Perfect match
- TrEMBL entries = Close match
- Local caching to avoid repeated API queries

**Output:**  
Final CSV containing:
- EC number
- Gene name
- Protein name
- UniProt accession
- Match type (Perfect / Close / No match)
- Precursor molecule link

Species-Specific Gene Mapping: Implementation Guide

The workflow implements a two-stage strategy for handling species specificity:

Stage 1 – Precursor Prediction:
Reverse biotransformation is performed using the RetroRules universal database or the human database, which are the only publicly available versions currently available.

Stage 2 – EC → Gene Mapping:
Organism-specific filtering is applied during UniProt-based gene annotation.

This design maximises metabolite coverage during precursor inference while enabling taxonomic contextualisation at the gene and pathway annotation stage.

Default Behaviour:- 

Without modification:
EC-to-UniProt mapping queries all organisms (cross-species).

Returned results:
Genes and proteins from any organism annotated with the corresponding EC number.

Recommended use case:
Exploratory metabolomics, cross-kingdom analyses, or datasets where organism context is unknown or mixed.

Implementing Species-Specific Filtering

Species restriction can be enabled by adding organism-level filters to the UniProt query.

Step 1: Choose the Target Organism

Add the following configuration block near the top of the UniProt mapping script(ec_to_gene.py) (after imports):

# ============================================================================
# ORGANISM FILTER CONFIGURATION
# ============================================================================

ORGANISM_ID = "9606"     # Human (Homo sapiens) – DEFAULT
# ORGANISM_ID = "83333"  # Escherichia coli K-12
# ORGANISM_ID = "559292" # Yeast (Saccharomyces cerevisiae S288C)
# ORGANISM_ID = "10090"  # Mouse (Mus musculus)
# ORGANISM_ID = "10116"  # Rat (Rattus norvegicus)
# ORGANISM_ID = "3702"   # Arabidopsis thaliana
# ORGANISM_ID = "7227"   # Drosophila melanogaster
# ORGANISM_ID = None     # All organisms (disable filtering)

print(f"🔬 Organism filter active: {ORGANISM_ID}")



For finding organism taxonomy IDs:

-Visit: https://www.uniprot.org/taxonomy

-Search for the species name

-Use the numeric taxonomy identifier (e.g. Homo sapiens → 9606)

Step 2: Modify the UniProt Query Function(in ec_to_gene.py)

Original implementation:-

def query_uniprot(ec_query, reviewed=True):
    q = f"ec:{ec_query}"
    if reviewed:
        q += " AND reviewed:true"


Updated implementation(with organism filtering):-

def query_uniprot(ec_query, reviewed=True, organism=None):
    q = f"ec:{ec_query}"
    if reviewed:
        q += " AND reviewed:true"
    if organism:
        q += f" AND organism_id:{organism}"

Step 3: Update Function Calls

In the process_ec_number function, add the organism argument to both UniProt queries:

# Reviewed proteins
gene, protein = query_uniprot(q, reviewed=True, organism=ORGANISM_ID)

# Unreviewed proteins (fallback)
gene, protein = query_uniprot(q, reviewed=False, organism=ORGANISM_ID)

Step 4: Clear Cache When Changing Organisms

The script caches UniProt responses to avoid repeated API calls. When switching organism filters, clear the cache:

import os
import tempfile

cache_file = os.path.join(tempfile.gettempdir(), "uniprot_ec_cache.json")

if os.path.exists(cache_file):
    os.remove(cache_file)
    print("Cache cleared for new organism filter")

Common Organism Taxonomy IDs:- 
Organism	Taxonomy ID
Human	9606	
Mouse	10090	
Rat	10116	
E. coli K-12	83333	
Yeast (S. cerevisiae)	559292	
Arabidopsis	3702	
Zebrafish	7955
C. elegans	6239	

-Full taxonomy list: https://www.uniprot.org/taxonomy

Expected Impact of Species Filtering:-

-Reduced match rate:
Species-specific filtering may reduce gene matches by approximately 30–70% compared to cross-species queries.

-Increased biological relevance:
Returned genes reflect the metabolic capacity of the target organism.

-Improved pathway accuracy:
Enables organism-appropriate pathway enrichment and functional interpretation.

Troubleshooting

Issue: No gene matches after applying organism filter
Solution:

-Verify the taxonomy ID using the UniProt taxonomy browser

-Ensure the unreviewed protein fallback is enabled

-Some EC numbers may not be annotated in the target organism, which reflects true biological absence

Issue: Cache not updating after organism change
Solution:
Manually delete the cache file located at:
{tempfile.gettempdir()}/uniprot_ec_cache.json-

---

### C. Validation Pipeline (For Known Literature Metabolites)

**Environment:** `py313_env`
```bash
conda activate py313_env
python validation_pipeline/smiles_prediction.py
```

**Input:**  
Cleaned literature metabolites CSV

**Output:**
- Predicted SMILES

**Note:** These SMILES can be used as input for the Original Pipeline (starting from Step 2) to proceed to further steps of gene prediction(Respective environment to be used).
~ After getting the final output folder, run:-
**Still in `py313_env`:**
python validation_pipeline/sorted_gene_list_ENRICHR.py

The resultant list of genes can be used directly in the ENRICHR WEB SERVER.




---

### D. Optional - Cytoscape Network Table

**Environment:** `py312_env`
```bash
conda activate py312_env
python metabolite-gene_pairs_data_organisation.py
```

**Input requirements:**
- `nodes` column = list of metabolites
- `edges` column = list of genes

**Output:**
- CSV compatible with Cytoscape network import(Preferably convert it to Tab-limited Text file)
- Clean metabolite-gene edge table

---

## 5. Input File Requirements

### Original Pipeline

1. **Input 1:** `.mgf` MS/MS file
2. **Input 2:** CSV from MS2Query
3. **Input 3:** CSV from RDKit validation
4. **Input 4:** Output CSV of EC predictions

Each stage outputs a file that becomes the input of the next stage.

### Validation Pipeline

- Literature metabolite names CSV
- Predicted SMILES used
- ALL other output files produced starting from Step 2 of the Original Pipeline


---

## 6. Troubleshooting

### General Issues

- Always activate the correct environment before running a script
- Check file paths in the scripts
- Verify input formats (column names match expectations)
- Confirm Python version matches requirements

### MS2Query Import Errors

If MS2Query fails to import:
1. Rerun prerequisites scripts
2. Verify Python version is 3.10.8–3.10.18
3. Check that ms2query==1.5.4 and matchms==0.26.4 are installed

### UniProt API Errors

- May be due to rate limits
- Local cache prevents repeated calls
- Check internet connection
- Verify cache file location in `tempfile.gettempdir()`

### Environment Conflicts

If packages conflict:
1. Delete the environment: `conda env remove -n <env_name>`
2. Recreate using the commands in Section 2
3. Ensure you're using the exact Python versions specified

---

## 7. Citation

If you use this pipeline in your research, please cite:
```
[Manuscript citation will be added upon publication]
```

---

## 8. License

This project is licensed under the Apache License, Version 2.0.  
---

## 9. Acknowledgments

This pipeline integrates tools and databases including:
- MS2Query
- MatchMS
- RDKit
- RetroRules
- BioTransformer
- UniProt

We thank the developers and maintainers of these resources.
