# Data Download Instructions

This directory should contain the source data files for analysis.

## Required Data Sources

### 1. ClinVar Variants

**Source:** https://ftp.ncbi.nlm.nih.gov/pub/clinvar/

**Required file:** `variant_summary.txt.gz`

**Download:**
```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip variant_summary.txt.gz
```

**Filtering criteria:**
- Clinical significance: Pathogenic, Likely pathogenic, Benign, Likely benign
- Review status: >= 2 stars
- Molecular consequence: missense_variant
- Assembly: GRCh38

### 2. REVEL Scores

**Source:** dbNSFP v4.4a

**Download:** https://sites.google.com/site/jpopgen/dbNSFP

**Required columns:**
- chr, pos, ref, alt
- REVEL_score

**Note:** dbNSFP is large (~35GB compressed). You may also use the REVEL-specific file if available.

### 3. Reference Genome

**Source:** GRCh38.p14

**Download from Ensembl:**
```bash
wget https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

**Index with samtools:**
```bash
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
```

## Provided Files

### gene_structural_atlas.csv

Pre-computed structural profiles for 20,242 human protein-coding genes.

**Columns:**
| Column | Description |
|--------|-------------|
| gene_symbol | HGNC gene symbol |
| ensembl_id | Ensembl gene ID |
| cds_length | Coding sequence length (bp) |
| stiffness_mean | Mean local stiffness |
| stiffness_sd | SD of local stiffness |
| harmony_mean | Mean harmonic balance |
| harmony_sd | SD of harmonic balance |
| complexity_mean | Mean complexity |
| conflict_mean | Mean conflict rate |
| structural_risk | Composite risk score |
| category | Fragile/Intermediate/Robust |

## Data Preparation Script

```python
import pandas as pd

# Load ClinVar
clinvar = pd.read_csv('variant_summary.txt', sep='\t')

# Filter
clinvar = clinvar[
    (clinvar['ClinicalSignificance'].str.contains('athogenic|enign', case=False)) &
    (clinvar['ReviewStatus'].str.contains('2|3|4')) &
    (clinvar['Assembly'] == 'GRCh38') &
    (clinvar['Type'] == 'single nucleotide variant')
]

# Merge with REVEL
# ... see scripts/prepare_data.py for complete pipeline
```

## Citation

When using ClinVar data, cite:
> Landrum MJ, et al. ClinVar: improvements to accessing data. 
> Nucleic Acids Res. 2020;48(D1):D835-D844.

When using REVEL, cite:
> Ioannidis NM, et al. REVEL: An ensemble method for predicting 
> the pathogenicity of rare missense variants. Am J Hum Genet. 2016.
