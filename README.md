# HexaGene ClinVar Analysis

Sequence-intrinsic biophysical scoring for variant pathogenicity prediction.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains analysis code for the manuscript:

> **Sequence-Intrinsic Thermodynamic and Mechanical Constraints Resolve Clinically Classified Missense Variants Unresolved by Conservation-Based Predictors**
>
> Bachani S. bioRxiv 2024. DOI: [10.1101/XXXX.XX.XX.XXXXXX]

## Key Findings

| Result | Value |
|--------|-------|
| Overall discrimination | AUC = 0.72, T = 13.20 |
| Grey zone (REVEL 0.4-0.6) | AUC = 0.67, T = 3.62, p = 5.1×10⁻⁴ |
| Orthogonality to REVEL | r = 0.27 |
| Fragile gene advantage | +40% (T = 10.30 vs 7.33) |

## Installation

```bash
git clone https://github.com/[username]/hexagene-clinvar.git
cd hexagene-clinvar
pip install -r requirements.txt
```

## Quick Start

```python
from src.scoring import HexaGeneScorer

# Initialize
scorer = HexaGeneScorer()

# Score a variant
result = scorer.score_variant(
    ref_hexamer='ATCGAT',
    alt_hexamer='ATCTAT',
    codon_position=2
)

print(f"SDS Score: {result['sds']:.2f}")
print(f"L1 (Nucleotide): {result['l1']:.2f}")
print(f"L7 (Harmony): {result['l7']:.2f}")
```

## Biophysical Features

| Feature | Symbol | Biophysical Basis | Pathogenic Direction |
|---------|--------|-------------------|---------------------|
| Nucleotide Physics | L1 | Base-stacking ΔG | Higher |
| Kinetic Flow | L2 | Ribosomal kinetics | Higher |
| Local Stiffness | L4 | Persistence length | Lower |
| Codon Position | L5 | Translation impact | Higher |
| Harmonic Balance | L7 | Helix geometry | Lower |
| Complexity | L8 | Shannon entropy | Lower |
| Conflict Rate | - | Boundary transitions | Higher |

## Repository Structure

```
hexagene-clinvar/
├── README.md
├── LICENSE
├── requirements.txt
├── src/
│   ├── __init__.py
│   ├── hexamer.py          # Hexamer encoding (6-bit binary)
│   ├── features.py         # L1-L8 feature calculations
│   ├── scoring.py          # Composite SDS scoring
│   ├── gene_atlas.py       # Gene structural profiles
│   └── analysis.py         # Statistical analyses
├── data/
│   ├── README.md           # ClinVar/dbNSFP download instructions
│   └── gene_structural_atlas.csv
├── notebooks/
│   ├── 01_data_preparation.ipynb
│   ├── 02_main_analysis.ipynb
│   └── 03_figure_generation.ipynb
└── scripts/
    ├── score_variants.py   # Batch scoring script
    └── run_analysis.py     # Reproduce paper results
```

## Data Sources

### ClinVar Variants
Download from: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/

Required files:
- `variant_summary.txt.gz` (January 2024 or later)

### REVEL Scores
Download from dbNSFP: https://sites.google.com/site/jpaborern/dbNSFP

Required version: v4.4a or later

### Reference Genome
- GRCh38.p14 from Ensembl or UCSC

## Reproducing Paper Results

```bash
# 1. Download data (see data/README.md)

# 2. Prepare dataset
python scripts/prepare_data.py

# 3. Run main analysis
python scripts/run_analysis.py

# 4. Generate figures
python scripts/make_figures.py
```

## Feature Calculations

### L1: Nucleotide Transition Severity

```python
def calculate_l1(ref_base, alt_base):
    """
    Transversion = 1.0 (higher ΔG perturbation)
    Transition = 0.5 (lower ΔG perturbation)
    """
    purines = {'A', 'G'}
    ref_pur = ref_base in purines
    alt_pur = alt_base in purines
    return 1.0 if ref_pur != alt_pur else 0.5
```

### L4: Local Stiffness

```python
def calculate_l4(hexamer):
    """GC content as persistence length proxy."""
    return sum(1 for b in hexamer if b in 'GC') / 6
```

### L7: Harmonic Balance

```python
def calculate_l7(hexamer):
    """Purine-pyrimidine symmetry of central 4 bases."""
    core = hexamer[1:5]
    pattern = ''.join('R' if b in 'AG' else 'Y' for b in core)
    symmetric = ['RRYY', 'YYRR', 'RYRY', 'YRYR']
    return 1.0 if pattern in symmetric else 0.0
```

### Composite Score (SDS)

```python
def calculate_sds(features):
    """
    Structural Damage Score = Σ(weight × feature × sign)
    """
    weights = {'l1': 9.51, 'l4': 7.70, 'l7': 7.36, ...}
    signs = {'l1': +1, 'l4': -1, 'l7': -1, ...}
    
    return sum(weights[f] * features[f] * signs[f] for f in features)
```

## Gene Structural Atlas

Pre-computed structural profiles for 20,242 human genes:

```python
import pandas as pd

atlas = pd.read_csv('data/gene_structural_atlas.csv')

# Filter fragile genes
fragile = atlas[atlas['category'] == 'Fragile']
print(f"Fragile genes: {len(fragile)}")
# Examples: TTN, BRCA2, DMD, RYR1
```

## Citation

```bibtex
@article{bachani2024sequence,
  title={Sequence-Intrinsic Thermodynamic and Mechanical Constraints 
         Resolve Clinically Classified Missense Variants Unresolved 
         by Conservation-Based Predictors},
  author={Bachani, Sharad},
  journal={bioRxiv},
  year={2024},
  doi={10.1101/XXXX.XX.XX.XXXXXX}
}
```

## License

- **Code:** MIT License
- **Data:** CC-BY 4.0
- **Commercial use:** Contact author for licensing

## Related

- [Zenodo Supplementary Materials](https://doi.org/10.5281/zenodo.XXXXXXX)
- [bioRxiv Preprint](https://doi.org/10.1101/XXXX.XX.XX.XXXXXX)

## Contact

- **Author:** Sharad Bachani
- **Affiliation:** Merlin Digital, Dubai, UAE
- **Email:** [email]
