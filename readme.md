# HexaGene Validation Dataset

Supplementary data and analysis scripts for:

**"HexaGene: A Sequence-Intrinsic Biophysical Score for Variant Pathogenicity 
Prediction in Machine Learning Grey Zones"**

## Quick Start
```bash
pip install pandas scipy scikit-learn matplotlib seaborn
python analysis_script.py
```

## Expected Output
```
HEXAGENE VALIDATION RESULTS
==================================================
Total Variants: 2000
Pathogenic: 1000, Benign: 1000

1. GLOBAL PERFORMANCE
   T-statistic: 13.20
   ROC AUC:     0.67

2. ORTHOGONALITY
   Correlation (r) with REVEL: 0.27

3. GREY ZONE RESCUE (REVEL 0.4-0.6)
   Variants in Zone: 145
   HexaGene T-stat: 3.62
   HexaGene AUC:    0.67
```

## Data Description

| Column | Description |
|--------|-------------|
| Chrom | Chromosome |
| Pos | Genomic position (GRCh38) |
| Ref/Alt | Reference and alternate alleles |
| ClinVar_Label | 1 = Pathogenic, 0 = Benign |
| REVEL_Score | REVEL v1.3 pathogenicity score |
| HexaGene_Structural_Score | HexaGene composite physics score |

## Citation

[bioRxiv preprint link - pending]

## License

MIT