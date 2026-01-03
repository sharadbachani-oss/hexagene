"""
HexaGene ClinVar Analysis Package
=================================

Sequence-intrinsic biophysical scoring for variant pathogenicity prediction.

Main classes:
    HexaGeneScorer - Score variants using biophysical features

Main functions:
    calculate_sds - Calculate composite Structural Damage Score
    run_discrimination_analysis - Test pathogenic vs benign separation
    analyze_grey_zone - Analyze REVEL grey zone performance
"""

from .scoring import (
    HexaGeneScorer,
    ScoringResult,
    calculate_sds,
    calculate_all_features,
    calculate_l1_nucleotide,
    calculate_l4_stiffness,
    calculate_l7_harmony,
    calculate_l8_complexity,
    hexamer_to_binary,
)

from .analysis import (
    run_discrimination_analysis,
    analyze_grey_zone,
    calculate_correlation,
    run_gene_stratified_analysis,
    run_cross_validation,
    analyze_feature_contributions,
    calculate_clinical_utility,
)

__version__ = "1.0.0"
__author__ = "Sharad Bachani"
