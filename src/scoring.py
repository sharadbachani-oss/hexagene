"""
HexaGene ClinVar Scoring Module
===============================

Sequence-intrinsic biophysical scoring for variant pathogenicity prediction.

This module implements the core scoring engine described in the manuscript.
"""

import math
from typing import Dict, Tuple, Optional
from dataclasses import dataclass


# =============================================================================
# CONSTANTS
# =============================================================================

PURINES = {'A', 'G'}
PYRIMIDINES = {'C', 'T'}

# Feature weights (from training set T-statistics)
FEATURE_WEIGHTS = {
    'l1_nucleotide': 9.51,
    'l2_kinetic': 2.33,
    'l4_stiffness': 7.70,
    'l5_position': 4.12,
    'l7_harmony': 7.36,
    'l8_complexity': 6.23,
    'conflict_rate': 3.89,
}

# Direction of pathogenic association (+1 = higher in pathogenic, -1 = lower)
FEATURE_SIGNS = {
    'l1_nucleotide': +1,
    'l2_kinetic': +1,
    'l4_stiffness': -1,
    'l5_position': +1,
    'l7_harmony': -1,
    'l8_complexity': -1,
    'conflict_rate': +1,
}

# Symmetric patterns for harmonic balance
SYMMETRIC_PATTERNS = {
    'RRYY', 'YYRR', 'RYRY', 'YRYR', 
    'RRYR', 'YRYY', 'YRRY', 'RYYR'
}


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ScoringResult:
    """Container for scoring results."""
    sds: float
    features: Dict[str, float]
    hexamer_ref: str
    hexamer_alt: str
    is_transversion: bool
    codon_position: int


# =============================================================================
# FEATURE CALCULATIONS
# =============================================================================

def is_purine(base: str) -> bool:
    """Check if base is a purine (A or G)."""
    return base.upper() in PURINES


def is_pyrimidine(base: str) -> bool:
    """Check if base is a pyrimidine (C or T)."""
    return base.upper() in PYRIMIDINES


def calculate_l1_nucleotide(ref_base: str, alt_base: str) -> float:
    """
    Calculate nucleotide transition severity (L1).
    
    Based on nearest-neighbor thermodynamics:
    - Transversions (Pu <-> Py): Higher ΔG perturbation = 1.0
    - Transitions (Pu <-> Pu or Py <-> Py): Lower ΔG = 0.5
    
    Parameters
    ----------
    ref_base : str
        Reference nucleotide (A, C, G, or T)
    alt_base : str
        Alternate nucleotide
    
    Returns
    -------
    float
        Transition severity score (0.5 or 1.0)
    """
    ref_is_purine = is_purine(ref_base)
    alt_is_purine = is_purine(alt_base)
    
    if ref_is_purine != alt_is_purine:
        return 1.0  # Transversion
    return 0.5  # Transition


def calculate_l2_kinetic(hexamer_ref: str, hexamer_alt: str) -> float:
    """
    Calculate kinetic flow score (L2).
    
    Measures disruption to codon-boundary transitions.
    
    Parameters
    ----------
    hexamer_ref : str
        Reference 6-mer sequence
    hexamer_alt : str
        Alternate 6-mer sequence
    
    Returns
    -------
    float
        Kinetic flow disruption score (0.0-1.0)
    """
    # Check boundary transitions (positions 2-3 and 3-4 in 0-indexed)
    ref_transitions = 0
    alt_transitions = 0
    
    for i in range(5):
        ref_same = is_purine(hexamer_ref[i]) == is_purine(hexamer_ref[i+1])
        alt_same = is_purine(hexamer_alt[i]) == is_purine(hexamer_alt[i+1])
        
        if ref_same:
            ref_transitions += 1
        if alt_same:
            alt_transitions += 1
    
    # Return normalized difference
    return abs(alt_transitions - ref_transitions) / 5.0


def calculate_l4_stiffness(hexamer: str) -> float:
    """
    Calculate local stiffness (L4).
    
    GC content as proxy for DNA persistence length.
    Higher GC = stiffer DNA (~53nm vs ~45nm persistence length).
    
    Parameters
    ----------
    hexamer : str
        6-mer sequence
    
    Returns
    -------
    float
        Stiffness score (0.0-1.0, proportion GC)
    """
    gc_count = sum(1 for b in hexamer.upper() if b in 'GC')
    return gc_count / 6.0


def calculate_l5_position(codon_position: int) -> float:
    """
    Calculate codon position impact (L5).
    
    Position 1: Almost always changes amino acid
    Position 2: Always changes amino acid
    Position 3: Often synonymous (wobble position)
    
    Parameters
    ----------
    codon_position : int
        Position within codon (1, 2, or 3)
    
    Returns
    -------
    float
        Position impact weight
    """
    weights = {1: 1.0, 2: 0.8, 3: 0.2}
    return weights.get(codon_position, 0.5)


def calculate_l7_harmony(hexamer: str) -> float:
    """
    Calculate harmonic balance (L7).
    
    Purine-pyrimidine symmetry of central 4 bases ("nuclear core").
    Symmetric patterns maintain optimal helix geometry.
    
    Parameters
    ----------
    hexamer : str
        6-mer sequence
    
    Returns
    -------
    float
        Harmony score (1.0 if symmetric, 0.0 otherwise)
    """
    core = hexamer[1:5].upper()  # Positions 2,3,4,5 (1-indexed)
    pattern = ''.join('R' if b in PURINES else 'Y' for b in core)
    
    return 1.0 if pattern in SYMMETRIC_PATTERNS else 0.0


def calculate_l8_complexity(hexamer: str) -> float:
    """
    Calculate compositional complexity (L8).
    
    Shannon entropy of base composition.
    Low complexity (repeats) associated with structural issues.
    
    Parameters
    ----------
    hexamer : str
        6-mer sequence
    
    Returns
    -------
    float
        Normalized Shannon entropy (0.0-1.0)
    """
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for b in hexamer.upper():
        if b in counts:
            counts[b] += 1
    
    entropy = 0.0
    for count in counts.values():
        if count > 0:
            p = count / 6.0
            entropy -= p * math.log2(p)
    
    # Normalize to [0, 1] (max entropy = 2.0)
    return entropy / 2.0


def calculate_conflict_rate(codon1: str, codon2: str) -> float:
    """
    Calculate neighbor codon conflict rate.
    
    Checks for disrupted elemental transitions at codon boundaries.
    
    Parameters
    ----------
    codon1 : str
        First codon (3 bases)
    codon2 : str
        Second codon (3 bases)
    
    Returns
    -------
    float
        Conflict indicator (0.0 or 1.0)
    """
    if len(codon1) < 3 or len(codon2) < 3:
        return 0.0
    
    boundary1 = codon1[2]  # Last base of codon 1
    boundary2 = codon2[0]  # First base of codon 2
    
    # Conflict if same type creates "speed bump"
    same_type = is_purine(boundary1) == is_purine(boundary2)
    return 1.0 if same_type else 0.0


# =============================================================================
# COMPOSITE SCORING
# =============================================================================

def calculate_all_features(
    hexamer_ref: str,
    hexamer_alt: str,
    codon_position: int,
    ref_base: Optional[str] = None,
    alt_base: Optional[str] = None
) -> Dict[str, float]:
    """
    Calculate all biophysical features for a variant.
    
    Parameters
    ----------
    hexamer_ref : str
        Reference 6-mer centered on variant
    hexamer_alt : str
        Alternate 6-mer with variant
    codon_position : int
        Position within codon (1, 2, or 3)
    ref_base : str, optional
        Reference base (extracted from hexamer if not provided)
    alt_base : str, optional
        Alternate base (extracted from hexamer if not provided)
    
    Returns
    -------
    dict
        Dictionary of feature name -> value
    """
    # Extract central bases if not provided
    if ref_base is None:
        ref_base = hexamer_ref[2]  # 0-indexed position 2 = center
    if alt_base is None:
        alt_base = hexamer_alt[2]
    
    features = {
        'l1_nucleotide': calculate_l1_nucleotide(ref_base, alt_base),
        'l2_kinetic': calculate_l2_kinetic(hexamer_ref, hexamer_alt),
        'l4_stiffness': calculate_l4_stiffness(hexamer_alt),
        'l5_position': calculate_l5_position(codon_position),
        'l7_harmony': calculate_l7_harmony(hexamer_alt),
        'l8_complexity': calculate_l8_complexity(hexamer_alt),
    }
    
    # Conflict rate (using hexamer as two overlapping codons)
    if len(hexamer_alt) >= 6:
        codon1 = hexamer_alt[0:3]
        codon2 = hexamer_alt[3:6]
        features['conflict_rate'] = calculate_conflict_rate(codon1, codon2)
    else:
        features['conflict_rate'] = 0.0
    
    return features


def calculate_sds(features: Dict[str, float]) -> float:
    """
    Calculate Structural Damage Score (SDS).
    
    SDS = Σ(weight × feature × sign)
    
    Parameters
    ----------
    features : dict
        Feature name -> value mapping
    
    Returns
    -------
    float
        Composite SDS score
    """
    sds = 0.0
    for feature, value in features.items():
        weight = FEATURE_WEIGHTS.get(feature, 1.0)
        sign = FEATURE_SIGNS.get(feature, 1)
        sds += weight * value * sign
    
    return sds


# =============================================================================
# MAIN SCORER CLASS
# =============================================================================

class HexaGeneScorer:
    """
    HexaGene Biophysical Scorer for Variant Pathogenicity.
    
    Computes sequence-intrinsic features and composite damage score.
    
    Examples
    --------
    >>> scorer = HexaGeneScorer()
    >>> result = scorer.score_variant(
    ...     hexamer_ref='ATCGAT',
    ...     hexamer_alt='ATCTAT',
    ...     codon_position=2
    ... )
    >>> print(f"SDS: {result.sds:.2f}")
    """
    
    def __init__(
        self,
        weights: Optional[Dict[str, float]] = None,
        signs: Optional[Dict[str, int]] = None
    ):
        """
        Initialize scorer.
        
        Parameters
        ----------
        weights : dict, optional
            Custom feature weights (uses defaults if not provided)
        signs : dict, optional
            Custom feature signs (uses defaults if not provided)
        """
        self.weights = weights or FEATURE_WEIGHTS.copy()
        self.signs = signs or FEATURE_SIGNS.copy()
    
    def score_variant(
        self,
        hexamer_ref: str,
        hexamer_alt: str,
        codon_position: int,
        ref_base: Optional[str] = None,
        alt_base: Optional[str] = None
    ) -> ScoringResult:
        """
        Score a single variant.
        
        Parameters
        ----------
        hexamer_ref : str
            Reference 6-mer centered on variant
        hexamer_alt : str
            Alternate 6-mer with variant
        codon_position : int
            Position within codon (1, 2, or 3)
        ref_base : str, optional
            Reference base
        alt_base : str, optional
            Alternate base
        
        Returns
        -------
        ScoringResult
            Dataclass with SDS score and all features
        """
        # Extract bases if needed
        if ref_base is None:
            ref_base = hexamer_ref[2]
        if alt_base is None:
            alt_base = hexamer_alt[2]
        
        # Calculate features
        features = calculate_all_features(
            hexamer_ref, hexamer_alt, codon_position, ref_base, alt_base
        )
        
        # Calculate composite score
        sds = 0.0
        for feature, value in features.items():
            weight = self.weights.get(feature, 1.0)
            sign = self.signs.get(feature, 1)
            sds += weight * value * sign
        
        # Check if transversion
        is_tv = is_purine(ref_base) != is_purine(alt_base)
        
        return ScoringResult(
            sds=sds,
            features=features,
            hexamer_ref=hexamer_ref,
            hexamer_alt=hexamer_alt,
            is_transversion=is_tv,
            codon_position=codon_position
        )
    
    def score_batch(
        self,
        variants: list
    ) -> list:
        """
        Score multiple variants.
        
        Parameters
        ----------
        variants : list of dict
            Each dict must have: hexamer_ref, hexamer_alt, codon_position
        
        Returns
        -------
        list of ScoringResult
        """
        results = []
        for v in variants:
            result = self.score_variant(
                hexamer_ref=v['hexamer_ref'],
                hexamer_alt=v['hexamer_alt'],
                codon_position=v['codon_position'],
                ref_base=v.get('ref_base'),
                alt_base=v.get('alt_base')
            )
            results.append(result)
        return results


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def hexamer_to_binary(hexamer: str) -> int:
    """
    Convert 6-mer to binary integer (0-63).
    
    Purine (A/G) = 1, Pyrimidine (C/T) = 0
    
    Parameters
    ----------
    hexamer : str
        6-mer DNA sequence
    
    Returns
    -------
    int
        Integer encoding (0-63)
    """
    binary = 0
    for i, base in enumerate(hexamer.upper()):
        if base in PURINES:
            binary |= (1 << (5 - i))
    return binary


def binary_to_pattern(binary: int) -> str:
    """
    Convert binary integer to R/Y pattern.
    
    Parameters
    ----------
    binary : int
        Integer encoding (0-63)
    
    Returns
    -------
    str
        6-character R/Y pattern
    """
    pattern = ''
    for i in range(6):
        if binary & (1 << (5 - i)):
            pattern += 'R'
        else:
            pattern += 'Y'
    return pattern


# =============================================================================
# CLI
# =============================================================================

if __name__ == '__main__':
    import sys
    
    print("HexaGene ClinVar Scorer")
    print("=" * 50)
    
    # Example usage
    scorer = HexaGeneScorer()
    
    # Test variant: G>A transition in ATCGAT context
    result = scorer.score_variant(
        hexamer_ref='ATCGAT',
        hexamer_alt='ATCAAT',
        codon_position=2
    )
    
    print(f"\nExample: ATCGAT -> ATCAAT (codon position 2)")
    print(f"SDS Score: {result.sds:.2f}")
    print(f"Is Transversion: {result.is_transversion}")
    print(f"\nFeatures:")
    for name, value in result.features.items():
        print(f"  {name}: {value:.3f}")
