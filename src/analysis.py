"""
Statistical Analysis Functions for HexaGene ClinVar Validation.

Provides functions to reproduce the key analyses from the manuscript.
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold
from typing import Dict, Tuple, List


def calculate_cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """
    Calculate Cohen's d effect size.
    
    Parameters
    ----------
    group1 : array-like
        First group values
    group2 : array-like
        Second group values
    
    Returns
    -------
    float
        Cohen's d effect size
    """
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    
    return (np.mean(group1) - np.mean(group2)) / pooled_std


def run_discrimination_analysis(
    df: pd.DataFrame,
    score_col: str = 'sds_score',
    label_col: str = 'is_pathogenic'
) -> Dict:
    """
    Run discrimination analysis comparing pathogenic vs benign.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with scores and labels
    score_col : str
        Column name for scores
    label_col : str
        Column name for binary labels
    
    Returns
    -------
    dict
        Discrimination statistics
    """
    pathogenic = df[df[label_col] == 1][score_col].values
    benign = df[df[label_col] == 0][score_col].values
    
    # T-test
    t_stat, p_value = stats.ttest_ind(pathogenic, benign, equal_var=False)
    
    # Effect size
    cohens_d = calculate_cohens_d(pathogenic, benign)
    
    # AUC
    auc = roc_auc_score(df[label_col], df[score_col])
    
    # ROC curve
    fpr, tpr, thresholds = roc_curve(df[label_col], df[score_col])
    
    return {
        'n_pathogenic': len(pathogenic),
        'n_benign': len(benign),
        'path_mean': np.mean(pathogenic),
        'path_std': np.std(pathogenic),
        'ben_mean': np.mean(benign),
        'ben_std': np.std(benign),
        't_statistic': t_stat,
        'p_value': p_value,
        'cohens_d': cohens_d,
        'auc': auc,
        'roc_curve': (fpr, tpr, thresholds)
    }


def analyze_grey_zone(
    df: pd.DataFrame,
    revel_col: str = 'revel_score',
    score_col: str = 'sds_score',
    label_col: str = 'is_pathogenic',
    grey_low: float = 0.4,
    grey_high: float = 0.6
) -> Dict:
    """
    Analyze discrimination within REVEL grey zone.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with REVEL scores, SDS scores, and labels
    revel_col : str
        Column name for REVEL scores
    score_col : str
        Column name for SDS scores
    label_col : str
        Column name for binary labels
    grey_low : float
        Lower bound of grey zone
    grey_high : float
        Upper bound of grey zone
    
    Returns
    -------
    dict
        Grey zone analysis results
    """
    # Filter to grey zone
    grey_zone = df[(df[revel_col] >= grey_low) & (df[revel_col] <= grey_high)]
    
    if len(grey_zone) < 10:
        return {'error': 'Insufficient grey zone variants'}
    
    # Run discrimination analysis
    results = run_discrimination_analysis(grey_zone, score_col, label_col)
    results['n_grey_zone'] = len(grey_zone)
    results['grey_zone_fraction'] = len(grey_zone) / len(df)
    
    return results


def calculate_correlation(
    df: pd.DataFrame,
    col1: str,
    col2: str,
    method: str = 'pearson'
) -> Tuple[float, float]:
    """
    Calculate correlation between two columns.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data
    col1 : str
        First column name
    col2 : str
        Second column name
    method : str
        'pearson' or 'spearman'
    
    Returns
    -------
    tuple
        (correlation coefficient, p-value)
    """
    valid = df[[col1, col2]].dropna()
    
    if method == 'pearson':
        r, p = stats.pearsonr(valid[col1], valid[col2])
    else:
        r, p = stats.spearmanr(valid[col1], valid[col2])
    
    return r, p


def run_gene_stratified_analysis(
    df: pd.DataFrame,
    gene_atlas: pd.DataFrame,
    score_col: str = 'sds_score',
    label_col: str = 'is_pathogenic',
    gene_col: str = 'gene'
) -> Dict:
    """
    Run analysis stratified by gene structural category.
    
    Parameters
    ----------
    df : pd.DataFrame
        Variant data
    gene_atlas : pd.DataFrame
        Gene structural atlas with 'category' column
    score_col : str
        Score column name
    label_col : str
        Label column name
    gene_col : str
        Gene column name
    
    Returns
    -------
    dict
        Results by gene category
    """
    # Merge gene categories
    df_merged = df.merge(
        gene_atlas[['gene_symbol', 'category']],
        left_on=gene_col,
        right_on='gene_symbol',
        how='left'
    )
    
    results = {}
    
    for category in ['Fragile', 'Intermediate', 'Robust']:
        subset = df_merged[df_merged['category'] == category]
        
        if len(subset) < 20:
            results[category] = {'error': 'Insufficient variants'}
            continue
        
        results[category] = run_discrimination_analysis(
            subset, score_col, label_col
        )
        results[category]['n_variants'] = len(subset)
    
    # Calculate fragile vs robust ratio
    if 'Fragile' in results and 'Robust' in results:
        if 't_statistic' in results['Fragile'] and 't_statistic' in results['Robust']:
            results['fragile_robust_ratio'] = (
                results['Fragile']['t_statistic'] / 
                results['Robust']['t_statistic']
            )
    
    return results


def run_cross_validation(
    df: pd.DataFrame,
    score_col: str = 'sds_score',
    label_col: str = 'is_pathogenic',
    n_folds: int = 5,
    random_state: int = 42
) -> Dict:
    """
    Run stratified k-fold cross-validation.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with scores and labels
    score_col : str
        Score column name
    label_col : str
        Label column name
    n_folds : int
        Number of CV folds
    random_state : int
        Random seed
    
    Returns
    -------
    dict
        Cross-validation results
    """
    X = df[score_col].values.reshape(-1, 1)
    y = df[label_col].values
    
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=random_state)
    
    fold_results = []
    
    for fold, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        train_auc = roc_auc_score(y[train_idx], X[train_idx].ravel())
        test_auc = roc_auc_score(y[test_idx], X[test_idx].ravel())
        
        # T-test on test fold
        test_path = X[test_idx][y[test_idx] == 1].ravel()
        test_ben = X[test_idx][y[test_idx] == 0].ravel()
        t_stat, _ = stats.ttest_ind(test_path, test_ben, equal_var=False)
        
        fold_results.append({
            'fold': fold + 1,
            'train_auc': train_auc,
            'test_auc': test_auc,
            'test_t_stat': t_stat
        })
    
    fold_df = pd.DataFrame(fold_results)
    
    return {
        'fold_results': fold_df,
        'mean_train_auc': fold_df['train_auc'].mean(),
        'mean_test_auc': fold_df['test_auc'].mean(),
        'std_test_auc': fold_df['test_auc'].std(),
        'train_test_gap': fold_df['train_auc'].mean() - fold_df['test_auc'].mean()
    }


def analyze_feature_contributions(
    df: pd.DataFrame,
    feature_cols: List[str],
    label_col: str = 'is_pathogenic'
) -> pd.DataFrame:
    """
    Analyze individual feature contributions.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data with features and labels
    feature_cols : list
        List of feature column names
    label_col : str
        Label column name
    
    Returns
    -------
    pd.DataFrame
        Feature statistics sorted by absolute T-statistic
    """
    results = []
    
    for col in feature_cols:
        path_vals = df[df[label_col] == 1][col].values
        ben_vals = df[df[label_col] == 0][col].values
        
        t_stat, p_val = stats.ttest_ind(path_vals, ben_vals, equal_var=False)
        d = calculate_cohens_d(path_vals, ben_vals)
        
        results.append({
            'feature': col,
            'path_mean': np.mean(path_vals),
            'path_std': np.std(path_vals),
            'ben_mean': np.mean(ben_vals),
            'ben_std': np.std(ben_vals),
            't_statistic': t_stat,
            'p_value': p_val,
            'cohens_d': d,
            'direction': 'Path higher' if t_stat > 0 else 'Ben higher'
        })
    
    results_df = pd.DataFrame(results)
    results_df['abs_t'] = results_df['t_statistic'].abs()
    results_df = results_df.sort_values('abs_t', ascending=False)
    
    return results_df


def calculate_clinical_utility(
    df: pd.DataFrame,
    score_col: str = 'sds_score',
    label_col: str = 'is_pathogenic',
    revel_col: str = 'revel_score',
    grey_low: float = 0.4,
    grey_high: float = 0.6,
    sds_thresholds: Tuple[float, float] = (1.5, 2.5)
) -> Dict:
    """
    Calculate clinical utility for grey zone reclassification.
    
    Parameters
    ----------
    df : pd.DataFrame
        Data
    score_col : str
        SDS score column
    label_col : str
        Label column
    revel_col : str
        REVEL score column
    grey_low, grey_high : float
        Grey zone bounds
    sds_thresholds : tuple
        (low, high) SDS thresholds for stratification
    
    Returns
    -------
    dict
        Clinical utility statistics
    """
    # Filter to grey zone
    grey = df[(df[revel_col] >= grey_low) & (df[revel_col] <= grey_high)]
    
    low_thresh, high_thresh = sds_thresholds
    
    # Stratify by SDS
    low_sds = grey[grey[score_col] < low_thresh]
    mid_sds = grey[(grey[score_col] >= low_thresh) & (grey[score_col] <= high_thresh)]
    high_sds = grey[grey[score_col] > high_thresh]
    
    results = {
        'n_grey_zone': len(grey),
        'low_sds': {
            'n': len(low_sds),
            'n_path': low_sds[label_col].sum(),
            'path_rate': low_sds[label_col].mean() if len(low_sds) > 0 else 0
        },
        'mid_sds': {
            'n': len(mid_sds),
            'n_path': mid_sds[label_col].sum(),
            'path_rate': mid_sds[label_col].mean() if len(mid_sds) > 0 else 0
        },
        'high_sds': {
            'n': len(high_sds),
            'n_path': high_sds[label_col].sum(),
            'path_rate': high_sds[label_col].mean() if len(high_sds) > 0 else 0
        }
    }
    
    # Calculate reclassification potential
    n_reclassifiable = len(low_sds) + len(high_sds)
    results['reclassification_potential'] = n_reclassifiable / len(grey) if len(grey) > 0 else 0
    
    return results
