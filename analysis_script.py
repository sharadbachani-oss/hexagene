#!/usr/bin/env python3
"""
HexaGene Validation Analysis Script
======================================================================
Reproduces the statistical validation metrics and figures for:
"HexaGene: An Orthogonal Structural Physics Layer for VUS Resolution"

Usage:
    python analysis_script.py

Input:
    hexagene_clinvar_public.csv (Supplementary Table S2)

Outputs:
    - Terminal: T-statistics, AUC, Correlation, VUS Performance
    - hexagene_distribution.png: Figure 1
    - hexagene_vus_rescue.png: Figure 2
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

DATA_FILE = "hexagene_clinvar_public.csv"

def run_public_validation():
    if not Path(DATA_FILE).exists():
        print(f"Error: {DATA_FILE} not found.")
        print("Please ensure the dataset is in the same directory.")
        return

    df = pd.read_csv(DATA_FILE)
    
    # Clean data (drop rows with missing benchmark scores)
    df = df.dropna(subset=['HexaGene_Structural_Score', 'REVEL_Score'])
    
    path = df[df['ClinVar_Label'] == 1]['HexaGene_Structural_Score']
    ben = df[df['ClinVar_Label'] == 0]['HexaGene_Structural_Score']
    
    print("="*60)
    print("HEXAGENE VALIDATION RESULTS")
    print("="*60)
    print(f"Total Variants: {len(df)}")
    print(f"Pathogenic: {len(path)}, Benign: {len(ben)}")
    
    # 1. Overall Discrimination
    t_stat, p_val = stats.ttest_ind(path, ben, equal_var=False)
    auc = roc_auc_score(df['ClinVar_Label'], df['HexaGene_Structural_Score'])
    
    print(f"\n1. GLOBAL PERFORMANCE")
    print(f"   T-statistic: {t_stat:.2f}")
    print(f"   P-value:     {p_val:.2e}")
    print(f"   ROC AUC:     {auc:.3f}")

    # 2. Orthogonality Check
    r, _ = stats.pearsonr(df['HexaGene_Structural_Score'], df['REVEL_Score'])
    print(f"\n2. ORTHOGONALITY")
    print(f"   Correlation (r) with REVEL: {r:.3f}")
    
    # 3. VUS Grey Zone Rescue
    # Define Grey Zone: REVEL 0.4 - 0.6
    grey = df[(df['REVEL_Score'] > 0.4) & (df['REVEL_Score'] < 0.6)]
    
    print(f"\n3. GREY ZONE RESCUE (REVEL 0.4-0.6)")
    print(f"   Variants in Zone: {len(grey)}")
    
    if len(grey) > 10:
        g_path = grey[grey['ClinVar_Label'] == 1]['HexaGene_Structural_Score']
        g_ben = grey[grey['ClinVar_Label'] == 0]['HexaGene_Structural_Score']
        
        g_t, g_p = stats.ttest_ind(g_path, g_ben, equal_var=False)
        try: g_auc = roc_auc_score(grey['ClinVar_Label'], grey['HexaGene_Structural_Score'])
        except: g_auc = 0.5
        
        print(f"   HexaGene T-stat: {g_t:.2f}")
        print(f"   HexaGene AUC:    {g_auc:.3f}")
        print(f"   P-value:         {g_p:.2e}")
    else:
        print("   (Insufficient data for grey zone statistics)")

    # --- PLOTTING ---
    plot_distribution(path, ben, t_stat)
    if len(grey) > 10:
        plot_vus_rescue(grey)

def plot_distribution(path, ben, t_score):
    plt.figure(figsize=(10, 6))
    sns.kdeplot(ben, fill=True, color='green', label='Benign', alpha=0.3)
    sns.kdeplot(path, fill=True, color='red', label='Pathogenic', alpha=0.3)
    plt.title(f"Figure 1: HexaGene Score Distribution (T={t_score:.1f})")
    plt.xlabel("Structural Risk Score")
    plt.legend()
    plt.tight_layout()
    plt.savefig("Figure1_Distribution.png", dpi=300)
    print(f"\nGenerated Figure1_Distribution.png")

def plot_vus_rescue(grey_df):
    plt.figure(figsize=(8, 6))
    sns.boxplot(x='ClinVar_Label', y='HexaGene_Structural_Score', data=grey_df, palette=['green', 'red'])
    plt.xticks([0, 1], ['Benign', 'Pathogenic'])
    plt.title("Figure 2: Discrimination in REVEL Grey Zone (0.4-0.6)")
    plt.ylabel("HexaGene Score")
    plt.tight_layout()
    plt.savefig("Figure2_VUS_Rescue.png", dpi=300)
    print(f"Generated Figure2_VUS_Rescue.png")

if __name__ == "__main__":
    run_public_validation()