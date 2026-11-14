# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 11:04:37 2025

@author: kylev
"""

from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# --------------------------------------
# Load and prepare data
# --------------------------------------
df = pd.read_csv("combined_shannon_index_genus.csv")

# Clean string columns
df['Group'] = df['Group'].str.strip().str.capitalize()
df['Locus'] = df['Locus'].str.strip().str.lower()

# Create clean display labels for plotting
locus_labels = {
    'rbcl': r'$\it{rbc}$L',
    'trnl': r'$\it{trn}$L'
}
df['Locus_display'] = df['Locus'].map(locus_labels).fillna(df['Locus'])

# Store test results
results = []

# --------------------------------------
# Visualize distribution of data (histogram + KDE)
# --------------------------------------
unique_loci = df['Locus'].unique()
n_loci = len(unique_loci)
fig, axes = plt.subplots(n_loci, 1, figsize=(8, 4 * n_loci))

if n_loci == 1:
    axes = [axes]  # ensure iterable

for i, locus in enumerate(unique_loci):
    sub_df = df[df['Locus'] == locus]
    ax = axes[i]

    locus_name = locus_labels.get(locus, locus)

    sns.histplot(data=sub_df, x='ShannonIndex', hue='Group', kde=True,
                 palette='Set2', element='step', stat='density',
                 common_norm=False, ax=ax)

    ax.set_title(f"Distribution of Shannon Index - {locus_name}")
    ax.set_xlabel("Shannon Index")
    ax.set_ylabel("Density")

plt.tight_layout()

# Save distribution plot
plt.savefig("Shannon_distribution_histogram_KDE.pdf",
            format="pdf", bbox_inches="tight")
plt.savefig("Shannon_distribution_histogram_KDE.png",
            format="png", dpi=600, bbox_inches="tight")

plt.show()

# --------------------------------------
# Mann–Whitney U Test per locus
# --------------------------------------
for locus in unique_loci:
    sub_df = df[df['Locus'] == locus]
    elephant = sub_df[sub_df['Group'] == 'Elephant']['ShannonIndex']
    impala = sub_df[sub_df['Group'] == 'Impala']['ShannonIndex']

    stat, p_val = mannwhitneyu(elephant, impala, alternative='two-sided')

    results.append({
        'Locus': locus,
        'Test': "Mann–Whitney U",
        'Statistic': stat,
        'p-value': p_val
    })

    locus_name = locus_labels.get(locus, locus)

    print(f"\n--- {locus_name} ---")
    print(f"Mann–Whitney U test: stat = {stat:.3f}, p = {p_val:.4f}")

# --------------------------------------
# Boxplot + Stripplot with p-values
# --------------------------------------
sns.set(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(12, 6))

ax = sns.boxplot(data=df, x='Locus_display', y='ShannonIndex', hue='Group',
                 palette='pastel', showfliers=True)

sns.stripplot(data=df, x='Locus_display', y='ShannonIndex', hue='Group',
              dodge=True, marker='o', alpha=0.7,
              linewidth=0.5, palette='dark', edgecolor='gray')

# Remove duplicate legend entries
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], title="Group")

# Add p-value annotations
for i, locus in enumerate(unique_loci):
    p_val = results[i]['p-value']
    label = f"p={p_val:.3f}" if p_val >= 0.05 else "p<0.05"

    locus_display = locus_labels.get(locus, locus)
    group_max = df[df['Locus'] == locus]['ShannonIndex'].max()
    y_offset = 0.05 * (df['ShannonIndex'].max() - df['ShannonIndex'].min())
    y_pos = group_max + y_offset

    ax.text(i, y_pos, label, ha='center', va='bottom',
            fontsize=10, color='black')

# Final labels
plt.ylabel("Shannon Diversity Index")
plt.xlabel("Locus")
plt.tight_layout()

# Save final alpha diversity figure
plt.savefig("Figure_1_alpha_diversity.pdf",
            format="pdf", bbox_inches="tight")
plt.savefig("Figure_1_alpha_diversity.png",
            format="png", dpi=600, bbox_inches="tight")

plt.show()