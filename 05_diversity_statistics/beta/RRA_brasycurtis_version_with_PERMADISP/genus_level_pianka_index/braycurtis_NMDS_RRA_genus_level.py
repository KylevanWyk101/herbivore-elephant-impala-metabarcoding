# -*- coding: utf-8 -*-
"""
Created on Sat Aug 23 11:51:44 2025
Updated with Genus-level Pianka's Index calculation

@author: kylev
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova, anosim
from sklearn.manifold import MDS
from sklearn.metrics.pairwise import cosine_similarity

# ------------------ Helper functions ------------------

def load_abundance(file, id_col="ASV_ID", tax_col="Genus"):
    """
    Read CSV with ASVs as rows, taxonomy column, and *_RRA columns as samples.
    Collapses to Genus level.
    Returns:
      - abundance: DataFrame with samples as rows, genera as columns,
      - group: label inferred from filename.
    """
    df = pd.read_csv(file)
    if id_col not in df.columns or tax_col not in df.columns:
        raise ValueError(f"{id_col} or {tax_col} not found in {file}")

    # Only keep RRA columns
    rra_cols = [c for c in df.columns if c.endswith("_RRA")]
    if not rra_cols:
        raise ValueError(f"No _RRA columns found in {file}")

    # Collapse to genus level
    df[tax_col] = df[tax_col].fillna("Unassigned")
    genus_df = df.groupby(tax_col)[rra_cols].sum()

    # Transpose: samples as rows, genera as columns
    abundance = genus_df.T
    abundance.index = [c.replace("_RRA", "") for c in genus_df.columns] if "_RRA" in genus_df.columns else abundance.index

    # Assign group
    group = "Elephant" if "E_" in file else "Impala" if "I_" in file else "Unknown"

    return abundance, group


def run_stats(dm, groups, locus):
    """Run PERMANOVA and ANOSIM, save results to CSV."""
    metadata = pd.DataFrame({"Group": groups}, index=dm.ids)

    if len(set(groups)) < 2:
        raise ValueError(f"Only one group found for {locus}: {set(groups)}")

    perma_res = permanova(dm, metadata, column="Group")
    print("\nPERMANOVA:\n", perma_res)
    perma_res.to_csv(f"permanova_{locus}_genus.csv")

    anosim_res = anosim(dm, metadata, column="Group")
    print("\nANOSIM:\n", anosim_res)
    anosim_res.to_csv(f"anosim_{locus}_genus.csv")

    return perma_res, anosim_res


def plot_nmds(dm, groups, title="NMDS", locus=""):
    """Run NMDS and plot scatter with group labels."""
    nmds = MDS(n_components=2, metric=False, dissimilarity="precomputed",
               random_state=42, max_iter=3000, n_init=4)
    coords = nmds.fit_transform(dm.data)
    scores = pd.DataFrame(coords, columns=["NMDS1", "NMDS2"], index=dm.ids)
    scores["Group"] = groups.values

    plt.figure(figsize=(6, 6))
    for grp, sub in scores.groupby("Group"):
        plt.scatter(sub["NMDS1"], sub["NMDS2"], label=grp, s=60, alpha=0.8)

    plt.xlabel("NMDS1")
    plt.ylabel("NMDS2")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"nmds_{locus}_genus.png", dpi=300)
    plt.close()
    print(f"NMDS stress for {locus}: {nmds.stress_:.4f}")


# ------------------ Pianka’s Index functions ------------------

def pairwise_pianka(df, normalize=True):
    """
    Compute pairwise Pianka overlaps among all samples in df.
    df : DataFrame, rows = samples, columns = Genera/features.
    """
    X = df.copy().astype(float)

    if normalize:
        X = X.div(X.sum(axis=1).replace(0, 1), axis=0)  # row proportions

    # Cosine similarity is identical to Pianka’s index
    mat = cosine_similarity(X.values)
    return pd.DataFrame(mat, index=X.index, columns=X.index)


def summarize_pianka(pianka_mat, groups):
    """Average Pianka index by group pairs."""
    pairs = []
    for i, si in enumerate(pianka_mat.index):
        for j, sj in enumerate(pianka_mat.index):
            if j > i:
                g1, g2 = groups[si], groups[sj]
                pairs.append((g1, g2, pianka_mat.loc[si, sj]))

    summary = (
        pd.DataFrame(pairs, columns=["Group1", "Group2", "Pianka"])
        .groupby(["Group1", "Group2"])
        .mean()
    )
    return summary


# ------------------ Main loop ------------------

datasets = {
    "rbcL": ["E_rbcl_diet_rra_only.csv", "I_rbcl_diet_rra_only.csv"],
    "trnL": ["E_trnl_diet_rra_only.csv", "I_trnl_diet_rra_only.csv"]
}

for locus, files in datasets.items():
    print(f"\n===== Processing locus: {locus} (Genus-level, Bray-Curtis + Pianka) =====")

    # Load abundance tables (already collapsed to Genus level)
    abundance_dfs = [load_abundance(f) for f in files]

    # Combine sample tables
    combined_df = pd.concat([df[0] for df in abundance_dfs], axis=0).fillna(0)

    # Build group labels
    groups = pd.Series(index=combined_df.index, dtype=str)
    for df, group in abundance_dfs:
        groups[df.index] = group

    print(f"Groups for {locus}:\n", groups)

    # ------------------ Bray-Curtis & ordination ------------------
    dm = beta_diversity("braycurtis", combined_df.values, ids=combined_df.index)

    # Run PERMANOVA & ANOSIM
    run_stats(dm, groups, locus)

    # NMDS plot
    plot_nmds(dm, groups, title=f"NMDS - {locus} (Bray-Curtis, Genus-level)", locus=locus)

    # ------------------ Pianka’s Index ------------------
    pianka_mat = pairwise_pianka(combined_df, normalize=True)
    pianka_mat.to_csv(f"pianka_matrix_{locus}_genus.csv")
    print(f"Pairwise Pianka’s index matrix for {locus} saved.")

    pianka_summary = summarize_pianka(pianka_mat, groups)
    pianka_summary.to_csv(f"pianka_summary_{locus}_genus.csv")
    print(f"Group-level Pianka summary for {locus} saved.")

print("\nGenus-level beta diversity + Pianka’s index analysis complete for all loci.")
