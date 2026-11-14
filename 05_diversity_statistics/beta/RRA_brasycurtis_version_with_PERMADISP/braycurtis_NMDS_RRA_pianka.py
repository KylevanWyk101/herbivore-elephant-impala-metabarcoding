import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova, anosim, permdisp
from skbio import DistanceMatrix

from sklearn.manifold import MDS
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import mannwhitneyu

# =====================================================
# === Helper functions
# =====================================================

def load_abundance(file, id_col="ASV_ID"):
    """
    Read CSV with ASVs as rows and *_RRA columns as samples.
    Returns:
      - abundance: DataFrame with samples as rows, ASVs as columns
      - group: label inferred from filename
    """
    df = pd.read_csv(file)

    if id_col not in df.columns:
        raise ValueError(f"{id_col} not found in {file}")

    df = df.set_index(id_col)

    # Only keep RRA columns
    rra_cols = [c for c in df.columns if c.endswith("_RRA")]

    if not rra_cols:
        raise ValueError(f"No _RRA columns found in {file}")

    # Transpose: samples as rows, ASVs as columns
    abundance = df[rra_cols].T
    abundance.index = [c.replace("_RRA", "") for c in rra_cols]

    # Assign group
    if "E_" in file:
        group = "Elephant"
    elif "I_" in file:
        group = "Impala"
    else:
        group = "Unknown"

    return abundance, group


def run_stats(dm, groups, locus):
    """
    Run PERMANOVA, ANOSIM, PERMDISP, and sqrt-transformed PERMDISP.
    
    PERMDISP tests whether the groups differ in multivariate dispersion.
    This is important because PERMANOVA can be sensitive to heterogeneity
    of dispersion between groups.
    """
    metadata = pd.DataFrame({"Group": groups}, index=dm.ids)

    print("\nMetadata Table:\n", metadata)

    if len(set(groups)) < 2:
        raise ValueError(f"Only one group found for {locus}: {set(groups)}")

    # PERMANOVA
    # Tests whether elephant and impala diet composition differs significantly.
    perma_res = permanova(dm, metadata, column="Group", permutations=999)
    print("\nPERMANOVA:\n", perma_res)
    perma_res.to_csv(f"permanova_{locus}_rra.csv")

    # ANOSIM
    # Tests whether between-group dissimilarities are greater than within-group dissimilarities.
    anosim_res = anosim(dm, metadata, column="Group", permutations=999)
    print("\nANOSIM:\n", anosim_res)
    anosim_res.to_csv(f"anosim_{locus}_rra.csv")

    # PERMDISP on original Bray-Curtis distances
    # Tests whether elephant and impala groups differ in within-group dispersion.
    permdisp_res = permdisp(dm, metadata, column="Group", permutations=999)
    print("\nPERMDISP - original Bray-Curtis:\n", permdisp_res)
    permdisp_res.to_csv(f"permdisp_{locus}_rra.csv")

    # PERMDISP on sqrt-transformed Bray-Curtis distances
    # This is a sensitivity check in case the PCoA step produces negative eigenvalue warnings.
    sqrt_dm = DistanceMatrix(
        np.sqrt(dm.data),
        ids=dm.ids
    )

    permdisp_sqrt_res = permdisp(
        sqrt_dm,
        metadata,
        column="Group",
        permutations=999
    )

    print("\nPERMDISP - sqrt-transformed Bray-Curtis:\n", permdisp_sqrt_res)
    permdisp_sqrt_res.to_csv(f"permdisp_sqrt_braycurtis_{locus}_rra.csv")

    return perma_res, anosim_res, permdisp_res, permdisp_sqrt_res


def plot_nmds(dm, groups, title="NMDS", locus=""):
    """Run NMDS and plot scatter with group labels."""
    nmds = MDS(
        n_components=2,
        metric=False,
        dissimilarity="precomputed",
        random_state=42,
        max_iter=3000,
        n_init=4
    )

    coords = nmds.fit_transform(dm.data)

    scores = pd.DataFrame(
        coords,
        columns=["NMDS1", "NMDS2"],
        index=dm.ids
    )

    scores["Group"] = groups.values

    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot by group, keeping coloured groups
    for grp, sub in scores.groupby("Group"):
        ax.scatter(
            sub["NMDS1"],
            sub["NMDS2"],
            label=grp,
            s=40,
            alpha=0.85
        )

    # Labels and title
    ax.set_xlabel("NMDS1", fontsize=12, color="black")
    ax.set_ylabel("NMDS2", fontsize=12, color="black")
    ax.set_title(title, fontsize=14, color="black")

    # Remove grid
    ax.grid(False)

    # Make axis/spine lines black
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)

    # Make ticks and tick labels black
    ax.tick_params(
        axis="both",
        colors="black",
        width=1.0,
        labelsize=10
    )

    # Legend
    ax.legend(frameon=True, fontsize=10)

    plt.tight_layout()

    # Save high-resolution files
    plt.savefig(f"nmds_{locus}_rra.png", dpi=600, bbox_inches="tight")
    plt.savefig(f"nmds_{locus}_rra.pdf", bbox_inches="tight")

    plt.close()

    print(f"NMDS stress for {locus}: {nmds.stress_:.4f}")


# =====================================================
# === Pianka’s Index functions
# =====================================================

def pairwise_pianka(df, normalize=True):
    """
    Compute pairwise Pianka overlaps among all samples in df.
    df: DataFrame, rows = samples, columns = ASVs/features.
    """
    X = df.copy().astype(float)

    # This normalizes the data to RRA per sample if it is not already
    if normalize:
        X = X.div(X.sum(axis=1).replace(0, 1), axis=0)

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


def test_pianka_groups(pianka_mat, groups, locus):
    """Test whether Elephant–Impala overlaps differ from within-group overlaps."""
    pairs = []

    for i, si in enumerate(pianka_mat.index):
        for j, sj in enumerate(pianka_mat.index):
            if j > i:
                g1, g2 = groups[si], groups[sj]
                pairs.append((g1, g2, pianka_mat.loc[si, sj]))

    df_pairs = pd.DataFrame(pairs, columns=["Group1", "Group2", "Pianka"])

    # Between-group: Elephant–Impala
    between = df_pairs[df_pairs["Group1"] != df_pairs["Group2"]]["Pianka"]
    print("Between-group pairs:", len(between))

    # Within-group: Elephant–Elephant + Impala–Impala
    within = df_pairs[df_pairs["Group1"] == df_pairs["Group2"]]["Pianka"]
    print("Within-group pairs:", len(within))

    # Mann–Whitney U test
    stat, pval = mannwhitneyu(within, between, alternative="greater")

    print(f"\nPianka index test ({locus}):")
    print(f"  Mean within-group overlap = {within.mean():.3f}")
    print(f"  Mean between-group overlap = {between.mean():.3f}")
    print(f"  Mann–Whitney U stat = {stat:.3f}, p = {pval:.3e}")

    with open(f"pianka_test_{locus}_rra.txt", "w") as f:
        f.write(f"Pianka index test ({locus})\n")
        f.write(f"Mean within-group overlap = {within.mean():.3f}\n")
        f.write(f"Mean between-group overlap = {between.mean():.3f}\n")
        f.write(f"Mann–Whitney U stat = {stat:.3f}, p = {pval:.4f}\n")


def plot_pianka_heatmap(pianka_mat, groups, locus):
    """Plot heatmap of Pianka index values with group annotations."""
    # Reorder samples by group for clearer blocks
    sorted_ids = sorted(pianka_mat.index, key=lambda x: groups[x])
    mat_sorted = pianka_mat.loc[sorted_ids, sorted_ids]

    fig, ax = plt.subplots(figsize=(8, 6))

    sns.heatmap(
        mat_sorted,
        cmap="viridis",
        vmin=0,
        vmax=1,
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "Pianka index"},
        ax=ax
    )

    ax.set_title(f"Pianka index heatmap - {locus}", fontsize=14, color="black")

    # Make outer plot lines black
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)

    ax.tick_params(
        axis="both",
        colors="black",
        labelsize=8
    )

    plt.tight_layout()

    # Save high-resolution files
    plt.savefig(f"pianka_heatmap_{locus}_rra.png", dpi=600, bbox_inches="tight")
    plt.savefig(f"pianka_heatmap_{locus}_rra.pdf", bbox_inches="tight")

    plt.close()

    print(f"Pianka heatmap for {locus} saved.")


# =====================================================
# === Main loop
# =====================================================

datasets = {
    "rbcL": ["E_rbcl_diet_rra_only.csv", "I_rbcl_diet_rra_only.csv"],
    "trnL": ["E_trnl_diet_rra_only.csv", "I_trnl_diet_rra_only.csv"]
}

for locus, files in datasets.items():
    print(f"\n===== Processing locus: {locus} (RRA, Bray-Curtis + Pianka) =====")

    # Load abundance tables
    abundance_dfs = [load_abundance(f) for f in files]

    # Combine sample tables
    combined_df = pd.concat(
        [df_group[0] for df_group in abundance_dfs],
        axis=0
    ).fillna(0)

    # Build group labels
    groups = pd.Series(index=combined_df.index, dtype=str)

    for abundance_df, group in abundance_dfs:
        groups[abundance_df.index] = group

    print(f"Groups for {locus}:\n", groups)

    # ------------------ Bray-Curtis & ordination ------------------
    dm = beta_diversity(
        "braycurtis",
        combined_df.values,
        ids=combined_df.index
    )

    # Run PERMANOVA, ANOSIM, PERMDISP, and sqrt-transformed PERMDISP
    perma_res, anosim_res, permdisp_res, permdisp_sqrt_res = run_stats(
        dm,
        groups,
        locus
    )

    # NMDS plot
    plot_nmds(
        dm,
        groups,
        title=f"NMDS - {locus} (Bray-Curtis, RRA)",
        locus=locus
    )

    # ------------------ Pianka’s Index ------------------
    pianka_mat = pairwise_pianka(combined_df, normalize=True)
    pianka_mat.to_csv(f"pianka_matrix_{locus}_rra.csv")
    print(f"Pairwise Pianka’s index matrix for {locus} saved.")

    pianka_summary = summarize_pianka(pianka_mat, groups)
    pianka_summary.to_csv(f"pianka_summary_{locus}_rra.csv")
    print(f"Group-level Pianka summary for {locus} saved.")

    # Test for significance
    test_pianka_groups(pianka_mat, groups, locus)

    # Heatmap
    plot_pianka_heatmap(pianka_mat, groups, locus)

print("\nBeta diversity + Pianka’s index analysis complete for all loci.")