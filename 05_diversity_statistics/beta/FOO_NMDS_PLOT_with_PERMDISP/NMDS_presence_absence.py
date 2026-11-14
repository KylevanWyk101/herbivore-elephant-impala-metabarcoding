import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova, anosim, permdisp, DistanceMatrix
from sklearn.manifold import MDS


# ------------------ Settings ------------------

SHOW_SAMPLE_LABELS = True      # Change to False if labels are too messy
ADD_JITTER = True              # Adds tiny random movement to reveal overlapping points
JITTER_AMOUNT = 0.015          # Increase slightly if points still overlap
POINT_SIZE = 35
POINT_ALPHA = 0.55


# ------------------ Helper functions ------------------

def load_abundance(file, id_col="ASV_ID"):
    """Load CSV, extract _RRA sample columns, transpose so samples are rows."""
    df = pd.read_csv(file)

    if id_col not in df.columns:
        raise ValueError(f"{id_col} not found in {file}")

    df = df.set_index(id_col)

    rra_cols = [c for c in df.columns if c.endswith("_RRA")]

    if not rra_cols:
        raise ValueError(f"No _RRA columns found in {file}")

    abundance = df[rra_cols].T
    abundance.index = [c.replace("_RRA", "") for c in rra_cols]

    if "E_" in file:
        group = "Elephant"
    elif "I_" in file:
        group = "Impala"
    else:
        group = "Unknown"

    print("\n----------------------------------------")
    print(f"Loaded file: {file}")
    print(f"Detected group: {group}")
    print(f"Number of _RRA sample columns found: {len(rra_cols)}")
    print("Sample IDs:")
    for sample in abundance.index:
        print(f"  {sample}")

    return abundance, group


def presence_absence(abundance_df):
    """Convert abundance table to presence/absence."""
    return (abundance_df > 0).astype(int)


def run_stats(dm, groups, locus):
    """Run PERMANOVA, ANOSIM, PERMDISP, and sqrt-transformed PERMDISP."""
    metadata = pd.DataFrame({"Group": groups}, index=dm.ids)

    if len(set(groups)) < 2:
        raise ValueError(
            f"Only one group found for {locus}: {set(groups)}. Need at least two groups."
        )

    perma_res = permanova(dm, metadata, column="Group", permutations=999)
    print("\nPERMANOVA:\n", perma_res)
    perma_res.to_csv(f"permanova_{locus}.csv")

    anosim_res = anosim(dm, metadata, column="Group", permutations=999)
    print("\nANOSIM:\n", anosim_res)
    anosim_res.to_csv(f"anosim_{locus}.csv")

    permdisp_res = permdisp(dm, metadata, column="Group", permutations=999)
    print("\nPERMDISP - original Jaccard:\n", permdisp_res)
    permdisp_res.to_csv(f"permdisp_{locus}.csv")

    sqrt_dm = DistanceMatrix(
        np.sqrt(dm.data),
        ids=dm.ids
    )

    permdisp_sqrt_res = permdisp(sqrt_dm, metadata, column="Group", permutations=999)
    print("\nPERMDISP - sqrt-transformed Jaccard:\n", permdisp_sqrt_res)
    permdisp_sqrt_res.to_csv(f"permdisp_sqrt_jaccard_{locus}.csv")

    return perma_res, anosim_res, permdisp_res, permdisp_sqrt_res


def detect_overlapping_points(scores, locus, decimals=3):
    """
    Check whether NMDS points overlap or nearly overlap.
    Rounds coordinates to detect points that occupy similar positions.
    """
    tmp = scores.copy()
    tmp["NMDS1_round"] = tmp["NMDS1"].round(decimals)
    tmp["NMDS2_round"] = tmp["NMDS2"].round(decimals)

    overlap_summary = (
        tmp.groupby(["NMDS1_round", "NMDS2_round"])
        .size()
        .reset_index(name="n_points")
        .query("n_points > 1")
        .sort_values("n_points", ascending=False)
    )

    overlap_summary.to_csv(f"overlapping_points_{locus}.csv", index=False)

    print(f"\nPossible overlapping / near-overlapping NMDS points for {locus}:")
    if overlap_summary.empty:
        print("No obvious overlapping points detected at this rounding level.")
    else:
        print(overlap_summary)

    return overlap_summary


def plot_nmds(dm, groups, title="NMDS", locus=""):
    """Run NMDS and plot with jitter and labels to reveal overlapping dots."""
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

    # Save original NMDS coordinates
    scores.to_csv(f"nmds_scores_{locus}.csv")

    print(f"\nNMDS coordinate table saved: nmds_scores_{locus}.csv")
    print(f"Number of points in NMDS for {locus}: {scores.shape[0]}")
    print(scores["Group"].value_counts())

    # Detect overlapping or near-overlapping points
    detect_overlapping_points(scores, locus, decimals=3)

    # Create plotting coordinates
    plot_scores = scores.copy()

    if ADD_JITTER:
        rng = np.random.default_rng(42)
        plot_scores["NMDS1_plot"] = plot_scores["NMDS1"] + rng.normal(
            0, JITTER_AMOUNT, size=len(plot_scores)
        )
        plot_scores["NMDS2_plot"] = plot_scores["NMDS2"] + rng.normal(
            0, JITTER_AMOUNT, size=len(plot_scores)
        )
    else:
        plot_scores["NMDS1_plot"] = plot_scores["NMDS1"]
        plot_scores["NMDS2_plot"] = plot_scores["NMDS2"]

    fig, ax = plt.subplots(figsize=(7, 7))

    for grp, sub in plot_scores.groupby("Group"):
        ax.scatter(
            sub["NMDS1_plot"],
            sub["NMDS2_plot"],
            label=f"{grp} (n={len(sub)})",
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            edgecolors="black",
            linewidths=0.4
        )

    if SHOW_SAMPLE_LABELS:
        for sample_id, row in plot_scores.iterrows():
            ax.text(
                row["NMDS1_plot"],
                row["NMDS2_plot"],
                sample_id,
                fontsize=6,
                alpha=0.75
            )

    ax.set_xlabel("NMDS1", fontsize=12, color="black")
    ax.set_ylabel("NMDS2", fontsize=12, color="black")

    if ADD_JITTER:
        ax.set_title(f"{title}\nJitter added only for visualisation", fontsize=13, color="black")
    else:
        ax.set_title(title, fontsize=14, color="black")

    ax.grid(False)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1.0)

    ax.tick_params(
        axis="both",
        colors="black",
        width=1.0,
        labelsize=10
    )

    ax.legend(frameon=True, fontsize=10)

    plt.tight_layout()

    plt.savefig(f"nmds_{locus}_jitter_labels.png", dpi=600, bbox_inches="tight")
    plt.savefig(f"nmds_{locus}_jitter_labels.pdf", bbox_inches="tight")

    plt.close()

    print(f"NMDS stress for {locus}: {nmds.stress_:.4f}")
    print(f"Saved jittered NMDS plot: nmds_{locus}_jitter_labels.png")
    print(f"Saved jittered NMDS plot: nmds_{locus}_jitter_labels.pdf")


# ------------------ Main loop ------------------

datasets = {
    "rbcL": ["E_rbcl_diet_rra_only.csv", "I_rbcl_diet_rra_only.csv"],
    "trnL": ["E_trnl_diet_rra_only.csv", "I_trnl_diet_rra_only.csv"]
}

for locus, files in datasets.items():
    print("\n========================================")
    print(f"Processing locus: {locus}")
    print("========================================")

    abundance_dfs = [load_abundance(f) for f in files]

    combined_df = pd.concat(
        [df_group[0] for df_group in abundance_dfs],
        axis=0
    )

    combined_df = combined_df.fillna(0)

    pa_df = presence_absence(combined_df)

    groups = pd.Series(index=pa_df.index, dtype=str)

    for abundance_df, group in abundance_dfs:
        groups.loc[abundance_df.index] = group

    print("\n----------------------------------------")
    print(f"FINAL SAMPLE COUNT USED FOR {locus}")
    print("----------------------------------------")
    print(groups.value_counts())
    print(f"Total samples used for {locus}: {len(groups)}")

    sample_summary = pd.DataFrame({
        "Sample_ID": groups.index,
        "Group": groups.values
    })

    sample_summary.to_csv(f"sample_counts_used_{locus}.csv", index=False)
    print(f"Saved sample list: sample_counts_used_{locus}.csv")

    zero_samples = pa_df.index[pa_df.sum(axis=1) == 0].tolist()

    if zero_samples:
        print("\nWARNING: These samples have zero detected taxa after presence/absence conversion:")
        for s in zero_samples:
            print(f"  {s}")
    else:
        print("\nNo all-zero samples detected.")

    dm = beta_diversity(
        "jaccard",
        pa_df.values,
        ids=pa_df.index
    )

    perma_res, anosim_res, permdisp_res, permdisp_sqrt_res = run_stats(
        dm,
        groups,
        locus
    )

    plot_nmds(
        dm,
        groups,
        title=f"NMDS - {locus} (Jaccard, Presence/Absence)",
        locus=locus
    )

print("\nBeta diversity analysis complete for all loci.")