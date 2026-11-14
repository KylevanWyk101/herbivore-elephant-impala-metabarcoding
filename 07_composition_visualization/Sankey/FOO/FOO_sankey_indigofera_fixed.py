import pandas as pd
import numpy as np
import plotly.graph_objects as go
import textwrap

# === LIFE FORM COLOR MAP ===
lifeform_colors = {
    "Graminoid": "#228B22",
    "Woody taxa": "#8B4513",
    "Succulent": "#90EE90",
    "Forb": "#FF69B4",
    "NA": "#B0B0B0"
}

# === COLLAPSE ASVs TO GENUS, THEN CALCULATE FOO% ===
def collapse_and_foo(df):
    """
    Collapse ASV-level counts to genus level by SUM, then convert genus x sample counts
    to presence/absence.

    Returns a genus-level DataFrame with:
      - Genus
      - Life form
      - sample-level presence/absence columns
    """

    # Detect sample columns
    sample_cols = [c for c in df.columns if c.endswith("_E") or c.endswith("_I")]

    if not sample_cols:
        raise ValueError("No sample columns found ending in _E or _I.")

    # Ensure Genus column exists
    if "Genus" not in df.columns:
        raise ValueError("Input dataframe must contain a 'Genus' column.")

    # Collapse ASVs to genus by summing counts across ASVs
    grouped = df.groupby("Genus")[sample_cols].sum().reset_index()

    # Add life form information
    if "Life form" in df.columns:
        life_form = df.groupby("Genus")["Life form"].first().reset_index()
        grouped = pd.merge(grouped, life_form, on="Genus", how="left")
    else:
        grouped["Life form"] = "NA"

    # Convert genus-level counts to presence/absence per sample
    for c in sample_cols:
        grouped[c] = (grouped[c] > 0).astype(int)

    return grouped


# === SANKEY DIAGRAM FOR DATASET ===
def sankey_for_dataset(df, dataset_name, out_png=None, out_html=None, flip=False):
    """
    Create Sankey diagram comparing Elephant and Impala diets at genus level.

    Important:
    - ASVs are first collapsed to genus level.
    - Genus-level sample columns are converted to presence/absence.
    - FOO is calculated as percentage of samples in which each genus occurs:
        FOO% = number of samples where genus is present / total samples in group * 100
    """

    # ---- Compute genus-level presence/absence table ----
    genus_pa = collapse_and_foo(df)

    # Map each genus to its life form
    life_form_map = genus_pa.set_index("Genus")["Life form"].to_dict()

    # Identify elephant and impala sample columns
    e_cols = [c for c in genus_pa.columns if c.endswith("_E")]
    i_cols = [c for c in genus_pa.columns if c.endswith("_I")]

    if not e_cols and not i_cols:
        raise ValueError("No sample columns found ending in _E or _I.")

    # Number of samples per herbivore group
    n_e = len(e_cols)
    n_i = len(i_cols)

    print(f"{dataset_name.upper()} elephant sample count: {n_e}")
    print(f"{dataset_name.upper()} impala sample count: {n_i}")

    # Calculate FOO counts
    genus_pa["E_FOO_count"] = genus_pa[e_cols].sum(axis=1) if n_e > 0 else 0
    genus_pa["I_FOO_count"] = genus_pa[i_cols].sum(axis=1) if n_i > 0 else 0

    # Calculate FOO percentages
    genus_pa["E_FOO"] = (genus_pa["E_FOO_count"] / n_e) * 100 if n_e > 0 else 0
    genus_pa["I_FOO"] = (genus_pa["I_FOO_count"] / n_i) * 100 if n_i > 0 else 0

    # Keep only genera that occur in at least one sample
    foo = genus_pa[
        (genus_pa["E_FOO_count"] > 0) | (genus_pa["I_FOO_count"] > 0)
    ].copy()

    # Save FOO CSV with both counts and percentages
    out_csv = f"{dataset_name.lower()}_foo_percent.csv"
    foo.to_csv(out_csv, index=False)
    print(f"Saved FOO% table: {out_csv}")

    # Ordering for layout
    total_foo = foo["E_FOO"] + foo["I_FOO"]

    # Ratio based on relative contribution to impala FOO%
    foo["ratio"] = np.where(total_foo > 0, foo["I_FOO"] / total_foo, 0.0)

    # Sort from elephant-associated to impala-associated
    foo_ordered = foo.sort_values("ratio", ascending=True).reset_index(drop=True)

    genera = foo_ordered["Genus"].tolist()

    # Wrap genus labels
    genera_labels = ['<br>'.join(textwrap.wrap(str(g), width=25)) for g in genera]

    # Node heights based on total FOO%
    foo_ordered["size"] = foo_ordered["E_FOO"] + foo_ordered["I_FOO"]
    total_size = foo_ordered["size"].sum()

    margin_y = 0.90
    min_node_height = 0.03
    y_gen = []
    cumul_y = 0.05

    if total_size > 0:
        for _, row in foo_ordered.iterrows():
            h = max((row["size"] / total_size) * margin_y, min_node_height)
            y = cumul_y + h / 2
            y_gen.append(y)
            cumul_y += h
    else:
        y_gen = np.linspace(0.05, 0.95, max(1, len(genera))).tolist()
        cumul_y = 0.95

    used_span = cumul_y - 0.05
    base_height = 800
    height = base_height

    if used_span > margin_y:
        scale_factor = margin_y / used_span
        y_gen = [(y - 0.05) * scale_factor + 0.05 for y in y_gen]
        height = int(base_height * (used_span / margin_y))

    # Prepare Sankey nodes and links
    src, tgt, val, lcolors = [], [], [], []
    e_conns, i_conns = [], []

    genus_colors = [
        lifeform_colors.get(life_form_map.get(g, "NA"), "#B0B0B0")
        for g in genera
    ]

    # Build node list depending on flip
    if not flip:
        nodes = ["Elephant", "Impala"] + genera_labels
        node_map = {n: i for i, n in enumerate(nodes)}

        for gi, row in foo_ordered.iterrows():
            g_idx = node_map[genera_labels[gi]]

            if row["E_FOO"] > 0:
                src.append(node_map["Elephant"])
                tgt.append(g_idx)
                val.append(float(row["E_FOO"]))
                lcolors.append("rgba(255,107,107,0.45)")
                e_conns.append((y_gen[gi], float(row["E_FOO"])))

            if row["I_FOO"] > 0:
                src.append(node_map["Impala"])
                tgt.append(g_idx)
                val.append(float(row["I_FOO"]))
                lcolors.append("rgba(78,205,196,0.45)")
                i_conns.append((y_gen[gi], float(row["I_FOO"])))

        node_colors = ["#FF6B6B", "#4ECDC4"] + genus_colors
        x = [0.01, 0.01] + [0.5] * len(genera)

    else:
        nodes = genera_labels + ["Elephant", "Impala"]
        node_map = {n: i for i, n in enumerate(nodes)}

        for gi, row in foo_ordered.iterrows():
            g_idx = node_map[genera_labels[gi]]

            if row["E_FOO"] > 0:
                src.append(g_idx)
                tgt.append(node_map["Elephant"])
                val.append(float(row["E_FOO"]))
                lcolors.append("rgba(255,107,107,0.45)")
                e_conns.append((y_gen[gi], float(row["E_FOO"])))

            if row["I_FOO"] > 0:
                src.append(g_idx)
                tgt.append(node_map["Impala"])
                val.append(float(row["I_FOO"]))
                lcolors.append("rgba(78,205,196,0.45)")
                i_conns.append((y_gen[gi], float(row["I_FOO"])))

        node_colors = genus_colors + ["#FF6B6B", "#4ECDC4"]
        x = [0.01] * len(genera) + [0.5, 0.5]

    # Center Elephant and Impala nodes vertically based on weighted link positions
    e_y = sum(y * v for y, v in e_conns) / sum(v for _, v in e_conns) if e_conns else 0.5
    i_y = sum(y * v for y, v in i_conns) / sum(v for _, v in i_conns) if i_conns else 0.5

    source_extra_spacing = 0.1
    mid_y = (e_y + i_y) / 2
    sep = abs(e_y - i_y) + source_extra_spacing
    half_sep = sep / 2

    e_y = mid_y - half_sep
    i_y = mid_y + half_sep

    e_y = max(0.05, min(0.95, e_y))
    i_y = max(0.05, min(0.95, i_y))

    # Final y positions
    if not flip:
        y = [e_y, i_y] + y_gen
    else:
        y = y_gen + [e_y, i_y]

    # Create Sankey figure
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            label=nodes,
            pad=25,
            thickness=15,
            x=x,
            y=y,
            color=node_colors
        ),
        link=dict(
            source=src,
            target=tgt,
            value=val,
            color=lcolors
        )
    )])

    # Life form legend
    present_lifeforms = sorted(set(life_form_map.get(g, "NA") for g in genera))

    for lf in present_lifeforms:
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode="markers",
            marker=dict(
                size=10,
                color=lifeform_colors.get(lf, "#B0B0B0")
            ),
            name=lf,
            showlegend=True
        ))

    fig.update_layout(
        title=f"{dataset_name.upper()} — Elephant vs Impala (Genus-level FOO%)",
        font=dict(size=12),
        width=1400,
        height=height,
        paper_bgcolor="white",
        plot_bgcolor="white",
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=-0.1,
            xanchor="center",
            x=0.5,
            title="Life form",
            bordercolor="black",
            borderwidth=1
        )
    )

    # Save files
    out_png = out_png or f"{dataset_name.lower()}_sankey_foo_percent.png"
    fig.write_image(out_png, width=1400, height=height, scale=2)
    print(f"Saved PNG: {out_png}")

    out_html = out_html or f"{dataset_name.lower()}_sankey_foo_percent.html"
    fig.write_html(out_html, include_plotlyjs="cdn", full_html=True)
    print(f"Saved HTML: {out_html}")


# === Example usage ===
rbcl_df = pd.read_csv("rbcl_combined_dietlist.csv")
trnl_df = pd.read_csv("trnl_combined_dietlist.csv")

sankey_for_dataset(
    rbcl_df,
    "rbcl",
    out_png="rbcl_sankey_foo_percent.png",
    out_html="rbcl_sankey_foo_percent.html"
)

sankey_for_dataset(
    trnl_df,
    "trnl",
    out_png="trnl_sankey_foo_percent.png",
    out_html="trnl_sankey_foo_percent.html",
    flip=True
)