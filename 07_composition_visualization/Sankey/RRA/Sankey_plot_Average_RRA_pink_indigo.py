import pandas as pd
import numpy as np
import plotly.graph_objects as go
import textwrap

# =====================================================
# === LIFE FORM COLOR MAP ===
# =====================================================
lifeform_colors = {
    "Graminoid": "#228B22",   # dark green
    "Woody taxa": "#8B4513",  # brown
    "Succulent": "#90EE90",   # light green
    "Forb": "#FF69B4",        # pink
    "NA": "#B0B0B0"           # gray
}

def collapse_and_rra(df):
    """Group by Genus and convert raw counts to RRA (%) per sample."""
    sample_cols = [c for c in df.columns if c.endswith("_E") or c.endswith("_I")]
    Genus = df.groupby("Genus")[sample_cols].sum()
    totals = Genus.sum(axis=0).replace(0, np.nan)  # avoid div-by-zero
    rra = Genus.div(totals, axis=1).fillna(0) * 100.0
    rra.reset_index(inplace=True)
    return rra

def sankey_for_dataset(df, dataset_name, out_png=None, out_html=None, flip=False):
    """
    Create a Sankey diagram comparing Elephant and Impala diets at Genus level.
    Genera with RRA > 0.5% in at least one group are shown.
    Saves as PNG, HTML, and RRA CSV.
    """
    # Compute RRA
    rra = collapse_and_rra(df)

    # Map life forms
    life_form_map = df.groupby('Genus')['Life form'].first().to_dict()

    # Calculate mean RRA for Elephant and Impala
    e_cols = [c for c in rra.columns if c.endswith("_E")]
    i_cols = [c for c in rra.columns if c.endswith("_I")]
    rra["E_mean"] = rra[e_cols].mean(axis=1) if e_cols else 0
    rra["I_mean"] = rra[i_cols].mean(axis=1) if i_cols else 0

    # Filter genera
    rra = rra[(rra["E_mean"] > 0.5) | (rra["I_mean"] > 0.5)]

    # === SAVE RRA TABLE ===
    out_csv = f"{dataset_name.lower()}_rra.csv"
    rra.to_csv(out_csv, index=False)
    print(f"Saved RRA table: {out_csv}")

    # Compute ordering ratio
    total_mean = rra["E_mean"] + rra["I_mean"]
    rra["ratio"] = rra["I_mean"] / total_mean.where(total_mean != 0, 0)
    # Sort for ordering
    rra_ordered = rra.sort_values("ratio", ascending=True)
    genera = rra_ordered["Genus"].tolist()

    # Wrap Genus labels
    genera_labels = ['<br>'.join(textwrap.wrap(g, width=25)) for g in genera]

    # Compute proportional node heights
    rra_ordered["size"] = rra_ordered["E_mean"] + rra_ordered["I_mean"]
    total_size = rra_ordered["size"].sum()
    margin_y = 0.90
    min_node_height = 0.03
    y_gen = []
    cumul_y = 0.05
    if total_size > 0:
        for _, row in rra_ordered.iterrows():
            h = max((row["size"] / total_size) * margin_y, min_node_height)
            y = cumul_y + h / 2
            y_gen.append(y)
            cumul_y += h
    else:
        y_gen = np.linspace(0.05, 0.95, len(genera)).tolist()
        cumul_y = 0.95

    # Adjust height if needed
    used_span = cumul_y - 0.05
    base_height = 800
    height = base_height
    if used_span > margin_y:
        scale_factor = margin_y / used_span
        y_gen = [(y - 0.05) * scale_factor + 0.05 for y in y_gen]
        height = base_height * (used_span / margin_y)

    # Prepare nodes and links
    src, tgt, val, lcolors = [], [], [], []
    e_conns, i_conns = [], []

    # Genus colors
    Genus_colors = [lifeform_colors.get(life_form_map.get(g, "NA"), "#B0B0B0") for g in genera]

    if not flip:
        nodes = ["Elephant", "Impala"] + genera_labels
        node_map = {n: i for i, n in enumerate(nodes)}
        for gi, (_, row) in enumerate(rra_ordered.iterrows()):
            g_idx = node_map[genera_labels[gi]]
            if row["E_mean"] > 0:
                src.append(node_map["Elephant"])
                tgt.append(g_idx)
                val.append(row["E_mean"])
                lcolors.append("rgba(255,107,107,0.45)")
                e_conns.append((y_gen[gi], row["E_mean"]))
            if row["I_mean"] > 0:
                src.append(node_map["Impala"])
                tgt.append(g_idx)
                val.append(row["I_mean"])
                lcolors.append("rgba(78,205,196,0.45)")
                i_conns.append((y_gen[gi], row["I_mean"]))
        node_colors = ["#FF6B6B", "#4ECDC4"] + Genus_colors
        x = [0.01, 0.01] + [0.5] * len(genera)
        y = [] + y_gen  # Placeholder for species y
    else:
        nodes = genera_labels + ["Elephant", "Impala"]
        node_map = {n: i for i, n in enumerate(nodes)}
        for gi, (_, row) in enumerate(rra_ordered.iterrows()):
            g_idx = node_map[genera_labels[gi]]
            if row["E_mean"] > 0:
                src.append(g_idx)
                tgt.append(node_map["Elephant"])
                val.append(row["E_mean"])
                lcolors.append("rgba(255,107,107,0.45)")
                e_conns.append((y_gen[gi], row["E_mean"]))
            if row["I_mean"] > 0:
                src.append(g_idx)
                tgt.append(node_map["Impala"])
                val.append(row["I_mean"])
                lcolors.append("rgba(78,205,196,0.45)")
                i_conns.append((y_gen[gi], row["I_mean"]))
        node_colors = Genus_colors + ["#FF6B6B", "#4ECDC4"]
        x = [0.01] * len(genera) + [0.5, 0.5]
        y = y_gen + []  # Placeholder for species y

    # Compute weighted y positions for Elephant and Impala
    e_y = sum(y * v for y, v in e_conns) / sum(v for _, v in e_conns) if e_conns else 0.5
    i_y = sum(y * v for y, v in i_conns) / sum(v for _, v in i_conns) if i_conns else 0.5

    # Adjust positions with extra spacing
    source_extra_spacing = 0.1
    mid_y = (e_y + i_y) / 2
    sep = abs(e_y - i_y) + source_extra_spacing
    half_sep = sep / 2
    e_y = mid_y - half_sep
    i_y = mid_y + half_sep

    # Clamp y positions
    e_y = max(0.05, min(0.95, e_y))
    i_y = max(0.05, min(0.95, i_y))

    # Set final y positions
    if not flip:
        y = [e_y, i_y] + y_gen
    else:
        y += [e_y, i_y]

    # Create figure
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            label=nodes,
            pad=25,
            thickness=15,
            x=x,
            y=y,
            color=node_colors
        ),
        link=dict(source=src, target=tgt, value=val, color=lcolors)
    )])

    # Add life form legend
    present_lifeforms = sorted(set(life_form_map.get(g, "NA") for g in genera))
    for lf in present_lifeforms:
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=10, color=lifeform_colors[lf]),
            name=lf,
            showlegend=True
        ))

    # Update layout
    fig.update_layout(
        title=f"{dataset_name.upper()} — Elephant vs Impala (Genus-level RRA, genera >0.5% shown)",
        font=dict(size=12),
        width=1400,
        height=height,
        paper_bgcolor="white",
        plot_bgcolor="white",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=-0.1,
            xanchor='center',
            x=0.5,
            title="Life form",
            bordercolor="black",
            borderwidth=1
        )
    )

    # Save as PNG
    out_png = out_png or f"{dataset_name.lower()}_sankey.png"
    fig.write_image(out_png, width=1400, height=height, scale=2)
    print(f"Saved PNG: {out_png}")

    # Save as HTML
    out_html = out_html or f"{dataset_name.lower()}_sankey.html"
    fig.write_html(out_html, include_plotlyjs='cdn', full_html=True)
    print(f"Saved HTML: {out_html}")

# ---- Example usage ----
rbcl_df = pd.read_csv("rbcl_combined_dietlist.csv")
trnl_df = pd.read_csv("trnl_combined_dietlist.csv")

sankey_for_dataset(rbcl_df, "rbcl", out_png="rbcl_sankey.png", out_html="rbcl_sankey.html")
sankey_for_dataset(trnl_df, "trnl", out_png="trnl_sankey.png", out_html="trnl_sankey.html", flip=True)