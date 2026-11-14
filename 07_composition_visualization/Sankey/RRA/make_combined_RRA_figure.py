"""
Combined rbcL + trnL RRA Sankey figure, exported directly as a print-ready
PNG at 1000 DPI (top of Wiley's 600-1000 DPI "line art" preferred range,
electronic_artwork_guidelines.pdf).

Replaces the old workflow of exporting each gene's Sankey to PNG, manually
screenshotting/arranging them plus silhouettes and a legend in PowerPoint,
then exporting to PDF. The node/link/RRA logic below matches
Sankey_plot_Average_RRA_pink_indigo.py exactly; that script and its outputs
are untouched by this one.

Run from this directory:
    python make_combined_RRA_figure.py
"""
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from PIL import Image, ImageFont

HERE = Path(__file__).resolve().parent
ASSETS = HERE.parent / "assets"

TARGET_DPI = 1000         # top of Wiley's 600-1000 DPI "line art" preferred range (electronic_artwork_guidelines.pdf)
PRINT_WIDTH_IN = 170 / 25.4      # 170mm double-column journal figure width

LEFT_DOMAIN = (0.0, 0.495)   # rbcL panel
RIGHT_DOMAIN = (0.505, 1.0)  # trnL panel
SANKEY_Y_DOMAIN = (0.015, 1.0)  # leaves room at the bottom for the shared legend

lifeform_colors = {
    "Graminoid": "#228B22",
    "Woody taxa": "#8B4513",
    "Succulent": "#90EE90",
    "Forb": "#FF69B4",
    "NA": "#B0B0B0",
}


def collapse_and_rra(df):
    """Group by Genus and convert raw counts to RRA (%) per sample."""
    sample_cols = [c for c in df.columns if c.endswith("_E") or c.endswith("_I")]
    genus = df.groupby("Genus")[sample_cols].sum()
    totals = genus.sum(axis=0).replace(0, np.nan)
    rra = genus.div(totals, axis=1).fillna(0) * 100.0
    rra.reset_index(inplace=True)
    return rra


def build_panel(df, flip):
    """Build one gene panel's Sankey node/link data. Same logic as
    sankey_for_dataset() in Sankey_plot_Average_RRA_pink_indigo.py, minus
    the per-panel title/legend/export (those are shared across both panels
    in the combined figure instead)."""
    rra = collapse_and_rra(df)
    life_form_map = df.groupby("Genus")["Life form"].first().to_dict()

    e_cols = [c for c in rra.columns if c.endswith("_E")]
    i_cols = [c for c in rra.columns if c.endswith("_I")]
    rra["E_mean"] = rra[e_cols].mean(axis=1) if e_cols else 0
    rra["I_mean"] = rra[i_cols].mean(axis=1) if i_cols else 0

    rra = rra[(rra["E_mean"] > 0.5) | (rra["I_mean"] > 0.5)]

    total_mean = rra["E_mean"] + rra["I_mean"]
    rra["ratio"] = rra["I_mean"] / total_mean.where(total_mean != 0, 0)
    rra_ordered = rra.sort_values("ratio", ascending=True).copy()
    genera = rra_ordered["Genus"].tolist()
    genera_labels = ["<br>".join(textwrap.wrap(g, width=25)) for g in genera]

    rra_ordered["size"] = rra_ordered["E_mean"] + rra_ordered["I_mean"]
    total_size = rra_ordered["size"].sum()
    margin_y, min_node_height = 0.85, 0.024  # genera packed slightly tighter
    y_gen, cumul_y = [], 0.05
    if total_size > 0:
        for _, row in rra_ordered.iterrows():
            h = max((row["size"] / total_size) * margin_y, min_node_height)
            y_gen.append(cumul_y + h / 2)
            cumul_y += h
    else:
        y_gen = np.linspace(0.05, 0.95, len(genera)).tolist()
        cumul_y = 0.95
    used_span = cumul_y - 0.05
    if used_span > margin_y:
        scale_factor = margin_y / used_span
        y_gen = [(y - 0.05) * scale_factor + 0.05 for y in y_gen]

    src, tgt, val, lcolors = [], [], [], []
    e_conns, i_conns = [], []
    genus_colors = [lifeform_colors.get(life_form_map.get(g, "NA"), "#B0B0B0") for g in genera]

    # Links are added in two grouped passes (all Elephant-bound, then all
    # Impala-bound -- both still in ratio-ascending/y_gen order within their
    # own group) rather than interleaved per-row. Plotly stacks the multiple
    # links converging on a shared target (Elephant/Impala) in the order they
    # appear in the link arrays; interleaving the two colours per-row was
    # producing a stray low-value link (e.g. a thin Impala flow) crossing
    # underneath a much larger one instead of staying grouped with its
    # neighbours in source order.
    if not flip:
        nodes = ["Elephant", "Impala"] + genera_labels
        node_map = {n: i for i, n in enumerate(nodes)}
        for gi, (_, row) in enumerate(rra_ordered.iterrows()):
            if row["E_mean"] > 0:
                g_idx = node_map[genera_labels[gi]]
                src.append(node_map["Elephant"]); tgt.append(g_idx)
                val.append(row["E_mean"]); lcolors.append("rgba(255,107,107,0.45)")
                e_conns.append((y_gen[gi], row["E_mean"]))
        for gi, (_, row) in enumerate(rra_ordered.iterrows()):
            if row["I_mean"] > 0:
                g_idx = node_map[genera_labels[gi]]
                src.append(node_map["Impala"]); tgt.append(g_idx)
                val.append(row["I_mean"]); lcolors.append("rgba(78,205,196,0.45)")
                i_conns.append((y_gen[gi], row["I_mean"]))
        node_colors = ["#FF6B6B", "#4ECDC4"] + genus_colors
        x = [0.39, 0.39] + [0.995] * len(genera)
    else:
        nodes = genera_labels + ["Elephant", "Impala"]
        node_map = {n: i for i, n in enumerate(nodes)}
        for gi, (_, row) in enumerate(rra_ordered.iterrows()):
            if row["E_mean"] > 0:
                g_idx = node_map[genera_labels[gi]]
                src.append(g_idx); tgt.append(node_map["Elephant"])
                val.append(row["E_mean"]); lcolors.append("rgba(255,107,107,0.45)")
                e_conns.append((y_gen[gi], row["E_mean"]))
        for gi, (_, row) in enumerate(rra_ordered.iterrows()):
            if row["I_mean"] > 0:
                g_idx = node_map[genera_labels[gi]]
                src.append(g_idx); tgt.append(node_map["Impala"])
                val.append(row["I_mean"]); lcolors.append("rgba(78,205,196,0.45)")
                i_conns.append((y_gen[gi], row["I_mean"]))
        node_colors = genus_colors + ["#FF6B6B", "#4ECDC4"]
        x = [0.005] * len(genera) + [0.61, 0.61]

    e_y = sum(y * v for y, v in e_conns) / sum(v for _, v in e_conns) if e_conns else 0.5
    i_y = sum(y * v for y, v in i_conns) / sum(v for _, v in i_conns) if i_conns else 0.5
    mid_y = (e_y + i_y) / 2
    sep = abs(e_y - i_y) + 0.1
    e_y = max(0.05, min(0.95, mid_y - sep / 2))
    i_y = max(0.05, min(0.95, mid_y + sep / 2))

    y = ([e_y, i_y] + y_gen) if not flip else (y_gen + [e_y, i_y])
    life_forms_present = set(life_form_map.get(g, "NA") for g in genera)

    return dict(
        nodes=nodes, node_colors=node_colors, x=x, y=y,
        src=src, tgt=tgt, val=val, lcolors=lcolors,
        e_y=e_y, i_y=i_y, life_forms=life_forms_present,
    )


def panel_paper_x(domain, x_frac):
    return domain[0] + x_frac * (domain[1] - domain[0])


def sankey_y_to_paper(y_frac, domain_y=SANKEY_Y_DOMAIN):
    # Sankey node y=0 is the top of its own domain, y=1 is the bottom.
    top, bottom = domain_y[1], domain_y[0]
    return top - y_frac * (top - bottom)


def add_silhouette(fig, path, x_frac, domain, y_frac, width_frac=0.095):
    with Image.open(path) as im:
        aspect = im.width / im.height
    fig.add_layout_image(dict(
        source=Image.open(path),
        xref="paper", yref="paper",
        x=panel_paper_x(domain, x_frac),
        y=sankey_y_to_paper(y_frac),
        sizex=width_frac, sizey=width_frac / aspect,
        xanchor="center", yanchor="middle", layer="above",
    ))


def main():
    rbcl_df = pd.read_csv(HERE / "rbcl_combined_dietlist.csv")
    trnl_df = pd.read_csv(HERE / "trnl_combined_dietlist.csv")

    left = build_panel(rbcl_df, flip=False)
    right = build_panel(trnl_df, flip=True)

    # Needed up front now: the custom legend below positions itself using the
    # canvas's pixel aspect ratio so its colour-key circles render as true
    # circles, not ellipses stretched by the wide canvas.
    export_width = round(PRINT_WIDTH_IN * TARGET_DPI)
    export_height = round(export_width * 0.62)  # ~ same aspect as the old 13.33x7.5in slide

    fig = go.Figure()
    fig.add_trace(go.Sankey(
        arrangement="fixed",
        domain=dict(x=list(LEFT_DOMAIN), y=list(SANKEY_Y_DOMAIN)),
        textfont=dict(size=98, family="Arial, sans-serif"),  # genus/animal labels: even bigger per supervisor feedback
        node=dict(label=left["nodes"], pad=192, thickness=58,
                  x=left["x"], y=left["y"], color=left["node_colors"]),
        link=dict(source=left["src"], target=left["tgt"],
                  value=left["val"], color=left["lcolors"]),
    ))
    fig.add_trace(go.Sankey(
        arrangement="fixed",
        domain=dict(x=list(RIGHT_DOMAIN), y=list(SANKEY_Y_DOMAIN)),
        textfont=dict(size=98, family="Arial, sans-serif"),  # genus/animal labels: even bigger per supervisor feedback
        node=dict(label=right["nodes"], pad=182, thickness=58,
                  x=right["x"], y=right["y"], color=right["node_colors"]),
        link=dict(source=right["src"], target=right["tgt"],
                  value=right["val"], color=right["lcolors"]),
    ))

    # One shared legend for the union of life forms actually present in either
    # panel, drawn as shapes+annotations rather than Plotly's native legend --
    # Plotly's legend symbol size does not scale with marker.size for these
    # dummy marker-only traces when exported through kaleido, so this is the
    # only way to get colour keys that are actually controllable in size.
    LEGEND_FONT_SIZE = 105
    LEGEND_Y_CENTER = 0.08  # close to the bottom of the plot
    LEGEND_SWATCH_H = 0.00966  # fraction of canvas height
    LEGEND_SWATCH_W = LEGEND_SWATCH_H * (export_height / export_width)  # keeps swatches circular on the wide canvas
    GROUP_GAP = 0.0298  # at this size, small per-word text-width estimate errors become proportionally invisible
    SWATCH_LABEL_GAP = 0.00095

    # Real font metrics instead of a flat per-character average -- a fixed
    # avg-glyph-width estimate is systematically wrong by a different amount
    # for every word depending on its actual letter shapes.
    _legend_font = ImageFont.truetype("/home/kyle-linux/miniforge3/fonts/DejaVuSans.ttf", LEGEND_FONT_SIZE)

    # DejaVu's metrics don't match kaleido's actual rendered font 1:1 -- this
    # flat correction factor (empirically calibrated on the sibling FOO plot)
    # brings the estimate back in line across the board.
    DEJAVU_CORRECTION = 0.83

    def text_width_frac(s):
        bbox = _legend_font.getbbox(s)
        return (bbox[2] - bbox[0]) / export_width * DEJAVU_CORRECTION

    legend_items = [("Life form", None)] + [(lf, lifeform_colors.get(lf, "#B0B0B0"))
                                             for lf in sorted(left["life_forms"] | right["life_forms"])]
    widths = []
    for text, color in legend_items:
        w = text_width_frac(text)
        if color:
            w += LEGEND_SWATCH_W + SWATCH_LABEL_GAP
        widths.append(w)
    total_w = sum(widths) + GROUP_GAP * (len(legend_items) - 1)

    cursor_x = 0.5 - total_w / 2
    pad_x = 0.0040
    pad_y = 0.0035
    fig.add_shape(type="rect", xref="paper", yref="paper",
                  x0=cursor_x - pad_x, x1=cursor_x + total_w + pad_x,
                  y0=LEGEND_Y_CENTER - LEGEND_SWATCH_H / 2 - pad_y,
                  y1=LEGEND_Y_CENTER + LEGEND_SWATCH_H / 2 + pad_y,
                  line=dict(color="black", width=2.0), fillcolor="white", layer="below")

    for (text, color), w in zip(legend_items, widths):
        if color:
            fig.add_shape(type="circle", xref="paper", yref="paper",
                          x0=cursor_x, x1=cursor_x + LEGEND_SWATCH_W,
                          y0=LEGEND_Y_CENTER - LEGEND_SWATCH_H / 2,
                          y1=LEGEND_Y_CENTER + LEGEND_SWATCH_H / 2,
                          line=dict(width=0), fillcolor=color)
            text_x = cursor_x + LEGEND_SWATCH_W + SWATCH_LABEL_GAP
        else:
            text_x = cursor_x
        fig.add_annotation(text=text, x=text_x, y=LEGEND_Y_CENTER, xref="paper", yref="paper",
                            showarrow=False, font=dict(size=LEGEND_FONT_SIZE, color="#2a3f5f", family="Arial, sans-serif"),
                            xanchor="left", yanchor="middle")
        cursor_x += w + GROUP_GAP

    SILHOUETTE_X = 0.32  # domain-relative x shared by the silhouettes and the gene labels above them -- pulled closer to the genus columns

    rbcl_label_x = panel_paper_x(LEFT_DOMAIN, SILHOUETTE_X)
    trnl_label_x = panel_paper_x(RIGHT_DOMAIN, 1 - SILHOUETTE_X)
    rbcl_label_y = 0.965
    trnl_label_y = 0.965

    fig.add_annotation(text="<b><i>rbc</i>L</b>", x=rbcl_label_x, y=rbcl_label_y, xref="paper", yref="paper",
                        showarrow=False, font=dict(size=150), xanchor="center")
    fig.add_annotation(text="<b><i>trn</i>L</b>", x=trnl_label_x, y=trnl_label_y, xref="paper", yref="paper",
                        showarrow=False, font=dict(size=150), xanchor="center")

    add_silhouette(fig, ASSETS / "elephant.png", SILHOUETTE_X, LEFT_DOMAIN, left["e_y"])
    add_silhouette(fig, ASSETS / "impala.png", SILHOUETTE_X, LEFT_DOMAIN, left["i_y"])
    add_silhouette(fig, ASSETS / "elephant.png", 1 - SILHOUETTE_X, RIGHT_DOMAIN, right["e_y"])
    add_silhouette(fig, ASSETS / "impala.png", 1 - SILHOUETTE_X, RIGHT_DOMAIN, right["i_y"])

    fig.update_layout(
        title=dict(text="Figure A8", x=0.5, y=0.995, yanchor="top", font=dict(size=130)),
        font=dict(size=29),
        paper_bgcolor="white",
        plot_bgcolor="white",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        margin=dict(l=48, r=48, t=263, b=144),
        showlegend=False,  # custom shape-based legend above replaces Plotly's native one
    )

    # Render at 2x the target pixel size, then downsample with a high-quality
    # filter before saving -- supersampling like this makes thin Sankey links
    # and small text noticeably crisper than rendering natively at the final
    # size, without changing the final pixel dimensions or DPI at all.
    SUPERSAMPLE = 2
    out_path = HERE / "Figure_RRA_composite.png"
    fig.write_image(str(out_path), width=export_width * SUPERSAMPLE,
                     height=export_height * SUPERSAMPLE, scale=1)

    # kaleido doesn't write a DPI tag; stamp the correct one so the pixel
    # dimensions map to PRINT_WIDTH_IN at TARGET_DPI when placed in the manuscript.
    with Image.open(out_path) as im:
        im = im.resize((export_width, export_height), Image.LANCZOS)
        im.save(out_path, dpi=(TARGET_DPI, TARGET_DPI), optimize=True)

    size_kb = out_path.stat().st_size / 1024
    print(f"Saved {out_path}")
    print(f"  {export_width}x{export_height}px @ {TARGET_DPI} DPI -> {size_kb:.0f} KB")


if __name__ == "__main__":
    main()
