"""
Combined rbcL + trnL NMDS ordination (Jaccard distance, presence/absence)
-- Figure 2 -- exported directly as a print-ready PNG.

Replaces the old workflow of saving each locus's NMDS plot separately,
screenshotting them into PowerPoint side by side, and exporting to PDF --
Figure_2.pdf ended up at 960x540pt (33.87x19.05cm, ~1.9x over Wiley's 18cm
width cap) with both panels downsampled to ~200 ppi despite the source
PNGs being rendered at 600 DPI.

No recomputation needed here: NMDS_presence_absence.py already saves its
final (pre-jitter) coordinates to nmds_scores_rbcL.csv / nmds_scores_trnL.csv
before it applies jitter for its own (different-looking) visualization --
those two CSVs are exactly what produced Figure_2.pdf's point positions.
That script, its CSVs, and its own outputs are untouched by this one.

Run from this directory (needs pandas + matplotlib + Pillow, e.g. the
"metabarcoding" conda env on this machine):
    python make_figure2_composite.py
"""
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from PIL import Image

HERE = Path(__file__).resolve().parent
FOO_DIR = HERE / "FOO_NMDS_PLOT_with_PERMDISP"

TARGET_DPI = 1000  # top of Wiley's 600-1000 DPI "line art" preferred range (electronic_artwork_guidelines.pdf)
FIG_WIDTH_IN = 170 / 25.4   # under Wiley's 180mm cap
FIG_HEIGHT_IN = 95 / 25.4
BODY_FONT_SIZE = 7          # Wiley requires figure text sized 6-10pt
TICK_FONT_SIZE = 6
HEADING_FONT_SIZE = 10

GROUP_COLORS = {"Elephant": "tab:blue", "Impala": "tab:orange"}
# Shared with make_figureA7_composite.py so both figures render dots
# identically -- the originals used different alpha/size per script (0.55 vs
# 0.85, 35 vs 40), which is fine at full size but at this panel size makes
# a) the two figures look inconsistent and b) overlapping semi-transparent
# points blend into confusing muddy colours. Higher alpha + smaller size
# (matching the shrunk panel, since `s` is an absolute area, not relative)
# fixes both.
POINT_SIZE = 12
POINT_ALPHA = 0.85
POINT_EDGE_COLOR = "black"
POINT_EDGE_WIDTH = 0.3


def draw_panel(ax, scores):
    for grp, sub in scores.groupby("Group"):
        ax.scatter(sub["NMDS1"], sub["NMDS2"], label=grp, color=GROUP_COLORS[grp],
                   s=POINT_SIZE, alpha=POINT_ALPHA, edgecolors=POINT_EDGE_COLOR,
                   linewidths=POINT_EDGE_WIDTH)
    ax.set_xlabel("NMDS1", fontsize=BODY_FONT_SIZE)
    ax.set_ylabel("NMDS2", fontsize=BODY_FONT_SIZE)
    ax.tick_params(axis='both', labelsize=TICK_FONT_SIZE, width=0.6)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
    ax.set_box_aspect(1)


def main():
    rbcl_scores = pd.read_csv(FOO_DIR / "nmds_scores_rbcL.csv", index_col=0)
    trnl_scores = pd.read_csv(FOO_DIR / "nmds_scores_trnL.csv", index_col=0)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(FIG_WIDTH_IN, FIG_HEIGHT_IN),
                                      constrained_layout=True)

    draw_panel(ax_a, rbcl_scores)
    draw_panel(ax_b, trnl_scores)

    # Match panel A's x-axis tick increment to panel B's -- each panel's data
    # range differs (different marker gene), so matplotlib's auto-locator
    # otherwise picks a different step per panel, making the two harder to
    # visually compare side by side.
    b_ticks = ax_b.get_xticks()
    ax_a.xaxis.set_major_locator(MultipleLocator(b_ticks[1] - b_ticks[0]))

    for ax, letter in [(ax_a, "A"), (ax_b, "B")]:
        ax.text(-0.08, 1.05, letter, transform=ax.transAxes, fontsize=HEADING_FONT_SIZE,
                fontweight='bold')

    # One shared legend instead of a duplicate in each panel -- both panels
    # show the same two categories/colours.
    handles, labels = ax_a.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.08),
               ncol=2, fontsize=BODY_FONT_SIZE, frameon=True)

    fig.suptitle("Figure 2", fontsize=HEADING_FONT_SIZE, y=1.05)

    # Render at 2x TARGET_DPI internally, then downsample with a high-quality
    # filter to the pixel size TARGET_DPI actually implies -- sharper than
    # rendering natively at the final size, same final DPI/physical size.
    SUPERSAMPLE = 2
    out_path = HERE / "Figure_2_composite.png"
    fig.savefig(out_path, dpi=TARGET_DPI * SUPERSAMPLE, bbox_inches='tight')

    with Image.open(out_path) as im:
        target_size = (round(im.width / SUPERSAMPLE), round(im.height / SUPERSAMPLE))
        im = im.resize(target_size, Image.LANCZOS)
        im.save(out_path, dpi=(TARGET_DPI, TARGET_DPI), optimize=True)
        w_mm, h_mm = im.width / TARGET_DPI * 25.4, im.height / TARGET_DPI * 25.4
        size_kb = out_path.stat().st_size / 1024
        print(f"Saved {out_path}")
        print(f"  {im.width}x{im.height}px @ {TARGET_DPI} DPI -> {w_mm:.1f} x {h_mm:.1f} mm, {size_kb:.0f} KB")


if __name__ == "__main__":
    main()
