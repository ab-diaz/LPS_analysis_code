#!/usr/bin/env python3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch


LPS = Path("results/lps_analysis")
OUTDIR = LPS / "selection_phase2" / "figures"

COLORS = {
    "Conserved": "#287a8e",
    "ISS-enriched": "#2f6fb0",
    "Earth-enriched": "#d26a45",
    "Variable": "#d6d6d6",
}
EDGE = "#2f2f2f"


MODULES = [
    {
        "title": "Lipid A biosynthesis",
        "subtitle": "Basal outer-membrane anchor",
        "x": 0.7,
        "y": 5.0,
        "genes": [["lpxA", "lpxC", "lpxD"], ["lpxH", "lpxB", "lpxK"], ["lpxL", "lpxM"]],
    },
    {
        "title": "KDO / lipid A-core linkage",
        "subtitle": "Conserved LPS backbone",
        "x": 4.2,
        "y": 5.0,
        "genes": [["kdsD", "kdsC", "kdsA"], ["kdsB", "kdtA"]],
    },
    {
        "title": "Core oligosaccharide assembly",
        "subtitle": "Inner/core OS construction",
        "x": 7.45,
        "y": 5.0,
        "genes": [["gmhA", "gmhB"], ["rfaC", "rfaF", "rfaG"], ["rfaP", "rfaQ"]],
    },
    {
        "title": "O-antigen precursor biosynthesis",
        "subtitle": "dTDP-L-rhamnose / variable surface sugars",
        "x": 1.2,
        "y": 2.05,
        "genes": [["rfbA", "rfbB"], ["rfbC", "rfbD"]],
    },
    {
        "title": "O-antigen ligation and export",
        "subtitle": "Surface assembly interface",
        "x": 5.9,
        "y": 2.05,
        "genes": [["waaL"], ["wzm", "wzt"]],
    },
]


def load_gene_status():
    prev = pd.read_csv(LPS / "lps_strict_prevalence_diff.tsv", sep="\t")
    fisher = pd.read_csv(LPS / "lps_strict_gene_fisher.tsv", sep="\t")
    fisher = fisher.set_index("gene_id")

    status = {}
    for _, row in prev.iterrows():
        gene = row["gene_id"]
        iss = float(row["iss_frac"])
        earth = float(row["earth_frac"])
        q = np.nan
        if gene in fisher.index:
            q = float(fisher.loc[gene, "q_value"])
        if q <= 0.05 and iss > earth:
            status[gene] = "ISS-enriched"
        elif q <= 0.05 and earth > iss:
            status[gene] = "Earth-enriched"
        elif iss >= 0.85 and earth >= 0.85:
            status[gene] = "Conserved"
        elif iss > 0 or earth > 0:
            status[gene] = "Variable"

    # Prefer the familiar waa/rfa aliases used in the manuscript.
    aliases = {
        "waaC": "rfaC",
        "waaF": "rfaF",
        "waaG": "rfaG",
        "waaP": "rfaP",
        "waaQ": "rfaQ",
        "kdtA": "kdtA",
    }
    for alias, canonical in aliases.items():
        if canonical in status:
            status[alias] = status[canonical]

    return status


def gene_box(ax, x, y, label, status):
    w, h = 0.62, 0.34
    color = COLORS.get(status, COLORS["Variable"])
    box = FancyBboxPatch(
        (x - w / 2, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.035",
        linewidth=0.8,
        edgecolor=EDGE,
        facecolor=color,
    )
    ax.add_patch(box)
    text_color = "white" if status in {"Conserved", "ISS-enriched", "Earth-enriched"} else "#222222"
    ax.text(x, y, label, ha="center", va="center", fontsize=8.2, color=text_color, weight="bold")


def draw_module(ax, module, status):
    x0, y0 = module["x"], module["y"]
    rows = module["genes"]
    width = max(len(r) for r in rows) * 0.82 + 0.55
    height = len(rows) * 0.52 + 0.94
    frame = FancyBboxPatch(
        (x0 - 0.18, y0 - height + 0.22),
        width,
        height,
        boxstyle="round,pad=0.08,rounding_size=0.09",
        linewidth=1.0,
        edgecolor="#808080",
        facecolor="#fbfbfb",
    )
    ax.add_patch(frame)
    ax.text(x0 + width / 2 - 0.18, y0 + 0.02, module["title"], ha="center", va="bottom", fontsize=10, weight="bold")
    ax.text(x0 + width / 2 - 0.18, y0 - 0.18, module["subtitle"], ha="center", va="top", fontsize=7.6, color="#4b4b4b")

    positions = {}
    start_y = y0 - 0.68
    for r_idx, row in enumerate(rows):
        row_y = start_y - r_idx * 0.52
        row_width = (len(row) - 1) * 0.78
        start_x = x0 + width / 2 - 0.18 - row_width / 2
        for c_idx, gene in enumerate(row):
            gx = start_x + c_idx * 0.78
            gene_box(ax, gx, row_y, gene, status.get(gene, "Variable"))
            positions[gene] = (gx, row_y)
            if c_idx > 0:
                prev_x = start_x + (c_idx - 1) * 0.78
                arrow(ax, prev_x + 0.34, row_y, gx - 0.34, row_y, lw=0.7)
    return positions


def arrow(ax, x1, y1, x2, y2, lw=1.1, color="#555555", rad=0):
    arr = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="-|>",
        mutation_scale=8,
        linewidth=lw,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
    )
    ax.add_patch(arr)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    status = load_gene_status()

    fig, ax = plt.subplots(figsize=(11.2, 6.2))
    ax.set_xlim(0, 11)
    ax.set_ylim(-0.25, 6.8)
    ax.axis("off")

    ax.text(
        5.5,
        6.50,
        "Simplified LPS biosynthesis architecture in ISS- and Earth-derived Pantoea genomes",
        ha="center",
        va="center",
        fontsize=13,
        weight="bold",
    )
    ax.text(
        5.5,
        6.16,
        "Basal lipid A and core LPS modules are broadly conserved; environment-associated signals concentrate in O-antigen synthesis, ligation, and export.",
        ha="center",
        va="center",
        fontsize=9,
        color="#444444",
    )

    positions = {}
    for module in MODULES:
        positions.update(draw_module(ax, module, status))

    # Flow between conserved backbone modules.
    arrow(ax, 3.25, 4.05, 4.18, 4.05, lw=1.2)
    arrow(ax, 6.60, 4.05, 7.35, 4.05, lw=1.2)

    # Downstream connection from core to O-antigen interface.
    arrow(ax, 8.9, 3.10, 7.35, 2.88, lw=1.0, rad=-0.10)
    arrow(ax, 3.35, 2.75, 5.72, 2.75, lw=1.0)

    # Module-level interpretation callouts.
    callout = FancyBboxPatch(
        (0.62, 0.05),
        4.55,
        0.36,
        boxstyle="round,pad=0.05,rounding_size=0.08",
        linewidth=0.8,
        edgecolor="#287a8e",
        facecolor="#eef7f8",
    )
    ax.add_patch(callout)
    ax.text(2.9, 0.23, "Lipid A and core LPS biosynthesis: conserved in both groups", ha="center", va="center", fontsize=8.1, color="#1f5d6d")

    callout2 = FancyBboxPatch(
        (5.62, 0.05),
        4.85,
        0.36,
        boxstyle="round,pad=0.05,rounding_size=0.08",
        linewidth=0.8,
        edgecolor="#7a7a7a",
        facecolor="#f7f7f7",
    )
    ax.add_patch(callout2)
    ax.text(8.05, 0.23, "Major Earth/ISS differences occur downstream in O-antigen genes", ha="center", va="center", fontsize=8.1, color="#333333")

    handles = [
        Patch(facecolor=COLORS["Conserved"], edgecolor=EDGE, label="Conserved in ISS and Earth"),
        Patch(facecolor=COLORS["ISS-enriched"], edgecolor=EDGE, label="ISS-enriched"),
        Patch(facecolor=COLORS["Earth-enriched"], edgecolor=EDGE, label="Earth-enriched"),
        Patch(facecolor=COLORS["Variable"], edgecolor=EDGE, label="Variable / not emphasized"),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.13), ncol=4, fontsize=8)

    png = OUTDIR / "figure6B_simplified_lps_schematic.png"
    pdf = OUTDIR / "figure6B_simplified_lps_schematic.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
