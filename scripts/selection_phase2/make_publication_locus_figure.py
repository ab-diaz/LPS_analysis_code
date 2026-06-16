#!/usr/bin/env python3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Patch
import numpy as np
import pandas as pd


BASE = Path("results/lps_analysis/selection_phase2")
OUT = BASE / "figures"
DATA = BASE / "pantoea_route1_strengthening"
LOC = DATA / "target_gene_bakta_locations.tsv"
NEIGH = DATA / "target_gene_bakta_neighborhoods_pm6.tsv"
SUMMARY = DATA / "oantigen_locus_architecture_summary.tsv"

GENE_COLORS = {
    "rfba": "#3f77b5",
    "rfbb": "#f28e2b",
    "rfbc": "#59a14f",
    "rfbd": "#b07aa1",
    "waal": "#e15759",
    "wzm": "#76b7b2",
    "wzt": "#ff9da7",
    "rfaz": "#9c755f",
}
GENE_LABELS = {
    "rfba": "rfbA",
    "rfbb": "rfbB",
    "rfbc": "rfbC",
    "rfbd": "rfbD",
    "waal": "waaL",
    "wzm": "wzm",
    "wzt": "wzt",
    "rfaz": "rfaZ",
}
TARGETS = set(GENE_COLORS)

REPRESENTATIVES = [
    ("ISS 2016", "Pantoea_piersonii_IIIF1SW-P2_B"),
    ("ISS 2021", "Pantoea_piersonii_F10_5S-D1_P5"),
    ("ISS 2022", "Pantoea_piersonii_F11_7S_D1_EB1"),
    ("Earth complete rfb", "Pantoea_deleyi_LMG24200"),
    ("Earth transport pair", "Pantoea_alfalfae_CQ10"),
    ("Earth rfaZ lineage", "Pantoea_latae_AS1"),
]


def choose_contig(sub, summary_row):
    if summary_row is not None and isinstance(summary_row.get("rfb_contig", ""), str) and summary_row.get("rfb_contig", ""):
        return summary_row["rfb_contig"]
    if not sub.empty:
        # Prefer a contig containing the most focal genes.
        counts = sub.groupby("contig")["gene"].nunique().sort_values(ascending=False)
        return counts.index[0]
    return None


def load_panel_a_rows():
    loc = pd.read_csv(LOC, sep="\t")
    neigh = pd.read_csv(NEIGH, sep="\t")
    summary = pd.read_csv(SUMMARY, sep="\t")
    rows = []
    for label, isolate in REPRESENTATIVES:
        sub = loc[loc["isolate"].eq(isolate)].copy()
        srow = summary[summary["isolate"].eq(isolate)]
        sdict = srow.iloc[0].to_dict() if len(srow) else None
        contig = choose_contig(sub, sdict)
        if contig is None:
            continue
        anchors = sub[sub["contig"].eq(contig)]
        rfb_anchors = anchors[anchors["gene"].isin(["rfba", "rfbb", "rfbc", "rfbd"])]
        if rfb_anchors["gene"].nunique() >= 3:
            focus = rfb_anchors
        else:
            focus = anchors
        anchor_loci = set(anchors["locus_tag"])
        # Include neighborhoods around focal genes on the selected contig.
        nsub = neigh[(neigh["isolate"].eq(isolate)) & (neigh["contig"].eq(contig))].copy()
        if nsub.empty:
            nsub = anchors.rename(columns={"product": "product"}).copy()
            nsub["anchor_locus"] = ""
        # Deduplicate by locus tag and keep local region spanning focal genes plus flanks.
        cols = ["contig", "start", "end", "strand", "locus_tag", "gene", "product"]
        nsub = nsub[cols].drop_duplicates("locus_tag").sort_values("start")
        if focus.empty:
            region = nsub
        else:
            lo = int(focus["start"].min()) - 2500
            hi = int(focus["end"].max()) + 2500
            region = nsub[(nsub["end"] >= lo) & (nsub["start"] <= hi)].copy()
            if region.empty:
                region = nsub
        target_by_locus = dict(zip(anchors["locus_tag"], anchors["gene"]))
        region["plot_gene"] = region.apply(
            lambda r: target_by_locus.get(r["locus_tag"], str(r["gene"]).lower() if pd.notna(r["gene"]) else ""),
            axis=1,
        )
        rows.append((label, isolate, contig, region))
    return rows


def draw_gene(ax, x1, x2, y, strand, color, alpha=1.0):
    if x2 < x1:
        x1, x2 = x2, x1
    width = max(x2 - x1, 120)
    if strand == "-":
        start, end = x2, x1
    else:
        start, end = x1, x2
    arrow = FancyArrowPatch(
        (start, y),
        (end, y),
        arrowstyle="-|>",
        mutation_scale=12,
        linewidth=0.8,
        facecolor=color,
        edgecolor="#222222",
        alpha=alpha,
        shrinkA=0,
        shrinkB=0,
    )
    # Add a body line so very short genes remain visible.
    ax.plot([x1, x2], [y, y], color=color, linewidth=8, alpha=alpha, solid_capstyle="butt")
    ax.add_patch(arrow)


def make_figure():
    OUT.mkdir(exist_ok=True)
    rows = load_panel_a_rows()
    summary = pd.read_csv(SUMMARY, sep="\t")

    fig = plt.figure(figsize=(11, 8.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.7, 1.0], wspace=0.22)
    ax = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    y_positions = np.arange(len(rows))[::-1]
    max_width = 0
    for y, (label, isolate, contig, region) in zip(y_positions, rows):
        if region.empty:
            continue
        start0 = int(region["start"].min()) - 300
        width = int(region["end"].max()) - start0 + 300
        max_width = max(max_width, width)
        ax.hlines(y, 0, width, color="#d9d9d9", linewidth=1)
        for _, gene in region.iterrows():
            g = str(gene["plot_gene"]).lower() if pd.notna(gene["plot_gene"]) else ""
            color = GENE_COLORS.get(g, "#c7c7c7")
            alpha = 1.0 if g in TARGETS else 0.55
            x1 = int(gene["start"]) - start0
            x2 = int(gene["end"]) - start0
            draw_gene(ax, x1, x2, y, gene["strand"], color, alpha)
        short = isolate.replace("Pantoea_piersonii_", "P. piersonii ").replace("Pantoea_", "P. ")
        ax.text(-420, y, f"{label}\n{short}", ha="right", va="center", fontsize=8)

    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_xlim(-1300, max_width + 500)
    ax.set_yticks([])
    ax.set_xlabel("Relative genomic position in local O-antigen/LPS region (bp)")
    ax.text(0.0, 1.02, "A", transform=ax.transAxes, ha="left", va="bottom", fontsize=13, fontweight="bold")
    ax.spines[["left", "right", "top"]].set_visible(False)

    legend_handles = [Patch(facecolor=GENE_COLORS[g], edgecolor="#222222", label=GENE_LABELS[g]) for g in GENE_COLORS]
    legend_handles.append(Patch(facecolor="#c7c7c7", edgecolor="#222222", label="neighboring CDS"))
    ax.legend(handles=legend_handles, loc="lower left", bbox_to_anchor=(0.0, -0.23), ncol=5, frameon=False, fontsize=8)

    metrics = [
        ("Co-localized\nrfbABCD", "rfbABCD_colocalized_15kb"),
        ("waaL", "has_waaL"),
        ("wzm", "has_wzm"),
        ("wzt", "has_wzt"),
        ("rfaZ", "has_rfaZ"),
    ]
    iss = summary[summary["group"].eq("ISS")]
    earth = summary[summary["group"].eq("Earth")]
    labels = [m[0] for m in metrics]
    iss_vals = [iss[m[1]].sum() / len(iss) for m in metrics]
    earth_vals = [earth[m[1]].sum() / len(earth) for m in metrics]
    iss_counts = [f"{int(iss[m[1]].sum())}/{len(iss)}" for m in metrics]
    earth_counts = [f"{int(earth[m[1]].sum())}/{len(earth)}" for m in metrics]

    x = np.arange(len(metrics))
    w = 0.36
    ax2.bar(x - w / 2, iss_vals, width=w, color="#b23a48", label="ISS")
    ax2.bar(x + w / 2, earth_vals, width=w, color="#3d6f9f", label="Earth")
    for xi, v, txt in zip(x - w / 2, iss_vals, iss_counts):
        ax2.text(xi, min(v + 0.035, 1.06), txt, ha="center", va="bottom", fontsize=8, rotation=90)
    for xi, v, txt in zip(x + w / 2, earth_vals, earth_counts):
        ax2.text(xi, min(v + 0.035, 1.06), txt, ha="center", va="bottom", fontsize=8, rotation=90)
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=30, ha="right")
    ax2.set_ylim(0, 1.15)
    ax2.set_ylabel("Fraction of genomes")
    ax2.text(0.0, 1.02, "B", transform=ax2.transAxes, ha="left", va="bottom", fontsize=13, fontweight="bold")
    ax2.legend(frameon=False, loc="upper right")
    ax2.spines[["top", "right"]].set_visible(False)

    fig.subplots_adjust(left=0.19, right=0.97, top=0.95, bottom=0.2)
    for ext in ["png", "pdf", "svg"]:
        out = OUT / f"figure5_oantigen_locus_architecture.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else None)
        print(out)


if __name__ == "__main__":
    make_figure()
