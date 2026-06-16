#!/usr/bin/env python3
"""Create main-text Figure 8 for the exploratory Pseudomonas fulva analysis.

The file name keeps the historical "figure7" stem, but the manuscript uses
this artwork as Figure 8 because Figure 7 is already the RELAX/GARD panel.
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parent
FIG_DIR = ROOT.parent / "figures"
META = ROOT / "genome_metadata_pfulva_expanded_with_mt2_iss.tsv"
MATRIX = ROOT / "pfulva_expanded_o_antigen_gene_presence.tsv"
OUT_PNG = FIG_DIR / "figure7_pseudomonas_fulva_exploratory_lps_oantigen.png"
OUT_PDF = FIG_DIR / "figure7_pseudomonas_fulva_exploratory_lps_oantigen.pdf"

GROUP_ORDER = ["ISS", "Earth_spacecraft_associated", "Earth_type_strain", "Earth_reference"]
GROUP_LABELS = {
    "ISS": "ISS MT-2",
    "Earth_spacecraft_associated": "Mars Odyssey",
    "Earth_type_strain": "Type strain",
    "Earth_reference": "Earth references",
}
GROUP_COLORS = {
    "ISS": "#2f6fbb",
    "Earth_spacecraft_associated": "#9b6b2f",
    "Earth_type_strain": "#5c5c5c",
    "Earth_reference": "#3b8c59",
}

FOCUSED_GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "wzm", "wzt"]
HEATMAP_GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "wbpL", "wbpM", "wbpV", "wbpY", "wbpZ", "wzm", "wzt", "wzx"]


def read_meta():
    with META.open(newline="", encoding="utf-8") as handle:
        return {r["isolate_id"]: r for r in csv.DictReader(handle, delimiter="\t")}


def read_matrix():
    with MATRIX.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        isolates = header[1:]
        matrix = {}
        for row in reader:
            matrix[row[0]] = {iso: int(v) for iso, v in zip(isolates, row[1:])}
    return isolates, matrix


def ordered_isolates(isolates, meta, matrix, genes):
    ordered = []
    for group in GROUP_ORDER:
        group_isolates = [iso for iso in isolates if meta[iso]["group"] == group]
        group_isolates.sort(
            key=lambda iso: (tuple(matrix[g][iso] for g in genes if g in matrix), iso),
            reverse=True,
        )
        ordered.extend(group_isolates)
    return ordered


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    meta = read_meta()
    isolates, matrix = read_matrix()
    focused_genes = [g for g in FOCUSED_GENES if g in matrix]
    heatmap_genes = [g for g in HEATMAP_GENES if g in matrix]
    isolates = ordered_isolates(isolates, meta, matrix, heatmap_genes)
    group_counts = Counter(meta[iso]["group"] for iso in isolates)

    fig = plt.figure(figsize=(11.6, 7.4))
    gs = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.02, 1.72, 0.18],
        left=0.08,
        right=0.985,
        top=0.88,
        bottom=0.14,
        hspace=0.34,
    )
    ax_freq = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[1, 0])
    ax_strip = fig.add_subplot(gs[2, 0])

    # Panel A: focused gene frequencies.
    x = np.arange(len(focused_genes))
    width = 0.18
    freq_groups = ["ISS", "Earth_reference", "Earth_type_strain", "Earth_spacecraft_associated"]
    offsets = np.linspace(-1.5 * width, 1.5 * width, len(freq_groups))
    for offset, group in zip(offsets, freq_groups):
        group_isolates = [iso for iso in isolates if meta[iso]["group"] == group]
        vals = [
            100 * sum(matrix[gene][iso] for iso in group_isolates) / len(group_isolates)
            if group_isolates
            else 0
            for gene in focused_genes
        ]
        ax_freq.bar(
            x + offset,
            vals,
            width=width,
            color=GROUP_COLORS[group],
            edgecolor="white",
            linewidth=0.5,
            label=f"{GROUP_LABELS[group]} (n={len(group_isolates)})",
        )
    ax_freq.set_ylim(0, 108)
    ax_freq.set_ylabel("Genomes with gene (%)")
    ax_freq.set_xticks(x)
    ax_freq.set_xticklabels(focused_genes, fontsize=9)
    ax_freq.text(0.0, 1.03, "A", transform=ax_freq.transAxes, ha="left", va="bottom", fontsize=12, fontweight="bold")
    ax_freq.spines[["top", "right"]].set_visible(False)
    ax_freq.grid(axis="y", color="#dddddd", linewidth=0.7)

    # Panel B: compact heatmap.
    heat = np.array([[matrix[gene][iso] for iso in isolates] for gene in heatmap_genes], dtype=int)
    cmap = ListedColormap(["#f1f1f1", "#1f4e79"])
    ax_heat.imshow(heat, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    ax_heat.set_yticks(np.arange(len(heatmap_genes)))
    ax_heat.set_yticklabels(heatmap_genes, fontsize=9)
    ax_heat.set_xticks([])
    ax_heat.text(0.0, 1.03, "B", transform=ax_heat.transAxes, ha="left", va="bottom", fontsize=12, fontweight="bold")
    ax_heat.set_xticks(np.arange(-0.5, len(isolates), 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, len(heatmap_genes), 1), minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=0.35)
    ax_heat.tick_params(which="minor", bottom=False, left=False)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    # Source-composition strip directly below the heatmap.
    group_row = np.array([[GROUP_ORDER.index(meta[iso]["group"]) for iso in isolates]])
    group_cmap = ListedColormap([GROUP_COLORS[g] for g in GROUP_ORDER])
    ax_strip.imshow(group_row, aspect="auto", interpolation="nearest", cmap=group_cmap, vmin=0, vmax=len(GROUP_ORDER) - 1)
    ax_strip.set_yticks([])
    ax_strip.set_xticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)

    start = 0
    for group in GROUP_ORDER:
        count = group_counts[group]
        if count == 0:
            continue
        end = start + count
        ax_heat.axvline(start - 0.5, color="#444444", linewidth=0.75)
        ax_strip.axvline(start - 0.5, color="#444444", linewidth=0.75)
        start = end
    ax_heat.axvline(len(isolates) - 0.5, color="#444444", linewidth=0.75)
    ax_strip.axvline(len(isolates) - 0.5, color="#444444", linewidth=0.75)

    presence_handles = [
        Patch(facecolor="#1f4e79", edgecolor="none", label="Gene detected"),
        Patch(facecolor="#f1f1f1", edgecolor="#999999", label="Not detected"),
    ]
    source_handles = [Patch(facecolor=GROUP_COLORS[g], edgecolor="none", label=f"{GROUP_LABELS[g]} (n={group_counts[g]})") for g in GROUP_ORDER]
    fig.legend(
        handles=presence_handles + source_handles,
        frameon=False,
        fontsize=8,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.045),
        ncol=6,
        columnspacing=1.0,
        handlelength=1.3,
    )

    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    print(OUT_PNG)
    print(OUT_PDF)


if __name__ == "__main__":
    main()
