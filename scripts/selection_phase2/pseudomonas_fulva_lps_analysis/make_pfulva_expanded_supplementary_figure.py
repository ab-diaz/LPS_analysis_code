#!/usr/bin/env python3
"""Make a supplementary figure for the expanded Pseudomonas fulva side analysis."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap


ROOT = Path(__file__).resolve().parent
FIG_DIR = ROOT.parent / "figures"
META = ROOT / "genome_metadata_pfulva_expanded_with_mt2_iss.tsv"
MATRIX = ROOT / "pfulva_expanded_o_antigen_gene_presence.tsv"
OUT_PNG = FIG_DIR / "supplementary_pfulva_lps_oantigen_exploratory.png"
OUT_PDF = FIG_DIR / "supplementary_pfulva_lps_oantigen_exploratory.pdf"

KEY_GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "wbpL", "wbpM", "wbpV", "wbpY", "wbpZ", "wzm", "wzt", "wzx"]
GROUP_ORDER = ["ISS", "Earth_spacecraft_associated", "Earth_type_strain", "Earth_reference"]
GROUP_LABELS = {
    "ISS": "ISS MT-2",
    "Earth_spacecraft_associated": "Mars Odyssey",
    "Earth_type_strain": "Type strain",
    "Earth_reference": "Earth refs",
}
GROUP_COLORS = {
    "ISS": "#2f6fbb",
    "Earth_spacecraft_associated": "#9b6b2f",
    "Earth_type_strain": "#555555",
    "Earth_reference": "#3b8c59",
}


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


def ordered_isolates(isolates, meta, matrix):
    ordered = []
    for group in GROUP_ORDER:
        group_isos = [iso for iso in isolates if meta[iso]["group"] == group]
        group_isos.sort(key=lambda iso: tuple(matrix[g][iso] for g in KEY_GENES if g in matrix), reverse=True)
        ordered.extend(group_isos)
    return ordered


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    meta = read_meta()
    isolates, matrix = read_matrix()
    genes = [g for g in KEY_GENES if g in matrix]
    isolates = ordered_isolates(isolates, meta, matrix)

    heat = np.array([[matrix[g][iso] for iso in isolates] for g in genes], dtype=int)
    group_counts = {group: sum(1 for iso in isolates if meta[iso]["group"] == group) for group in GROUP_ORDER}
    group_freq = {
        group: [100 * sum(matrix[g][iso] for iso in isolates if meta[iso]["group"] == group) / group_counts[group] for g in genes]
        for group in GROUP_ORDER
        if group_counts[group]
    }

    fig = plt.figure(figsize=(14, 7.2), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.85], height_ratios=[0.16, 1.0], wspace=0.18, hspace=0.05)
    ax_bar = fig.add_subplot(gs[1, 0])
    ax_group = fig.add_subplot(gs[0, 1])
    ax_heat = fig.add_subplot(gs[1, 1])

    x = np.arange(len(genes))
    width = 0.18
    offsets = np.linspace(-1.5 * width, 1.5 * width, len(group_freq))
    for offset, (group, vals) in zip(offsets, group_freq.items()):
        ax_bar.bar(x + offset, vals, width=width, label=GROUP_LABELS[group], color=GROUP_COLORS[group], edgecolor="white", linewidth=0.5)
    ax_bar.set_ylim(0, 105)
    ax_bar.set_ylabel("Genomes with gene (%)")
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(genes, rotation=45, ha="right", fontsize=9)
    ax_bar.set_title("A  Gene presence by group", loc="left", fontsize=12, fontweight="bold")
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.grid(axis="y", color="#dddddd", linewidth=0.7)
    ax_bar.legend(frameon=False, fontsize=8, ncol=1, loc="upper right")

    cmap = ListedColormap(["#f2f2f2", "#1f4e79"])
    ax_heat.imshow(heat, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    ax_heat.set_yticks(np.arange(len(genes)))
    ax_heat.set_yticklabels(genes, fontsize=9)
    ax_heat.set_xticks([])
    ax_heat.set_title("B  O-antigen/LPS gene profile across genomes", loc="left", fontsize=12, fontweight="bold")
    ax_heat.set_xlabel("Genomes grouped by source")
    ax_heat.set_xticks(np.arange(-0.5, len(isolates), 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, len(genes), 1), minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=0.4)
    ax_heat.tick_params(which="minor", bottom=False, left=False)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    group_row = np.array([[GROUP_ORDER.index(meta[iso]["group"]) for iso in isolates]])
    group_cmap = ListedColormap([GROUP_COLORS[g] for g in GROUP_ORDER])
    ax_group.imshow(group_row, aspect="auto", interpolation="nearest", cmap=group_cmap, vmin=0, vmax=len(GROUP_ORDER)-1)
    ax_group.set_axis_off()

    start = 0
    for group in GROUP_ORDER:
        count = group_counts[group]
        if count == 0:
            continue
        end = start + count
        ax_heat.axvline(start - 0.5, color="#444444", linewidth=0.8)
        ax_group.text((start + end - 1) / 2, -0.9, f"{GROUP_LABELS[group]} (n={count})", ha="center", va="bottom", fontsize=8)
        start = end
    ax_heat.axvline(len(isolates) - 0.5, color="#444444", linewidth=0.8)

    fig.suptitle("Exploratory Pseudomonas fulva LPS/O-antigen gene-profile comparison", fontsize=14, fontweight="bold", y=0.98)
    fig.text(
        0.5,
        0.015,
        "Blue cells indicate gene detection by Bakta/eggNOG-based screening. ISS MT-2 isolates are from one collection date and are shown as exploratory cross-genus context.",
        ha="center",
        fontsize=8,
    )
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    print(OUT_PNG)
    print(OUT_PDF)


if __name__ == "__main__":
    main()
