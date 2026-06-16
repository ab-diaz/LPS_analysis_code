#!/usr/bin/env python3
"""Plot PhyloPhlAn/IQ-TREE marker-gene ML tree with O-antigen/LPS heatmap."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo


ROOT = Path(__file__).resolve().parent
FOCAL_MATRIX = ROOT / "phylogeny_aware_inputs" / "focal_gene_matrix.tsv"
TREE = (
    ROOT
    / "phylogenomic_context_phylophlan"
    / "output"
    / "pantoea_92_phylophlan_marker_tree"
    / "input_proteomes.tre.treefile"
)
IQTREE = TREE.with_suffix(".iqtree")
ALIGNMENT = TREE.parent / "input_proteomes_concatenated.aln"
FIG = ROOT / "figures"
OUT = ROOT / "phylogenomic_context_phylophlan"

TARGET_GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL", "wzm", "wzt", "rfaZ"]

SOURCE_COLORS = {"ISS": "#1f77b4", "Earth": "#6b7280"}
YEAR_COLORS = {
    "2016": "#7b3294",
    "2018": "#008837",
    "2021": "#e66101",
    "2022": "#5e3c99",
    "Earth": "#d1d5db",
}
TAXON_COLORS = {
    "P. piersonii": "#0f766e",
    "P. septica": "#a855f7",
    "Other Pantoea": "#cbd5e1",
}
GENE_PRESENT = "#2563eb"
GENE_ABSENT = "#f3f4f6"


def display_tip_name(name: str) -> str:
    """Compact labels keep all tips visible without making the panel unusable."""
    if name.startswith("Pantoea_piersonii_"):
        return name.replace("Pantoea_piersonii_", "Pp_")
    if name.startswith("Pantoea_septica_"):
        return name.replace("Pantoea_septica_", "Ps_")
    if name.startswith("Pantoea_"):
        return name.replace("Pantoea_", "P_")
    return name


def infer_iqtree_stats() -> dict[str, str]:
    stats = {"alignment_sequences": "", "alignment_sites": "", "model": "", "log_likelihood": ""}
    if not IQTREE.exists():
        return stats
    for line in IQTREE.read_text(errors="ignore").splitlines():
        if line.startswith("Input data:"):
            # Input data: 92 sequences with 74089 amino-acid sites
            parts = line.split()
            stats["alignment_sequences"] = parts[2]
            stats["alignment_sites"] = parts[5]
        elif line.startswith("Model of substitution:"):
            stats["model"] = line.split(":", 1)[1].strip()
        elif line.startswith("Log-likelihood of the tree:"):
            stats["log_likelihood"] = line.split(":", 1)[1].strip().split()[0]
    return stats


def ladderize_for_display(tree):
    """Orient clades so the P. piersonii/ISS-rich context is visually compact."""
    meta = pd.read_csv(FOCAL_MATRIX, sep="\t", dtype=str).set_index("genome")

    def score(clade):
        tips = [t.name for t in clade.get_terminals()]
        iss = sum(1 for t in tips if meta.loc[t, "source"] == "ISS") if tips else 0
        piersonii = sum(1 for t in tips if meta.loc[t, "species"] == "Pantoea_piersonii") if tips else 0
        return (iss, piersonii, len(tips))

    for clade in tree.find_clades(order="postorder"):
        if clade.clades:
            clade.clades.sort(key=score)
    return tree


def tree_depths(tree):
    depths = tree.depths()
    max_depth = max(depths[t] for t in tree.get_terminals())
    return depths, max_depth


def assign_y_positions(tree):
    terminals = tree.get_terminals()
    y_pos = {terminal: i * 10 + 5 for i, terminal in enumerate(terminals)}
    for clade in tree.get_nonterminals(order="postorder"):
        y_pos[clade] = float(np.mean([y_pos[c] for c in clade.clades]))
    return y_pos


def draw_tree(ax, tree):
    depths, max_depth = tree_depths(tree)
    y_pos = assign_y_positions(tree)

    for clade in tree.find_clades(order="preorder"):
        x = depths[clade]
        y = y_pos[clade]
        if clade.clades:
            child_ys = [y_pos[c] for c in clade.clades]
            ax.plot([x, x], [min(child_ys), max(child_ys)], color="#374151", lw=0.8)
            for child in clade.clades:
                ax.plot([x, depths[child]], [y_pos[child], y_pos[child]], color="#374151", lw=0.8)

    label_offset = max_depth * 0.006
    for terminal in tree.get_terminals():
        ax.text(
            depths[terminal] + label_offset,
            y_pos[terminal],
            display_tip_name(terminal.name),
            va="center",
            ha="left",
            fontsize=3.2,
            color="#111827",
            clip_on=False,
        )

    ax.set_ylim(0, len(tree.get_terminals()) * 10)
    ax.set_xlim(-max_depth * 0.01, max_depth * 1.18)
    ax.set_xlabel("Substitutions/site")
    ax.set_yticks([])
    ax.tick_params(axis="y", which="both", left=False, right=False, labelleft=False)
    ax.text(0.0, 1.01, "A", transform=ax.transAxes, ha="left", va="bottom", fontsize=12, fontweight="bold")
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    return [t.name for t in tree.get_terminals()]


def plot() -> None:
    tree = Phylo.read(TREE, "newick")
    # IQ-TREE writes an unrooted tree with an arbitrary drawn root. Rooting the
    # display on a distant non-P. piersonii taxon avoids splitting the ISS-rich
    # P. piersonii context across the top and bottom of the figure.
    tree.root_with_outgroup("Pantoea_alhagi_LTYR-11Z")
    tree = ladderize_for_display(tree)
    tips = [t.name for t in tree.get_terminals()]

    meta = pd.read_csv(FOCAL_MATRIX, sep="\t", dtype=str).set_index("genome")
    missing = sorted(set(tips) - set(meta.index))
    if missing:
        raise SystemExit(f"Tree tips missing from metadata: {missing[:10]}")
    for gene in TARGET_GENES:
        meta[gene] = meta[gene].astype(int)
    meta["taxon_context"] = np.where(
        meta["species"].eq("Pantoea_piersonii"),
        "P. piersonii",
        np.where(meta["species"].eq("Pantoea_septica"), "P. septica", "Other Pantoea"),
    )
    ordered_meta = meta.loc[tips]

    fig = plt.figure(figsize=(18, 13))
    gs = fig.add_gridspec(
        nrows=1,
        ncols=6,
        width_ratios=[6.1, 0.28, 0.28, 0.28, 2.4, 1.15],
        wspace=0.04,
    )

    ax_tree = fig.add_subplot(gs[0, 0])
    ordered = draw_tree(ax_tree, tree)
    heat_extent = (-0.5, 0.5, 0, len(ordered) * 10)
    gene_extent = (-0.5, len(TARGET_GENES) - 0.5, 0, len(ordered) * 10)

    ax_source = fig.add_subplot(gs[0, 1], sharey=ax_tree)
    source_rgb = np.array([matplotlib.colors.to_rgb(SOURCE_COLORS[s]) for s in ordered_meta["source"]])
    ax_source.imshow(source_rgb.reshape(len(ordered), 1, 3), aspect="auto", origin="lower", extent=heat_extent)
    ax_source.set_xticks([0])
    ax_source.set_xticklabels(["Source"], rotation=90, fontsize=8)
    ax_source.tick_params(axis="y", left=False, labelleft=False)

    ax_year = fig.add_subplot(gs[0, 2], sharey=ax_tree)
    year_labels = ["Earth" if s == "Earth" else str(y) for s, y in zip(ordered_meta["source"], ordered_meta["year"])]
    year_rgb = np.array([matplotlib.colors.to_rgb(YEAR_COLORS.get(v, "#ffffff")) for v in year_labels])
    ax_year.imshow(year_rgb.reshape(len(ordered), 1, 3), aspect="auto", origin="lower", extent=heat_extent)
    ax_year.set_xticks([0])
    ax_year.set_xticklabels(["Year"], rotation=90, fontsize=8)
    ax_year.tick_params(axis="y", left=False, labelleft=False)

    ax_taxon = fig.add_subplot(gs[0, 3], sharey=ax_tree)
    taxon_rgb = np.array([matplotlib.colors.to_rgb(TAXON_COLORS[t]) for t in ordered_meta["taxon_context"]])
    ax_taxon.imshow(taxon_rgb.reshape(len(ordered), 1, 3), aspect="auto", origin="lower", extent=heat_extent)
    ax_taxon.set_xticks([0])
    ax_taxon.set_xticklabels(["Taxon"], rotation=90, fontsize=8)
    ax_taxon.tick_params(axis="y", left=False, labelleft=False)

    ax_gene = fig.add_subplot(gs[0, 4], sharey=ax_tree)
    gene_matrix = ordered_meta[TARGET_GENES].astype(int).values
    gene_cmap = matplotlib.colors.ListedColormap([GENE_ABSENT, GENE_PRESENT])
    ax_gene.imshow(gene_matrix, aspect="auto", origin="lower", cmap=gene_cmap, vmin=0, vmax=1, extent=gene_extent)
    ax_gene.set_xticks(range(len(TARGET_GENES)))
    ax_gene.set_xticklabels(TARGET_GENES, rotation=90, fontsize=8)
    ax_gene.tick_params(axis="y", left=False, labelleft=False)
    ax_gene.text(0.0, 1.01, "B", transform=ax_gene.transAxes, ha="left", va="bottom", fontsize=12, fontweight="bold")
    ax_gene.set_yticks([])
    ax_gene.set_xticks(np.arange(-0.5, len(TARGET_GENES), 1), minor=True)
    ax_gene.set_yticks(np.arange(0, len(ordered) * 10 + 1, 10), minor=True)
    ax_gene.grid(which="minor", color="white", linewidth=0.25)
    ax_gene.tick_params(which="minor", bottom=False, left=False)

    ax_leg = fig.add_subplot(gs[0, 5])
    ax_leg.axis("off")
    legend_items = [
        ("Source", None),
        ("ISS", SOURCE_COLORS["ISS"]),
        ("Earth", SOURCE_COLORS["Earth"]),
        ("", None),
        ("ISS year", None),
        ("2016", YEAR_COLORS["2016"]),
        ("2018", YEAR_COLORS["2018"]),
        ("2021", YEAR_COLORS["2021"]),
        ("2022", YEAR_COLORS["2022"]),
        ("", None),
        ("Taxon", None),
        ("P. piersonii", TAXON_COLORS["P. piersonii"]),
        ("P. septica", TAXON_COLORS["P. septica"]),
        ("Other Pantoea", TAXON_COLORS["Other Pantoea"]),
        ("", None),
        ("Gene", None),
        ("Present", GENE_PRESENT),
        ("Absent", GENE_ABSENT),
    ]
    y0 = 0.98
    for label, color in legend_items:
        if label == "":
            y0 -= 0.035
            continue
        if color is None:
            ax_leg.text(0.0, y0, label, fontsize=9, fontweight="bold", va="top")
        else:
            ax_leg.add_patch(plt.Rectangle((0.0, y0 - 0.022), 0.12, 0.022, color=color, transform=ax_leg.transAxes))
            ax_leg.text(0.16, y0 - 0.011, label, fontsize=8, va="center")
        y0 -= 0.038

    stats = infer_iqtree_stats()
    FIG.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    for out in [
        FIG / "figure2_phylogenomic_context_phylophlan_ml.png",
        OUT / "figure2_phylogenomic_context_phylophlan_ml.png",
    ]:
        fig.savefig(out, dpi=300, bbox_inches="tight")
    for out in [
        FIG / "figure2_phylogenomic_context_phylophlan_ml.pdf",
        OUT / "figure2_phylogenomic_context_phylophlan_ml.pdf",
    ]:
        fig.savefig(out, bbox_inches="tight")
    plt.close(fig)

    table = ordered_meta.reset_index().rename(columns={"index": "genome"})
    table.insert(0, "tree_order", range(1, len(table) + 1))
    keep = ["tree_order", "genome", "source", "year", "species"] + TARGET_GENES + [
        "strict_lps_profile_cluster",
        "target_profile_cluster",
    ]
    table[keep].to_csv(OUT / "phylophlan_ml_tree_order_with_focal_genes.tsv", sep="\t", index=False)

    summary = [
        "# PhyloPhlAn marker-gene ML context summary",
        "",
        f"- Tree file: `{TREE}`",
        f"- Genomes in tree: {len(tips)}",
        f"- Alignment: `{ALIGNMENT}`",
        f"- Alignment sequences: {stats.get('alignment_sequences', '')}",
        f"- Alignment sites: {stats.get('alignment_sites', '')}",
        f"- IQ-TREE model: {stats.get('model', '')}",
        f"- Log-likelihood: {stats.get('log_likelihood', '')}",
        "- Display root: `Pantoea_alhagi_LTYR-11Z`",
        "- Figure output: `figures/figure2_phylogenomic_context_phylophlan_ml.png/pdf`",
        "- Supplementary tree-order table: `phylogenomic_context_phylophlan/phylophlan_ml_tree_order_with_focal_genes.tsv`",
    ]
    (OUT / "phylophlan_ml_context_summary.md").write_text("\n".join(summary) + "\n")


if __name__ == "__main__":
    plot()
