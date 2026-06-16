#!/usr/bin/env python3
"""Create a Mash-based phylogenomic context figure for the Pantoea dataset."""

from __future__ import annotations

import csv
import math
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


ROOT = Path(__file__).resolve().parent
IN = ROOT / "phylogeny_aware_inputs"
FIG = ROOT / "figures"
OUT = ROOT / "phylogenomic_context_mash"

GENOME_INDEX = IN / "genome_file_index.tsv"
FOCAL_MATRIX = IN / "focal_gene_matrix.tsv"
FASTA_LIST = OUT / "pantoea_92_genome_fasta_paths.txt"
MASH_PREFIX = OUT / "pantoea_92_mash"
MASH_DIST = OUT / "pantoea_92_mash_dist.tsv"

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


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_fasta_list() -> list[dict[str, str]]:
    rows = read_tsv(GENOME_INDEX)
    missing = [r["genome"] for r in rows if not r["fna"] or not Path(r["fna"]).exists()]
    if missing:
        raise SystemExit(f"Missing FASTA files for {len(missing)} genomes: {missing[:10]}")
    OUT.mkdir(parents=True, exist_ok=True)
    FASTA_LIST.write_text("\n".join(r["fna"] for r in rows) + "\n")
    return rows


def run_mash() -> None:
    msh = MASH_PREFIX.with_suffix(".msh")
    if not msh.exists():
        subprocess.run(
            ["mash", "sketch", "-s", "10000", "-k", "21", "-o", str(MASH_PREFIX), "-l", str(FASTA_LIST)],
            check=True,
        )
    if not MASH_DIST.exists():
        with MASH_DIST.open("w") as out:
            subprocess.run(["mash", "dist", str(msh), str(msh)], check=True, stdout=out)


def fasta_path_to_genome(rows: list[dict[str, str]]) -> dict[str, str]:
    mapping = {}
    for row in rows:
        path = str(Path(row["fna"]))
        mapping[path] = row["genome"]
        mapping[Path(path).name] = row["genome"]
        mapping[Path(path).stem] = row["genome"]
    return mapping


def parse_mash_distance(rows: list[dict[str, str]]) -> pd.DataFrame:
    genomes = [r["genome"] for r in rows]
    mapper = fasta_path_to_genome(rows)
    matrix = pd.DataFrame(np.nan, index=genomes, columns=genomes, dtype=float)
    for genome in genomes:
        matrix.loc[genome, genome] = 0.0

    with MASH_DIST.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            q, r, dist, pval, shared = line.rstrip("\n").split("\t")[:5]
            qg = mapper.get(q) or mapper.get(Path(q).name) or mapper.get(Path(q).stem)
            rg = mapper.get(r) or mapper.get(Path(r).name) or mapper.get(Path(r).stem)
            if not qg or not rg:
                raise ValueError(f"Could not map Mash paths: {q}, {r}")
            value = float(dist)
            matrix.loc[qg, rg] = value
            matrix.loc[rg, qg] = value

    if matrix.isna().any().any():
        missing = int(matrix.isna().sum().sum())
        raise ValueError(f"Mash matrix has {missing} missing values")
    matrix.to_csv(OUT / "pantoea_92_mash_distance_matrix.tsv", sep="\t")
    return matrix


def build_order(distance: pd.DataFrame) -> tuple[list[str], np.ndarray]:
    condensed = squareform(distance.values, checks=False)
    z = linkage(condensed, method="average")
    leaves = dendrogram(z, no_plot=True, labels=list(distance.index))["ivl"]
    return leaves, z


def write_newick(node, labels: list[str], distances: np.ndarray, parent_height: float | None = None) -> str:
    # scipy cluster node ids: leaves 0..n-1, internal n..2n-2
    n = len(labels)
    if node < n:
        length = 0.0 if parent_height is None else parent_height
        return f"{labels[node]}:{length:.8f}"
    row = distances[node - n]
    left, right, height = int(row[0]), int(row[1]), float(row[2]) / 2.0
    branch = 0.0 if parent_height is None else max(parent_height - height, 0.0)
    return f"({write_newick(left, labels, distances, height)},{write_newick(right, labels, distances, height)}):{branch:.8f}"


def save_tree_files(distance: pd.DataFrame, z: np.ndarray) -> None:
    labels = list(distance.index)
    root = len(labels) + z.shape[0] - 1
    newick = write_newick(root, labels, z, None) + ";\n"
    (OUT / "pantoea_92_mash_upgma_tree.nwk").write_text(newick)

    tree_order = dendrogram(z, no_plot=True, labels=labels)["ivl"]
    pd.DataFrame({"genome": tree_order, "tree_order": range(1, len(tree_order) + 1)}).to_csv(
        OUT / "pantoea_92_mash_tree_order.tsv", sep="\t", index=False
    )


def plot_context(distance: pd.DataFrame, z: np.ndarray) -> None:
    meta = pd.read_csv(FOCAL_MATRIX, sep="\t", dtype=str).set_index("genome")
    for gene in TARGET_GENES:
        meta[gene] = meta[gene].astype(int)
    meta["taxon_context"] = np.where(
        meta["species"].eq("Pantoea_piersonii"),
        "P. piersonii",
        np.where(meta["species"].eq("Pantoea_septica"), "P. septica", "Other Pantoea"),
    )

    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(
        nrows=1,
        ncols=6,
        width_ratios=[4.6, 0.28, 0.28, 0.28, 2.4, 1.15],
        wspace=0.04,
    )

    ax_tree = fig.add_subplot(gs[0, 0])
    dendro = dendrogram(
        z,
        orientation="left",
        labels=list(distance.index),
        ax=ax_tree,
        leaf_font_size=4.8,
        color_threshold=0,
        above_threshold_color="#374151",
        link_color_func=lambda _: "#374151",
    )
    ordered = dendro["ivl"]
    ax_tree.set_xlabel("Mash distance")
    ax_tree.set_title("A. Genome-wide relatedness", loc="left", fontsize=12, fontweight="bold")
    ax_tree.tick_params(axis="y", labelsize=4.8)
    for spine in ["top", "right"]:
        ax_tree.spines[spine].set_visible(False)

    # Dendrogram leaf y positions are 5, 15, 25, ...
    # Heatmaps need the same y-scale, otherwise they collapse into the bottom
    # tenth of the shared axis.
    heat_extent = (-0.5, 0.5, 0, len(ordered) * 10)
    gene_extent = (-0.5, len(TARGET_GENES) - 0.5, 0, len(ordered) * 10)
    ordered_meta = meta.loc[ordered]

    ax_source = fig.add_subplot(gs[0, 1], sharey=ax_tree)
    source_rgb = np.array([matplotlib.colors.to_rgb(SOURCE_COLORS[s]) for s in ordered_meta["source"]])
    ax_source.imshow(source_rgb.reshape(len(ordered), 1, 3), aspect="auto", origin="lower", extent=heat_extent)
    ax_source.set_xticks([0])
    ax_source.set_xticklabels(["Source"], rotation=90, fontsize=8)
    ax_source.tick_params(axis="y", left=False, labelleft=False)
    ax_source.set_title("", fontsize=1)

    ax_year = fig.add_subplot(gs[0, 2], sharey=ax_tree)
    year_labels = ["Earth" if s == "Earth" else str(yv) for s, yv in zip(ordered_meta["source"], ordered_meta["year"])]
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
    ax_gene.set_title("B. O-antigen/LPS genes", loc="left", fontsize=12, fontweight="bold")
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

    fig.suptitle(
        "Phylogenomic context of the ISS-associated O-antigen/LPS profile",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.02,
        0.01,
        "Genome-wide relatedness was estimated with Mash distances from 92 Pantoea assemblies. "
        "The tree summarizes genomic relatedness for context and is not interpreted as a strict population-genetic sampling design.",
        fontsize=8,
    )

    FIG.mkdir(parents=True, exist_ok=True)
    for out in [FIG / "figure2_phylogenomic_context_mash.png", OUT / "figure2_phylogenomic_context_mash.png"]:
        fig.savefig(out, dpi=300, bbox_inches="tight")
    for out in [FIG / "figure2_phylogenomic_context_mash.pdf", OUT / "figure2_phylogenomic_context_mash.pdf"]:
        fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def write_summary(rows: list[dict[str, str]], distance: pd.DataFrame) -> None:
    meta = pd.read_csv(FOCAL_MATRIX, sep="\t", dtype=str)
    iss = meta.loc[meta["source"] == "ISS", "genome"].tolist()
    earth = meta.loc[meta["source"] == "Earth", "genome"].tolist()

    def mean_pair(a: list[str], b: list[str], within: bool = False) -> float:
        vals = []
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                if within and j <= i:
                    continue
                vals.append(distance.loc[x, y])
        return float(np.mean(vals)) if vals else math.nan

    lines = [
        "# Mash phylogenomic context summary",
        "",
        f"- Genomes included: {len(rows)}",
        f"- ISS genomes: {len(iss)}",
        f"- Earth genomes: {len(earth)}",
        "- Mash parameters: sketch size 10000, k-mer size 21",
        f"- Mean within-ISS Mash distance: {mean_pair(iss, iss, within=True):.6f}",
        f"- Mean within-Earth Mash distance: {mean_pair(earth, earth, within=True):.6f}",
        f"- Mean ISS-vs-Earth Mash distance: {mean_pair(iss, earth):.6f}",
        "",
        "Outputs:",
        "- `pantoea_92_mash_distance_matrix.tsv`",
        "- `pantoea_92_mash_upgma_tree.nwk`",
        "- `pantoea_92_mash_tree_order.tsv`",
        "- `figure2_phylogenomic_context_mash.png/pdf`",
    ]
    (OUT / "mash_phylogenomic_context_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    rows = write_fasta_list()
    run_mash()
    distance = parse_mash_distance(rows)
    ordered, z = build_order(distance)
    save_tree_files(distance, z)
    plot_context(distance, z)
    write_summary(rows, distance)
    print(OUT / "mash_phylogenomic_context_summary.md")
    print(FIG / "figure2_phylogenomic_context_mash.png")
    print(FIG / "figure2_phylogenomic_context_mash.pdf")


if __name__ == "__main__":
    main()
