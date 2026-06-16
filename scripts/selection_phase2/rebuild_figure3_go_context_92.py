#!/usr/bin/env python3
from pathlib import Path
import csv
import re
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path("results")
LPS_TERMS = Path("results/lps_analysis/lps_go_term_category_map.tsv")
TREE_META = Path("phylogenomic_context_phylophlan/phylophlan_ml_tree_order_with_focal_genes.tsv")
OUTDIR = Path("figure3_go_context_92")
FIGDIR = Path("figures")


def read_terms(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["go_id"], df["go_name"]))


def read_term_categories(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["go_id"], df["category"]))


def category_for_go_name(name: str) -> str:
    low = str(name).lower()
    if "lipid a" in low or "lipid-a" in low:
        return "Lipid A biosynthesis"
    if "o-antigen" in low or "o antigen" in low or "oantigen" in low:
        return "O-antigen pathways"
    if "core oligosaccharide" in low or "lipopolysaccharide core" in low:
        return "Core oligosaccharide assembly"
    return "Other LPS-associated GO terms"


def parse_emapper_go(path: Path, target_gos: set[str]) -> dict[str, set[str]]:
    query_idx = None
    gos_idx = None
    gene_to_gos = defaultdict(set)
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                query_idx = header.index("query") if "query" in header else None
                gos_idx = header.index("GOs") if "GOs" in header else None
                continue
            if line.startswith("#") or query_idx is None or gos_idx is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(query_idx, gos_idx):
                continue
            gene = parts[query_idx].strip()
            gos_field = parts[gos_idx].strip()
            if not gene or not gos_field or gos_field == "-":
                continue
            for go in re.split(r"[,;|]", gos_field):
                go = go.strip()
                if go in target_gos:
                    gene_to_gos[gene].add(go)
    return gene_to_gos


def main() -> None:
    OUTDIR.mkdir(exist_ok=True)
    FIGDIR.mkdir(exist_ok=True)

    meta = pd.read_csv(TREE_META, sep="\t")
    meta = meta[["genome", "source", "species"]].drop_duplicates()
    if meta["genome"].nunique() != 92:
        raise RuntimeError(f"Expected 92 genomes, found {meta['genome'].nunique()}")

    terms = read_terms(LPS_TERMS)
    go_to_category = read_term_categories(LPS_TERMS)
    categories = [
        "Lipid A biosynthesis",
        "Core oligosaccharide assembly",
        "O-antigen pathways",
        "LPS modification enzymes",
    ]

    rows = []
    missing = []
    for rec in meta.itertuples(index=False):
        ann = ROOT / rec.genome / "eggnog_out" / f"{rec.genome}.emapper.annotations"
        if not ann.exists():
            missing.append(str(ann))
            continue
        gene_to_gos = parse_emapper_go(ann, set(terms))
        gene_to_categories = defaultdict(set)
        for gene, gos in gene_to_gos.items():
            for go in gos:
                gene_to_categories[go_to_category[go]].add(gene)
        row = {
            "genome": rec.genome,
            "source": rec.source,
            "species": rec.species,
            "total_unique_lps_go_genes": len(gene_to_gos),
        }
        for cat in categories:
            row[cat] = len(gene_to_categories.get(cat, set()))
        rows.append(row)

    if missing:
        raise RuntimeError("Missing annotation files:\n" + "\n".join(missing[:20]))

    counts = pd.DataFrame(rows).sort_values(["source", "species", "genome"])
    counts.to_csv(OUTDIR / "figure3_lps_go_category_counts_92.tsv", sep="\t", index=False)

    long = counts.melt(
        id_vars=["genome", "source", "species"],
        value_vars=categories,
        var_name="category",
        value_name="genes",
    )
    summary = (
        long.groupby(["source", "category"], as_index=False)
        .agg(mean=("genes", "mean"), sem=("genes", "sem"), n=("genome", "nunique"))
    )
    summary.to_csv(OUTDIR / "figure3_lps_go_category_summary_92.tsv", sep="\t", index=False)

    palette = {"Earth": "#4C78A8", "ISS": "#7B61A9"}
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.0)
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.7), gridspec_kw={"width_ratios": [1.25, 1.0]})

    cat_order = [
        "LPS modification enzymes",
        "Core oligosaccharide assembly",
        "Lipid A biosynthesis",
    ]
    label_map = {
        "LPS modification enzymes": "LPS modification\nenzymes",
        "Core oligosaccharide assembly": "Core oligosaccharide\nassembly",
        "Lipid A biosynthesis": "Lipid A\nbiosynthesis",
    }

    ax = axes[0]
    sns.barplot(
        data=long,
        y="category",
        x="genes",
        hue="source",
        order=cat_order,
        hue_order=["Earth", "ISS"],
        palette=palette,
        errorbar="se",
        ax=ax,
    )
    ax.set_yticklabels([label_map[x] for x in cat_order])
    ax.set_xlabel("Genes per genome with LPS-related GO terms")
    ax.set_ylabel("")
    ax.text(0.0, 1.02, "A", transform=ax.transAxes, ha="left", va="bottom", fontweight="bold")
    ax.legend(frameon=False, title="", loc="lower right")
    ax.grid(axis="y", visible=False)

    ax2 = axes[1]
    sns.boxplot(
        data=counts,
        x="source",
        y="total_unique_lps_go_genes",
        order=["Earth", "ISS"],
        hue="source",
        palette=palette,
        width=0.55,
        fliersize=0,
        legend=False,
        ax=ax2,
    )
    sns.stripplot(
        data=counts,
        x="source",
        y="total_unique_lps_go_genes",
        order=["Earth", "ISS"],
        color="#222222",
        size=2.3,
        alpha=0.5,
        jitter=0.18,
        ax=ax2,
    )
    n_by_source = counts.groupby("source")["genome"].nunique().to_dict()
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels([f"Earth\nn={n_by_source.get('Earth', 0)}", f"ISS\nn={n_by_source.get('ISS', 0)}"])
    ax2.set_xlabel("")
    ax2.set_ylabel("Total genes per genome with LPS-related GO terms")
    ax2.text(0.0, 1.02, "B", transform=ax2.transAxes, ha="left", va="bottom", fontweight="bold")
    ax2.grid(axis="x", visible=False)

    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(FIGDIR / f"figure3_lps_go_category_context.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)

    print(counts.groupby("source")["genome"].nunique().to_string())
    print(FIGDIR / "figure3_lps_go_category_context.png")
    print(OUTDIR / "figure3_lps_go_category_counts_92.tsv")


if __name__ == "__main__":
    main()
