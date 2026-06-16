#!/usr/bin/env python3
import argparse
import csv
import os
import re
from collections import defaultdict


CATEGORIES = [
    "Lipid A biosynthesis",
    "Core oligosaccharide assembly",
    "O-antigen pathways",
    "LPS modification enzymes",
]


def classify_gene(gene_id: str) -> str:
    g = gene_id.strip().lower()
    if not g:
        return "LPS modification enzymes"

    # Lipid A biosynthesis core enzymes
    if re.match(r"^lpx[a-z0-9]*$", g):
        # lpxT is a lipid A modification enzyme, not core biosynthesis
        if g == "lpxt":
            return "LPS modification enzymes"
        return "Lipid A biosynthesis"

    # Core oligosaccharide assembly
    if re.match(r"^(rfa|waa)", g):
        return "Core oligosaccharide assembly"
    if re.match(r"^(kds|kdt)", g):
        return "Core oligosaccharide assembly"
    if re.match(r"^(hld|gmh)", g):
        return "Core oligosaccharide assembly"

    # O-antigen pathways
    if re.match(r"^(rfb|wbb|wzm|wzt|wzx|wzy|wbp|wzz)", g):
        return "O-antigen pathways"

    # LPS modification / transport / regulation
    if re.match(r"^(ept|arn|pmr|pag)", g):
        return "LPS modification enzymes"
    if re.match(r"^(lpt|msb)", g):
        return "LPS modification enzymes"
    if g in {"rfah", "upps", "diaa"}:
        return "LPS modification enzymes"

    return "LPS modification enzymes"


def load_gene_list(path):
    genes = []
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return genes
        gene_idx = header.index("gene") if "gene" in header else 0
        for row in reader:
            if len(row) <= gene_idx:
                continue
            gene = row[gene_idx].strip()
            if gene:
                genes.append(gene)
    return sorted(set(genes))


def load_presence_long(path):
    rows = []
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return rows
        isolate_idx = header.index("isolate")
        group_idx = header.index("group")
        gene_idx = header.index("gene_id")
        present_idx = header.index("present")
        for row in reader:
            if len(row) <= max(isolate_idx, group_idx, gene_idx, present_idx):
                continue
            rows.append(
                {
                    "isolate": row[isolate_idx].strip(),
                    "group": row[group_idx].strip(),
                    "gene_id": row[gene_idx].strip(),
                    "present": row[present_idx].strip(),
                }
            )
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Build per-isolate LPS category counts using KO-annotated LPS gene list."
    )
    parser.add_argument("--gene-list", required=True, help="TSV with LPS genes (e.g., lps_prevalence_diff.tsv)")
    parser.add_argument("--presence-long", required=True, help="Per-isolate presence table (lps_presence_long.tsv)")
    parser.add_argument("--out-mapping", required=True, help="Output mapping TSV (gene -> category)")
    parser.add_argument("--out-table", required=True, help="Output per-isolate table")
    parser.add_argument("--out-summary", required=True, help="Output group summary table")
    args = parser.parse_args()

    genes = load_gene_list(args.gene_list)
    gene_to_cat = {g: classify_gene(g) for g in genes}

    with open(args.out_mapping, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["gene_id", "category"])
        for g in sorted(gene_to_cat.keys()):
            writer.writerow([g, gene_to_cat[g]])

    presence_rows = load_presence_long(args.presence_long)
    counts = defaultdict(lambda: {cat: 0 for cat in CATEGORIES})
    totals = defaultdict(int)
    group_map = {}

    gene_set = set(genes)
    for row in presence_rows:
        gene = row["gene_id"]
        if gene not in gene_set:
            continue
        if row["present"] != "1":
            continue
        isolate = row["isolate"]
        group = row["group"]
        group_map[isolate] = group
        cat = gene_to_cat.get(gene, "LPS modification enzymes")
        counts[isolate][cat] += 1
        totals[isolate] += 1

    isolates = sorted(counts.keys(), key=lambda i: (group_map.get(i, ""), i))

    with open(args.out_table, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            [
                "isolate",
                "group",
                "lipid_a_genes",
                "core_genes",
                "o_antigen_genes",
                "lps_modification_genes",
                "total_unique_genes",
            ]
        )
        for iso in isolates:
            writer.writerow(
                [
                    iso,
                    group_map.get(iso, ""),
                    counts[iso]["Lipid A biosynthesis"],
                    counts[iso]["Core oligosaccharide assembly"],
                    counts[iso]["O-antigen pathways"],
                    counts[iso]["LPS modification enzymes"],
                    totals[iso],
                ]
            )

    summary = defaultdict(lambda: {cat: 0 for cat in CATEGORIES})
    group_sizes = defaultdict(int)
    for iso in isolates:
        grp = group_map.get(iso, "Unknown")
        group_sizes[grp] += 1
        for cat in CATEGORIES:
            summary[grp][cat] += counts[iso][cat]

    with open(args.out_summary, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(
            [
                "group",
                "n_isolates",
                "mean_lipid_a_genes",
                "mean_core_genes",
                "mean_o_antigen_genes",
                "mean_lps_modification_genes",
                "mean_total_unique_genes",
            ]
        )
        for grp in sorted(summary.keys()):
            n = group_sizes[grp]
            writer.writerow(
                [
                    grp,
                    n,
                    summary[grp]["Lipid A biosynthesis"] / n if n else 0,
                    summary[grp]["Core oligosaccharide assembly"] / n if n else 0,
                    summary[grp]["O-antigen pathways"] / n if n else 0,
                    summary[grp]["LPS modification enzymes"] / n if n else 0,
                    (summary[grp]["Lipid A biosynthesis"]
                     + summary[grp]["Core oligosaccharide assembly"]
                     + summary[grp]["O-antigen pathways"]
                     + summary[grp]["LPS modification enzymes"]) / n if n else 0,
                ]
            )


if __name__ == "__main__":
    main()
