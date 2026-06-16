#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re
from collections import defaultdict


DEFAULT_CATEGORIES = [
    "Lipid A biosynthesis",
    "Core oligosaccharide assembly",
    "O-antigen pathways",
    "LPS modification enzymes",
]


def read_lps_terms(path):
    terms = {}
    with open(path, "r", newline="") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            go_id = parts[0].strip() if parts else ""
            if not go_id.startswith("GO:"):
                continue
            name = parts[1].strip() if len(parts) > 1 else ""
            terms[go_id] = name
    return terms


def suggest_category(go_name):
    name = go_name.lower()
    if "lipid a" in name or "lipid-a" in name:
        return "Lipid A biosynthesis"
    if "o-antigen" in name or "o antigen" in name or "oantigen" in name:
        return "O-antigen pathways"
    if "core oligosaccharide" in name or "lipopolysaccharide core" in name or "core" in name:
        return "Core oligosaccharide assembly"
    return "LPS modification enzymes"


def write_mapping(path, terms):
    with open(path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["go_id", "go_name", "category"])
        for go_id in sorted(terms.keys()):
            writer.writerow([go_id, terms[go_id], suggest_category(terms[go_id])])


def load_mapping(path):
    mapping = {}
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return mapping
        try:
            go_idx = header.index("go_id")
            cat_idx = header.index("category")
        except ValueError:
            raise ValueError("Mapping file must include go_id and category columns.")
        for row in reader:
            if len(row) <= max(go_idx, cat_idx):
                continue
            go_id = row[go_idx].strip()
            cat = row[cat_idx].strip()
            if go_id and cat:
                mapping[go_id] = cat
    return mapping


def load_group_map(path, isolate_col="isolate", group_col="group"):
    group_map = {}
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return group_map
        isolate_idx = header.index(isolate_col)
        group_idx = header.index(group_col)
        for row in reader:
            if len(row) <= max(isolate_idx, group_idx):
                continue
            isolate = row[isolate_idx].strip()
            group = row[group_idx].strip()
            if isolate and group:
                group_map[isolate] = group
    return group_map


def find_annotation_files(root):
    pattern = os.path.join(root, "*", "eggnog_out", "*.emapper.annotations")
    return glob.glob(pattern)


def find_go_column(header):
    candidates = ["GOs", "GO_terms", "GO", "GOterm", "GO_terms"]
    for c in candidates:
        if c in header:
            return header.index(c)
    for i, col in enumerate(header):
        if col.lower().startswith("go"):
            return i
    return None


def parse_emapper_go(path, target_go_terms):
    gos_idx = None
    query_idx = None
    gene_to_go = defaultdict(set)
    with open(path, "r", newline="") as handle:
        for line in handle:
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                try:
                    query_idx = header.index("query")
                except ValueError:
                    return gene_to_go
                gos_idx = find_go_column(header)
                continue
            if line.startswith("#"):
                continue
            if gos_idx is None or query_idx is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(query_idx, gos_idx):
                continue
            query = parts[query_idx].strip()
            if not query:
                continue
            gos_field = parts[gos_idx].strip()
            if not gos_field or gos_field == "-":
                continue
            for go in re.split(r"[;,|]", gos_field):
                go = go.strip()
                if go in target_go_terms:
                    gene_to_go[query].add(go)
    return gene_to_go


def main():
    parser = argparse.ArgumentParser(
        description="Build per-genome LPS GO-term category table from eggNOG annotations."
    )
    parser.add_argument("--annotations-root", required=True, help="Root directory with genome folders")
    parser.add_argument("--lps-terms", required=True, help="TSV with GO IDs (LPS_terms.tsv)")
    parser.add_argument("--group-map", required=True, help="TSV with isolate and group columns")
    parser.add_argument("--mapping-out", required=True, help="Output mapping TSV (GO -> category)")
    parser.add_argument("--out-table", required=True, help="Output per-isolate table")
    parser.add_argument("--out-summary", required=True, help="Output group summary table")
    args = parser.parse_args()

    terms = read_lps_terms(args.lps_terms)
    if not terms:
        raise ValueError("No GO terms found in LPS_terms.tsv.")
    write_mapping(args.mapping_out, terms)
    go_to_cat = load_mapping(args.mapping_out)

    group_map = load_group_map(args.group_map)
    files = find_annotation_files(args.annotations_root)

    rows = []
    for path in files:
        isolate = os.path.basename(path).split(".emapper.annotations")[0]
        group = group_map.get(isolate, "")
        gene_to_go = parse_emapper_go(path, set(terms.keys()))
        gene_to_cat = defaultdict(set)
        for gene, gos in gene_to_go.items():
            for go in gos:
                cat = go_to_cat.get(go, "LPS modification enzymes")
                gene_to_cat[cat].add(gene)
        counts = {cat: len(gene_to_cat.get(cat, set())) for cat in DEFAULT_CATEGORIES}
        total_genes = len(gene_to_go)
        rows.append(
            {
                "isolate": isolate,
                "group": group,
                "lipid_a_genes": counts["Lipid A biosynthesis"],
                "core_genes": counts["Core oligosaccharide assembly"],
                "o_antigen_genes": counts["O-antigen pathways"],
                "lps_modification_genes": counts["LPS modification enzymes"],
                "total_unique_genes": total_genes,
            }
        )

    rows.sort(key=lambda r: (r["group"], r["isolate"]))

    with open(args.out_table, "w", newline="") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=[
                "isolate",
                "group",
                "lipid_a_genes",
                "core_genes",
                "o_antigen_genes",
                "lps_modification_genes",
                "total_unique_genes",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    summary = defaultdict(lambda: defaultdict(int))
    counts_per_group = defaultdict(int)
    for row in rows:
        grp = row["group"] or "Unknown"
        counts_per_group[grp] += 1
        for key in [
            "lipid_a_genes",
            "core_genes",
            "o_antigen_genes",
            "lps_modification_genes",
            "total_unique_genes",
        ]:
            summary[grp][key] += row[key]

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
            n = counts_per_group[grp]
            writer.writerow(
                [
                    grp,
                    n,
                    summary[grp]["lipid_a_genes"] / n if n else 0,
                    summary[grp]["core_genes"] / n if n else 0,
                    summary[grp]["o_antigen_genes"] / n if n else 0,
                    summary[grp]["lps_modification_genes"] / n if n else 0,
                    summary[grp]["total_unique_genes"] / n if n else 0,
                ]
            )


if __name__ == "__main__":
    main()
