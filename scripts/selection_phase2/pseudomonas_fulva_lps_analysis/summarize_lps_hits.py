#!/usr/bin/env python3
import csv
import os
import re
from collections import defaultdict


ROOT = os.path.dirname(os.path.abspath(__file__))
HITS = os.environ.get("HITS_FILE", os.path.join(ROOT, "lps_eggnog_hits.tsv"))
META = os.environ.get("META_FILE", os.path.join(ROOT, "genome_metadata.tsv"))
OUT_PREFIX = os.environ.get("OUT_PREFIX", "pfulva")


def classify_gene(gene):
    g = (gene or "").strip().lower()

    if re.match(r"^lpx[a-z0-9]*$", g):
        return "Lipid A biosynthesis" if g != "lpxt" else "LPS modification/transport"

    if re.match(r"^(rfa|waa|kds|kdt|hld|gmh)", g):
        return "Core oligosaccharide assembly"

    if re.match(r"^(rfb|wbb|wzm|wzt|wzx|wzy|wbp|wzz)", g):
        return "O-antigen pathways"

    if re.match(r"^(ept|arn|pmr|pag|lpt|msb)", g):
        return "LPS modification/transport"

    if g in {"rfah", "upps", "diaa"}:
        return "LPS modification/transport"

    return "Other LPS-related"


def load_metadata():
    meta = {}
    with open(META, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            meta[row["isolate_id"]] = row
    return meta


def main():
    metadata = load_metadata()
    rows = []

    with open(HITS, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row["preferred"].strip()
            category = classify_gene(gene)
            meta = metadata.get(row["isolate"], {})
            out = {
                "isolate": row["isolate"],
                "group": meta.get("group", ""),
                "source_note": meta.get("source_note", ""),
                "query": row["query"],
                "gene": gene,
                "category": category,
                "description": row["description"],
                "kegg_ko": row["kegg_ko"],
                "kegg_pathway": row["kegg_pathway"],
            }
            rows.append(out)

    categorized_path = os.path.join(ROOT, f"{OUT_PREFIX}_lps_hits_categorized.tsv")
    with open(categorized_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "isolate",
                "group",
                "source_note",
                "query",
                "gene",
                "category",
                "description",
                "kegg_ko",
                "kegg_pathway",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    categories = [
        "Lipid A biosynthesis",
        "Core oligosaccharide assembly",
        "O-antigen pathways",
        "LPS modification/transport",
        "Other LPS-related",
    ]
    counts = defaultdict(lambda: {c: 0 for c in categories})
    genes = defaultdict(set)

    for row in rows:
        iso = row["isolate"]
        gene = row["gene"]
        cat = row["category"]
        genes[(iso, cat)].add(gene)

    for (iso, cat), gene_set in genes.items():
        counts[iso][cat] = len(gene_set)

    counts_path = os.path.join(ROOT, f"{OUT_PREFIX}_lps_category_counts.tsv")
    with open(counts_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "isolate",
                "group",
                "source_note",
                "lipid_a_genes",
                "core_genes",
                "o_antigen_genes",
                "lps_modification_transport_genes",
                "other_lps_related_genes",
                "total_unique_lps_genes",
            ]
        )
        for iso in sorted(metadata):
            total = sum(counts[iso].values())
            writer.writerow(
                [
                    iso,
                    metadata[iso]["group"],
                    metadata[iso]["source_note"],
                    counts[iso]["Lipid A biosynthesis"],
                    counts[iso]["Core oligosaccharide assembly"],
                    counts[iso]["O-antigen pathways"],
                    counts[iso]["LPS modification/transport"],
                    counts[iso]["Other LPS-related"],
                    total,
                ]
            )

    oantigen_rows = [r for r in rows if r["category"] == "O-antigen pathways"]
    o_path = os.path.join(ROOT, f"{OUT_PREFIX}_o_antigen_hits.tsv")
    with open(o_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(oantigen_rows)

    all_o_genes = sorted({r["gene"] for r in oantigen_rows})
    matrix_path = os.path.join(ROOT, f"{OUT_PREFIX}_o_antigen_gene_presence.tsv")
    with open(matrix_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene"] + sorted(metadata))
        for gene in all_o_genes:
            writer.writerow(
                [gene]
                + [
                    "1" if any(r["isolate"] == iso and r["gene"] == gene for r in oantigen_rows) else "0"
                    for iso in sorted(metadata)
                ]
            )

    for path in [categorized_path, counts_path, o_path, matrix_path]:
        print(path)


if __name__ == "__main__":
    main()
