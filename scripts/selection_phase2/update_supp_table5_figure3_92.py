#!/usr/bin/env python3
from pathlib import Path

import pandas as pd


SUPP = Path("supplementary_tables_final/Supplementary_Table_5_phylogeny_GO_context.tsv")
COUNTS = Path("figure3_go_context_92/figure3_lps_go_category_counts_92.tsv")
SUMMARY = Path("figure3_go_context_92/figure3_lps_go_category_summary_92.tsv")


def blank_row(columns):
    return {c: "" for c in columns}


def main() -> None:
    supp = pd.read_csv(SUPP, sep="\t", dtype=str).fillna("")
    tree = supp[supp["section"] == "phylophlan_tree_order_with_focal_genes"].copy()
    columns = supp.columns.tolist()

    counts = pd.read_csv(COUNTS, sep="\t")
    count_rows = []
    for _, rec in counts.iterrows():
        row = blank_row(columns)
        row.update(
            {
                "section": "lps_go_category_counts_by_genome_expanded_92",
                "source_file": str(COUNTS),
                "isolate": rec["genome"],
                "genome": rec["genome"],
                "group": rec["source"],
                "source": rec["source"],
                "species": rec["species"],
                "lipid_a_genes": rec["Lipid A biosynthesis"],
                "core_genes": rec["Core oligosaccharide assembly"],
                "o_antigen_genes": rec["O-antigen pathways"],
                "lps_modification_genes": rec["LPS modification enzymes"],
                "total_unique_genes": rec["total_unique_lps_go_genes"],
            }
        )
        count_rows.append(row)

    summary = pd.read_csv(SUMMARY, sep="\t")
    wide = summary.pivot(index="source", columns="category", values="mean").reset_index()
    n_map = summary.groupby("source")["n"].max().to_dict()
    summary_rows = []
    for _, rec in wide.iterrows():
        row = blank_row(columns)
        src = rec["source"]
        row.update(
            {
                "section": "lps_go_category_summary_by_group_expanded_92",
                "source_file": str(SUMMARY),
                "group": src,
                "source": src,
                "n_isolates": n_map.get(src, ""),
                "mean_lipid_a_genes": rec["Lipid A biosynthesis"],
                "mean_core_genes": rec["Core oligosaccharide assembly"],
                "mean_o_antigen_genes": rec["O-antigen pathways"],
                "mean_lps_modification_genes": rec["LPS modification enzymes"],
                "mean_total_unique_genes": "",
            }
        )
        total = counts.loc[counts["source"] == src, "total_unique_lps_go_genes"].mean()
        row["mean_total_unique_genes"] = total
        summary_rows.append(row)

    out = pd.concat([tree, pd.DataFrame(count_rows), pd.DataFrame(summary_rows)], ignore_index=True)
    out = out[columns]
    out.to_csv(SUPP, sep="\t", index=False)
    print(SUPP)
    print(out["section"].value_counts().to_string())


if __name__ == "__main__":
    main()
