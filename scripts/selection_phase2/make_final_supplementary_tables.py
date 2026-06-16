#!/usr/bin/env python3
"""Build one final TSV file per supplementary table.

Several supplementary tables are backed by heterogeneous analysis outputs. For
those, the final TSV uses a `section` column plus a union of all columns so that
one deliverable file can contain multiple related result blocks without losing
traceability.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path("results/lps_analysis/selection_phase2")
LPS = Path("results/lps_analysis")
OUT = ROOT / "supplementary_tables_final"


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [
        str(c)
        .replace("\u00a0", " ")
        .strip()
        .lower()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("-", "_")
        for c in df.columns
    ]
    return df


def read_table(path: Path, sep: str = "\t") -> pd.DataFrame:
    if path.stat().st_size <= 1:
        return pd.DataFrame([{"note": "No rows detected in source file."}])
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    return clean_columns(df)


def add_section(df: pd.DataFrame, section: str, source_file: str) -> pd.DataFrame:
    df = df.copy()
    df.insert(0, "source_file", source_file)
    df.insert(0, "section", section)
    return df


def concat_sections(sections: list[pd.DataFrame]) -> pd.DataFrame:
    cols = []
    for df in sections:
        for c in df.columns:
            if c not in cols:
                cols.append(c)
    return pd.concat([df.reindex(columns=cols) for df in sections], ignore_index=True)


def write(df: pd.DataFrame, filename: str) -> Path:
    path = OUT / filename
    df.to_csv(path, sep="\t", index=False)
    return path


def table1() -> Path:
    df = read_table(LPS / "lps_go_term_category_map.tsv")
    df.insert(0, "supplementary_table", "Supplementary Table 1")
    return write(df, "Supplementary_Table_1_LPS_GO_terms.tsv")


def table2() -> Path:
    df = pd.read_csv(LPS / "SupplTbl2.csv", sep=";", dtype=str, keep_default_na=False)
    df = clean_columns(df)
    df.insert(0, "supplementary_table", "Supplementary Table 2")
    return write(df, "Supplementary_Table_2_O_antigen_LPS_keywords.tsv")


def table3() -> Path:
    meta = read_table(LPS / "metadata_with_ncbi.tsv")
    tree = read_table(ROOT / "phylogenomic_context_phylophlan" / "phylophlan_ml_tree_order_with_focal_genes.tsv")
    tree_small = tree[
        [
            "genome",
            "source",
            "year",
            "species",
            "rfba",
            "rfbb",
            "rfbc",
            "rfbd",
            "waal",
            "wzm",
            "wzt",
            "rfaz",
            "strict_lps_profile_cluster",
            "target_profile_cluster",
        ]
    ].rename(columns={"genome": "sample_id", "source": "analysis_source_group", "year": "iss_sampling_year"})
    merged = meta.merge(tree_small, on="sample_id", how="left")
    merged.insert(0, "supplementary_table", "Supplementary Table 3")
    return write(merged, "Supplementary_Table_3_Pantoea_genome_metadata.tsv")


def table4() -> Path:
    files = [
        ("key_result_summary", "pantoea_route1_strengthening/route1_key_results.tsv"),
        ("cluster_collapsed_strict_lps_fisher", "pantoea_route1_strengthening/cluster_collapsed_strict_lps_fisher.tsv"),
        ("leave_one_iss_year_out_fisher", "pantoea_route1_strengthening/leave_one_iss_year_out_fisher.tsv"),
        ("iss_year_persistence", "pantoea_route1_strengthening/iss_year_persistence_target_genes.tsv"),
        ("phylogeny_aware_association", "phylogeny_aware_association/phylogeny_aware_gene_association.tsv"),
    ]
    sections = [add_section(read_table(ROOT / rel), section, rel) for section, rel in files]
    return write(concat_sections(sections), "Supplementary_Table_4_association_founder_phylogeny.tsv")


def table5() -> Path:
    files = [
        ("phylophlan_tree_order_with_focal_genes", "phylogenomic_context_phylophlan/phylophlan_ml_tree_order_with_focal_genes.tsv"),
        ("lps_go_category_counts_by_isolate", "../lps_go_category_counts.tsv"),
        ("lps_go_category_summary_by_group", "../lps_go_category_summary.tsv"),
    ]
    sections = []
    for section, rel in files:
        path = (ROOT / rel).resolve() if rel.startswith("..") else ROOT / rel
        sections.append(add_section(read_table(path), section, rel))
    return write(concat_sections(sections), "Supplementary_Table_5_phylogeny_GO_context.tsv")


def table6() -> Path:
    files = [
        ("by_genome_gene_locus_validation", "pantoea_route1_strengthening/oantigen_locus_validation_by_genome_gene.tsv"),
        ("mapped_hit_locus_validation", "pantoea_route1_strengthening/oantigen_locus_validation_table.tsv"),
        ("surface_glycan_functional_annotation", "pantoea_route1_strengthening/surface_glycan_locus_functional_annotation.tsv"),
        ("locus_evolutionary_context", "pantoea_route1_strengthening/locus_evolutionary_context_by_region.tsv"),
        ("iss_locus_flanking_conservation", "pantoea_route1_strengthening/iss_locus_flanking_conservation.tsv"),
        ("locus_mobility_gene_hits", "pantoea_route1_strengthening/locus_mobility_gene_hits.tsv"),
        ("architecture_summary", "pantoea_route1_strengthening/oantigen_locus_architecture_summary.tsv"),
        ("figure5_representative_rationale", "pantoea_route1_strengthening/figure5_representative_genome_rationale.tsv"),
    ]
    sections = [add_section(read_table(ROOT / rel), section, rel) for section, rel in files]
    return write(concat_sections(sections), "Supplementary_Table_6_locus_validation_context.tsv")


def pfulva_presence_long(path: Path) -> pd.DataFrame:
    wide = read_table(path)
    return wide.melt(id_vars=["gene"], var_name="isolate", value_name="present")


def summary_md_rows(path: Path) -> pd.DataFrame:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append({"summary_line": line})
    return pd.DataFrame(rows)


def table7() -> Path:
    files = [
        ("pfulva_metadata", "pseudomonas_fulva_lps_analysis/genome_metadata_pfulva_expanded_with_mt2_iss.tsv", None),
        ("pfulva_o_antigen_presence_long", "pseudomonas_fulva_lps_analysis/pfulva_expanded_o_antigen_gene_presence.tsv", pfulva_presence_long),
        ("pfulva_lps_category_counts", "pseudomonas_fulva_lps_analysis/pfulva_expanded_lps_category_counts.tsv", None),
        ("pfulva_o_antigen_hits", "pseudomonas_fulva_lps_analysis/pfulva_expanded_o_antigen_hits.tsv", None),
        ("pfulva_results_summary", "pseudomonas_fulva_lps_analysis/pfulva_expanded_results_summary.md", summary_md_rows),
    ]
    sections = []
    for section, rel, loader in files:
        path = ROOT / rel
        df = loader(path) if loader else read_table(path)
        sections.append(add_section(df, section, rel))
    return write(concat_sections(sections), "Supplementary_Table_7_Pseudomonas_fulva_exploratory.tsv")


def main() -> None:
    OUT.mkdir(exist_ok=True)
    paths = [table1(), table2(), table3(), table4(), table5(), table6(), table7()]
    manifest = OUT / "README.md"
    manifest.write_text(
        "# Final supplementary table files\n\n"
        + "\n".join(f"- `{p.name}`" for p in paths)
        + "\n\nNotes:\n"
        "- Tables with heterogeneous outputs use `section` and `source_file` columns.\n"
        "- Supplementary Table 3 is built from the available 92-genome metadata table and PhyloPhlAn focal-gene/source annotations. It does not contain genome size/GC/CDS fields for every genome, so the manuscript description should not promise those fields unless a richer metadata file is supplied.\n",
        encoding="utf-8",
    )
    for p in paths:
        print(p)
    print(manifest)


if __name__ == "__main__":
    main()
