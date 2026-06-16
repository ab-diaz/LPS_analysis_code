#!/usr/bin/env python3
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
DATA_ROOT = Path("results")
OUTDIR = ROOT / "pantoea_route1_strengthening"
LOCATIONS = OUTDIR / "target_gene_bakta_locations.tsv"
ARCH = OUTDIR / "oantigen_locus_architecture_summary.tsv"
MATRIX = ROOT / "phylogeny_aware_inputs" / "focal_gene_matrix.tsv"
LPS_HITS = Path("results/lps_analysis/lps_eggnog_hits.tsv")
OUT_TSV = OUTDIR / "oantigen_locus_validation_table.tsv"
OUT_CSV = OUTDIR / "oantigen_locus_validation_table.csv"
OUT_BY_GENOME_GENE_TSV = OUTDIR / "oantigen_locus_validation_by_genome_gene.tsv"
OUT_BY_GENOME_GENE_CSV = OUTDIR / "oantigen_locus_validation_by_genome_gene.csv"
OUT_SUMMARY = OUTDIR / "oantigen_locus_validation_table_summary.md"

TARGET_ORDER = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL", "wzm", "wzt", "rfaZ"]
TARGET_LOWER = [g.lower() for g in TARGET_ORDER]


def read_emapper_annotations(isolate, wanted_loci):
    path = DATA_ROOT / isolate / "eggnog_out" / f"{isolate}.emapper.annotations"
    if not path.exists():
        return {}
    header = None
    rows = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("##") or not line.strip():
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#") or header is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            rec = dict(zip(header, parts))
            query = rec.get("query", "")
            if query in wanted_loci:
                rows[query] = rec
    return rows


def clean_value(value):
    if pd.isna(value):
        return ""
    value = str(value)
    return "" if value in {"nan", "-"} else value


def nearest_target_genes(row, locations):
    sub = locations[
        locations["isolate"].eq(row["isolate"])
        & locations["contig"].eq(row["contig"])
        & ~locations["locus_tag"].eq(row["locus_tag"])
    ].copy()
    if sub.empty:
        return ""

    anchor_mid = (int(row["start"]) + int(row["end"])) / 2
    sub["mid"] = (sub["start"].astype(int) + sub["end"].astype(int)) / 2
    sub["distance_bp"] = (sub["mid"] - anchor_mid).abs().astype(int)
    sub = sub.sort_values(["distance_bp", "start"]).head(4)
    return ";".join(
        f"{r.gene}({r.locus_tag},{int(r.distance_bp)}bp)"
        for r in sub.itertuples(index=False)
    )


def main():
    loc = pd.read_csv(LOCATIONS, sep="\t")
    arch = pd.read_csv(ARCH, sep="\t")
    matrix = pd.read_csv(MATRIX, sep="\t")
    lps_hits = pd.read_csv(LPS_HITS, sep="\t") if LPS_HITS.exists() else pd.DataFrame()

    loc["gene_lower"] = loc["gene"].str.lower()
    loc = loc[loc["gene_lower"].isin(TARGET_LOWER)].copy()
    loc["gene_standard"] = loc["gene_lower"].map(dict(zip(TARGET_LOWER, TARGET_ORDER)))

    meta_cols = [
        "genome",
        "source",
        "year",
        "collection_date",
        "year_source",
        "species",
        "target_profile_cluster",
        "strict_lps_profile_cluster",
    ]
    meta = matrix[meta_cols].rename(columns={"genome": "isolate"})
    arch_cols = [
        "isolate",
        "complete_rfbABCD",
        "rfbABCD_colocalized_15kb",
        "rfb_contig",
        "rfb_order_by_coordinate",
        "rfb_span_bp",
        "wzm_wzt_colocalized_10kb",
    ]
    arch = arch[arch_cols]

    wanted_by_isolate = loc.groupby("isolate")["locus_tag"].apply(lambda s: set(s)).to_dict()
    emapper_rows = {}
    for isolate, wanted in wanted_by_isolate.items():
        emapper_rows[isolate] = read_emapper_annotations(isolate, wanted)

    lps_lookup = {}
    if not lps_hits.empty:
        for r in lps_hits.itertuples(index=False):
            lps_lookup[(r.isolate, r.query)] = {
                "preferred": clean_value(getattr(r, "preferred", "")),
                "description": clean_value(getattr(r, "description", "")),
                "kegg_ko": clean_value(getattr(r, "kegg_ko", "")),
                "kegg_pathway": clean_value(getattr(r, "kegg_pathway", "")),
            }

    records = []
    for row in loc.sort_values(["isolate", "contig", "start", "end"]).itertuples(index=False):
        base = row._asdict()
        isolate = base["isolate"]
        locus_tag = base["locus_tag"]
        ann = emapper_rows.get(isolate, {}).get(locus_tag, {})
        fallback = lps_lookup.get((isolate, locus_tag), {})
        preferred = clean_value(ann.get("Preferred_name", "")) or fallback.get("preferred", "")
        description = clean_value(ann.get("Description", "")) or fallback.get("description", "") or clean_value(base.get("eggnog_description", ""))
        kegg_ko = clean_value(ann.get("KEGG_ko", "")) or fallback.get("kegg_ko", "")
        kegg_pathway = clean_value(ann.get("KEGG_Pathway", "")) or fallback.get("kegg_pathway", "")

        records.append(
            {
                "isolate": isolate,
                "source": "",
                "year": "",
                "collection_date": "",
                "year_source": "",
                "species": "",
                "gene": base["gene_standard"],
                "gene_call_used_for_detection": base["gene"],
                "bakta_locus_tag": locus_tag,
                "contig": base["contig"],
                "start": int(base["start"]),
                "end": int(base["end"]),
                "strand": base["strand"],
                "bakta_product": clean_value(base["product"]),
                "eggnog_preferred_name": preferred,
                "eggnog_description": description,
                "kegg_ko": kegg_ko,
                "kegg_pathway": kegg_pathway,
                "complete_rfbABCD_in_genome": "",
                "rfbABCD_colocalized_15kb": "",
                "rfb_contig": "",
                "rfb_order_by_coordinate": "",
                "rfb_span_bp": "",
                "hit_belongs_to_compact_rfbABCD_module": "",
                "wzm_wzt_colocalized_10kb": "",
                "nearest_target_genes_same_contig": nearest_target_genes(base, loc),
            }
        )

    table = pd.DataFrame(records)
    table = table.merge(meta, on="isolate", how="left", suffixes=("", "_meta"))
    for col in ["source", "year", "collection_date", "year_source", "species"]:
        table[col] = table[f"{col}_meta"].combine_first(table[col])
        table = table.drop(columns=[f"{col}_meta"])
    table = table.merge(arch, on="isolate", how="left", suffixes=("", "_arch"))
    for col in [
        "complete_rfbABCD",
        "rfbABCD_colocalized_15kb",
        "rfb_contig",
        "rfb_order_by_coordinate",
        "rfb_span_bp",
        "wzm_wzt_colocalized_10kb",
    ]:
        target = {
            "complete_rfbABCD": "complete_rfbABCD_in_genome",
            "rfbABCD_colocalized_15kb": "rfbABCD_colocalized_15kb",
            "rfb_contig": "rfb_contig",
            "rfb_order_by_coordinate": "rfb_order_by_coordinate",
            "rfb_span_bp": "rfb_span_bp",
            "wzm_wzt_colocalized_10kb": "wzm_wzt_colocalized_10kb",
        }[col]
        if f"{col}_arch" in table.columns:
            table[target] = table[f"{col}_arch"].combine_first(table[target])
            table = table.drop(columns=[f"{col}_arch"])
        elif col in table.columns and target != col:
            table[target] = table[col].combine_first(table[target])
            table = table.drop(columns=[col])

    table["hit_belongs_to_compact_rfbABCD_module"] = (
        table["gene"].isin(["rfbA", "rfbB", "rfbC", "rfbD"])
        & table["rfbABCD_colocalized_15kb"].fillna(0).astype(float).eq(1)
        & table["contig"].eq(table["rfb_contig"])
    ).astype(int)

    ordered_cols = [
        "isolate",
        "source",
        "year",
        "collection_date",
        "year_source",
        "species",
        "gene",
        "gene_call_used_for_detection",
        "bakta_locus_tag",
        "contig",
        "start",
        "end",
        "strand",
        "bakta_product",
        "eggnog_preferred_name",
        "eggnog_description",
        "kegg_ko",
        "kegg_pathway",
        "complete_rfbABCD_in_genome",
        "rfbABCD_colocalized_15kb",
        "rfb_contig",
        "rfb_order_by_coordinate",
        "rfb_span_bp",
        "hit_belongs_to_compact_rfbABCD_module",
        "wzm_wzt_colocalized_10kb",
        "nearest_target_genes_same_contig",
        "target_profile_cluster",
        "strict_lps_profile_cluster",
    ]
    table = table[ordered_cols]
    table = table.fillna("")
    table.to_csv(OUT_TSV, sep="\t", index=False)
    table.to_csv(OUT_CSV, index=False)

    genome_gene_rows = []
    grouped_hits = table.groupby(["isolate", "gene"], dropna=False)
    arch_by_iso = arch.set_index("isolate").to_dict(orient="index")
    for meta_row in matrix.sort_values(["source", "species", "genome"]).itertuples(index=False):
        isolate = meta_row.genome
        arch_row = arch_by_iso.get(isolate, {})
        for gene in TARGET_ORDER:
            hits = grouped_hits.get_group((isolate, gene)) if (isolate, gene) in grouped_hits.groups else pd.DataFrame()
            present = int(getattr(meta_row, gene) == 1)
            genome_gene_rows.append(
                {
                    "isolate": isolate,
                    "source": clean_value(meta_row.source),
                    "year": clean_value(meta_row.year),
                    "collection_date": clean_value(meta_row.collection_date),
                    "year_source": clean_value(meta_row.year_source),
                    "species": clean_value(meta_row.species),
                    "gene": gene,
                    "present": present,
                    "mapped_hit_count": int(len(hits)),
                    "bakta_locus_tags": ";".join(hits["bakta_locus_tag"].astype(str)) if len(hits) else "",
                    "contigs": ";".join(dict.fromkeys(hits["contig"].astype(str))) if len(hits) else "",
                    "coordinate_ranges": ";".join(
                        f"{r.contig}:{int(r.start)}-{int(r.end)}({r.strand})"
                        for r in hits.itertuples(index=False)
                    )
                    if len(hits)
                    else "",
                    "bakta_products": ";".join(dict.fromkeys(hits["bakta_product"].astype(str))) if len(hits) else "",
                    "eggnog_preferred_names": ";".join(dict.fromkeys(hits["eggnog_preferred_name"].astype(str))) if len(hits) else "",
                    "eggnog_descriptions": ";".join(dict.fromkeys(hits["eggnog_description"].astype(str))) if len(hits) else "",
                    "kegg_kos": ";".join(dict.fromkeys(hits["kegg_ko"].astype(str))) if len(hits) else "",
                    "kegg_pathways": ";".join(dict.fromkeys(hits["kegg_pathway"].astype(str))) if len(hits) else "",
                    "complete_rfbABCD_in_genome": clean_value(arch_row.get("complete_rfbABCD", "")),
                    "rfbABCD_colocalized_15kb": clean_value(arch_row.get("rfbABCD_colocalized_15kb", "")),
                    "rfb_contig": clean_value(arch_row.get("rfb_contig", "")),
                    "rfb_order_by_coordinate": clean_value(arch_row.get("rfb_order_by_coordinate", "")),
                    "rfb_span_bp": clean_value(arch_row.get("rfb_span_bp", "")),
                    "belongs_to_compact_rfbABCD_module": int(
                        gene in {"rfbA", "rfbB", "rfbC", "rfbD"}
                        and present
                        and str(clean_value(arch_row.get("rfbABCD_colocalized_15kb", ""))) in {"1", "1.0"}
                    ),
                    "wzm_wzt_colocalized_10kb": clean_value(arch_row.get("wzm_wzt_colocalized_10kb", "")),
                    "nearest_target_genes_same_contig": ";".join(
                        dict.fromkeys(hits["nearest_target_genes_same_contig"].astype(str))
                    )
                    if len(hits)
                    else "",
                    "target_profile_cluster": clean_value(meta_row.target_profile_cluster),
                    "strict_lps_profile_cluster": clean_value(meta_row.strict_lps_profile_cluster),
                }
            )

    genome_gene_table = pd.DataFrame(genome_gene_rows).fillna("")
    genome_gene_table.to_csv(OUT_BY_GENOME_GENE_TSV, sep="\t", index=False)
    genome_gene_table.to_csv(OUT_BY_GENOME_GENE_CSV, index=False)

    by_gene = table.groupby(["gene", "source"]).size().unstack(fill_value=0).reindex(TARGET_ORDER).fillna(0).astype(int)
    summary_lines = [
        "# O-antigen/LPS locus validation table",
        "",
        f"- Rows: {len(table)} focal gene hits",
        f"- Isolates with at least one focal hit: {table['isolate'].nunique()}",
        f"- Output TSV: `{OUT_TSV.relative_to(ROOT)}`",
        f"- Output CSV: `{OUT_CSV.relative_to(ROOT)}`",
        f"- Genome-by-gene TSV, including absent calls: `{OUT_BY_GENOME_GENE_TSV.relative_to(ROOT)}`",
        f"- Genome-by-gene CSV, including absent calls: `{OUT_BY_GENOME_GENE_CSV.relative_to(ROOT)}`",
        "",
        "Each row is one mapped focal gene hit. Coordinates and products come from Bakta; preferred names, descriptions, and KEGG fields come from eggNOG-mapper where available.",
        "",
        "## Hit counts by gene and source",
        "",
        "| gene | Earth | ISS |",
        "| --- | ---: | ---: |",
    ]
    for gene, row in by_gene.iterrows():
        summary_lines.append(f"| {gene} | {int(row.get('Earth', 0))} | {int(row.get('ISS', 0))} |")

    compact_count = int(table["hit_belongs_to_compact_rfbABCD_module"].astype(int).sum())
    summary_lines.extend(
        [
            "",
            f"- Hits belonging to a compact co-localized rfbABCD module: {compact_count}",
            f"- Hits with KEGG KO annotation: {(table['kegg_ko'] != '').sum()}",
            f"- Hits with eggNOG preferred name: {(table['eggnog_preferred_name'] != '').sum()}",
            f"- Genome-by-gene rows: {len(genome_gene_table)} ({matrix['genome'].nunique()} genomes x {len(TARGET_ORDER)} genes)",
            "",
        ]
    )
    OUT_SUMMARY.write_text("\n".join(summary_lines))
    print(f"Wrote {OUT_TSV}")
    print(f"Wrote {OUT_CSV}")
    print(f"Wrote {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
