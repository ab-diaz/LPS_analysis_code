#!/usr/bin/env python3
from pathlib import Path
import re

import pandas as pd
from Bio import SeqIO


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "pantoea_route1_strengthening"
GENOME_INDEX = ROOT / "phylogeny_aware_inputs" / "genome_file_index.tsv"
LOCATIONS = OUTDIR / "target_gene_bakta_locations.tsv"
ARCH = OUTDIR / "oantigen_locus_architecture_summary.tsv"
OUT_CONTEXT = OUTDIR / "locus_evolutionary_context_by_region.tsv"
OUT_ISS_FLANKS = OUTDIR / "iss_locus_flanking_conservation.tsv"
OUT_MOBILITY = OUTDIR / "locus_mobility_gene_hits.tsv"
OUT_SUMMARY = OUTDIR / "locus_evolutionary_context_summary.md"

TARGET_LOWER = {"rfba", "rfbb", "rfbc", "rfbd", "waal", "wzm", "wzt", "rfaz"}
PADDING_BP = 10_000
EDGE_WARNING_BP = 2_000
MAX_TARGET_GAP_SAME_REGION_BP = 20_000
MOBILITY_RE = re.compile(
    r"transposase|integrase|recombinase|insertion sequence|\\bIS\\b|phage|prophage|"
    r"plasmid|conjugative|relaxase|mobilization|mobilisation|resolvase|invertase",
    re.IGNORECASE,
)


def clean(value):
    if pd.isna(value):
        return ""
    value = str(value)
    return "" if value.lower() == "nan" else value


def gc_fraction(seq):
    seq = str(seq).upper()
    called = sum(1 for b in seq if b in {"A", "C", "G", "T"})
    if called == 0:
        return None
    return (seq.count("G") + seq.count("C")) / called


def read_bakta_tsv(path):
    header = None
    rows = []
    with Path(path).open(errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("#Sequence Id"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if line.startswith("#"):
                continue
            if header is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            rows.append(dict(zip(header, parts)))
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=["contig", "type", "start", "end", "strand", "locus_tag", "gene", "product"])
    df = df.rename(
        columns={
            "Sequence Id": "contig",
            "Type": "type",
            "Start": "start",
            "Stop": "end",
            "Strand": "strand",
            "Locus Tag": "locus_tag",
            "Gene": "gene",
            "Product": "product",
            "DbXrefs": "dbxrefs",
        }
    )
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    return df


def feature_label(row):
    gene = clean(row.get("gene", ""))
    product = clean(row.get("product", ""))
    locus = clean(row.get("locus_tag", ""))
    if gene:
        return f"{gene}:{product}" if product else gene
    return f"{locus}:{product}" if product else locus


def feature_product_label(row):
    gene = clean(row.get("gene", ""))
    product = clean(row.get("product", ""))
    if gene and product:
        return f"{gene}:{product}"
    return product or gene


def compact_labels(rows, limit=8):
    labels = [feature_label(row) for _, row in rows.iterrows()]
    labels = [x for x in labels if x]
    if len(labels) > limit:
        return ";".join(labels[:limit]) + f";...(+{len(labels)-limit})"
    return ";".join(labels)


def compact_product_labels(rows, limit=8):
    labels = [feature_product_label(row) for _, row in rows.iterrows()]
    labels = [x for x in labels if x]
    if len(labels) > limit:
        return ";".join(labels[:limit]) + f";...(+{len(labels)-limit})"
    return ";".join(labels)


def split_target_clusters(contig_hits, max_gap=MAX_TARGET_GAP_SAME_REGION_BP):
    hits = contig_hits.sort_values(["start", "end"]).copy()
    clusters = []
    current = []
    current_end = None
    for _, hit in hits.iterrows():
        start = int(hit["start"])
        end = int(hit["end"])
        if current and current_end is not None and start - current_end > max_gap:
            clusters.append(pd.DataFrame(current))
            current = []
        current.append(hit)
        current_end = max(current_end or end, end)
    if current:
        clusters.append(pd.DataFrame(current))
    return clusters


def main():
    index = pd.read_csv(GENOME_INDEX, sep="\t")
    loc = pd.read_csv(LOCATIONS, sep="\t")
    arch = pd.read_csv(ARCH, sep="\t")

    loc["gene_lower"] = loc["gene"].str.lower()
    loc = loc[loc["gene_lower"].isin(TARGET_LOWER)].copy()

    arch_by_iso = arch.set_index("isolate").to_dict(orient="index")
    meta_by_iso = index.set_index("genome").to_dict(orient="index")

    context_rows = []
    mobility_rows = []
    iss_flank_rows = []

    for isolate, iso_hits in loc.groupby("isolate"):
        meta = meta_by_iso.get(isolate, {})
        fna = Path(clean(meta.get("fna", "")))
        bakta_tsv = Path(clean(meta.get("bakta_tsv", "")))
        if not fna.exists() or not bakta_tsv.exists():
            continue

        seqs = {rec.id: rec.seq for rec in SeqIO.parse(fna, "fasta")}
        genome_len = sum(len(seq) for seq in seqs.values())
        genome_gc = gc_fraction("".join(str(seq) for seq in seqs.values()))
        features = read_bakta_tsv(bakta_tsv)
        cds = features[features["type"].str.lower().eq("cds")].copy()

        arch_row = arch_by_iso.get(isolate, {})
        rfb_contig = clean(arch_row.get("rfb_contig", ""))

        for contig, all_contig_hits in iso_hits.groupby("contig"):
            if contig not in seqs:
                continue
            seq = seqs[contig]
            contig_len = len(seq)
            target_clusters = split_target_clusters(all_contig_hits)

            for cluster_index, contig_hits in enumerate(target_clusters, start=1):
                locus_start = int(contig_hits["start"].min())
                locus_end = int(contig_hits["end"].max())
                region_start = max(1, locus_start - PADDING_BP)
                region_end = min(contig_len, locus_end + PADDING_BP)

                locus_seq = seq[locus_start - 1 : locus_end]
                region_seq = seq[region_start - 1 : region_end]
                local_gc = gc_fraction(locus_seq)
                region_gc = gc_fraction(region_seq)
                delta_gc = None if genome_gc is None or local_gc is None else local_gc - genome_gc

                contig_cds = cds[cds["contig"].eq(contig)].copy()
                neighborhood = contig_cds[
                    (contig_cds["start"].astype(int) <= region_end)
                    & (contig_cds["end"].astype(int) >= region_start)
                ].copy()
                target_loci = set(contig_hits["locus_tag"])
                flanks = neighborhood[~neighborhood["locus_tag"].isin(target_loci)].copy()
                upstream = flanks[flanks["end"].astype(int) < locus_start].sort_values("end", ascending=False).head(5)
                downstream = flanks[flanks["start"].astype(int) > locus_end].sort_values("start").head(5)
                mobility = neighborhood[
                    neighborhood[["gene", "product"]]
                    .fillna("")
                    .agg(" ".join, axis=1)
                    .str.contains(MOBILITY_RE, na=False)
                ].copy()

                for _, mob in mobility.iterrows():
                    mobility_rows.append(
                        {
                            "isolate": isolate,
                            "source": clean(meta.get("source", "")),
                            "year": clean(meta.get("year", "")),
                            "species": clean(meta.get("species", "")),
                            "contig": contig,
                            "target_cluster_index_on_contig": cluster_index,
                            "region_start": region_start,
                            "region_end": region_end,
                            "mobility_locus_tag": clean(mob.get("locus_tag", "")),
                            "mobility_gene": clean(mob.get("gene", "")),
                            "mobility_product": clean(mob.get("product", "")),
                            "mobility_start": int(mob["start"]),
                            "mobility_end": int(mob["end"]),
                        }
                    )

                genes_detected = ",".join(sorted(contig_hits["gene_lower"].unique()))
                target_gene_set = set(contig_hits["gene_lower"])
                is_primary_rfb_region = int(len(target_gene_set & {"rfba", "rfbb", "rfbc", "rfbd"}) >= 2)
                if rfb_contig and contig != rfb_contig:
                    is_primary_rfb_region = 0

                min_edge = min(locus_start - 1, contig_len - locus_end)
                row = {
                    "isolate": isolate,
                    "source": clean(meta.get("source", "")),
                    "year": clean(meta.get("year", "")),
                    "species": clean(meta.get("species", "")),
                    "contig": contig,
                    "target_cluster_index_on_contig": cluster_index,
                    "target_cluster_max_gap_bp": MAX_TARGET_GAP_SAME_REGION_BP,
                    "contig_length_bp": contig_len,
                    "genome_length_bp": genome_len,
                    "target_genes_in_region": genes_detected,
                    "target_hit_count": len(contig_hits),
                    "locus_start": locus_start,
                    "locus_end": locus_end,
                    "locus_size_bp": locus_end - locus_start + 1,
                    "padded_region_start": region_start,
                    "padded_region_end": region_end,
                    "padded_region_size_bp": region_end - region_start + 1,
                    "genome_gc": genome_gc,
                    "locus_gc": local_gc,
                    "padded_region_gc": region_gc,
                    "delta_locus_gc_vs_genome": delta_gc,
                    "distance_to_left_contig_edge_bp": locus_start - 1,
                    "distance_to_right_contig_edge_bp": contig_len - locus_end,
                    "min_distance_to_contig_edge_bp": min_edge,
                    "near_contig_edge_2kb": int(min_edge < EDGE_WARNING_BP),
                    "mobility_gene_count_padded_region": len(mobility),
                    "mobility_gene_products_padded_region": compact_labels(mobility, limit=10),
                    "upstream_flanking_genes_nearest5": compact_labels(upstream, limit=5),
                    "downstream_flanking_genes_nearest5": compact_labels(downstream, limit=5),
                    "upstream_flanking_products_nearest5": compact_product_labels(upstream, limit=5),
                    "downstream_flanking_products_nearest5": compact_product_labels(downstream, limit=5),
                    "complete_rfbABCD_in_genome": clean(arch_row.get("complete_rfbABCD", "")),
                    "rfbABCD_colocalized_15kb": clean(arch_row.get("rfbABCD_colocalized_15kb", "")),
                    "rfb_order_by_coordinate": clean(arch_row.get("rfb_order_by_coordinate", "")),
                    "rfb_span_bp": clean(arch_row.get("rfb_span_bp", "")),
                    "wzm_wzt_colocalized_10kb": clean(arch_row.get("wzm_wzt_colocalized_10kb", "")),
                    "is_primary_rfb_region": is_primary_rfb_region,
                }
                context_rows.append(row)

                if row["source"] == "ISS" and is_primary_rfb_region:
                    iss_flank_rows.append(
                        {
                            "isolate": isolate,
                            "year": row["year"],
                            "contig": contig,
                            "target_cluster_index_on_contig": cluster_index,
                            "target_genes_in_region": genes_detected,
                            "rfb_order_by_coordinate": row["rfb_order_by_coordinate"],
                            "locus_start": locus_start,
                            "locus_end": locus_end,
                            "min_distance_to_contig_edge_bp": min_edge,
                            "delta_locus_gc_vs_genome": delta_gc,
                            "mobility_gene_count_padded_region": len(mobility),
                            "upstream_flanking_genes_nearest5": row["upstream_flanking_genes_nearest5"],
                            "downstream_flanking_genes_nearest5": row["downstream_flanking_genes_nearest5"],
                            "upstream_flanking_products_nearest5": row["upstream_flanking_products_nearest5"],
                            "downstream_flanking_products_nearest5": row["downstream_flanking_products_nearest5"],
                        }
                    )

    context = pd.DataFrame(context_rows).sort_values(["source", "species", "isolate", "contig", "locus_start"])
    mobility = pd.DataFrame(mobility_rows).sort_values(["source", "species", "isolate", "contig", "mobility_start"]) if mobility_rows else pd.DataFrame()
    iss_flanks = pd.DataFrame(iss_flank_rows).sort_values(["year", "isolate", "contig"]) if iss_flank_rows else pd.DataFrame()

    context.to_csv(OUT_CONTEXT, sep="\t", index=False)
    mobility.to_csv(OUT_MOBILITY, sep="\t", index=False)
    iss_flanks.to_csv(OUT_ISS_FLANKS, sep="\t", index=False)

    primary = context[context["is_primary_rfb_region"].eq(1)].copy()
    iss_primary = primary[primary["source"].eq("ISS")].copy()
    earth_primary = primary[primary["source"].eq("Earth")].copy()

    def fmt_pct(x):
        if pd.isna(x):
            return "NA"
        return f"{x*100:.2f}%"

    def median_text(df, col, scale=1.0, suffix=""):
        if df.empty:
            return "NA"
        return f"{df[col].median()*scale:.2f}{suffix}"

    conserved_up = iss_flanks["upstream_flanking_products_nearest5"].value_counts().head(5) if not iss_flanks.empty else pd.Series(dtype=int)
    conserved_down = iss_flanks["downstream_flanking_products_nearest5"].value_counts().head(5) if not iss_flanks.empty else pd.Series(dtype=int)

    lines = [
        "# Locus evolutionary context summary",
        "",
        f"- Regions analyzed: {len(context)} target-gene regions across {context['isolate'].nunique() if not context.empty else 0} genomes",
        f"- Primary rfb-like regions: {len(primary)}",
        f"- ISS primary rfb-like regions: {len(iss_primary)}",
        f"- Earth primary rfb-like regions: {len(earth_primary)}",
        f"- Padding used for flanking/mobility inspection: +/- {PADDING_BP:,} bp",
        "",
        "## GC context",
        "",
        f"- ISS primary median genome GC: {fmt_pct(iss_primary['genome_gc'].median()) if not iss_primary.empty else 'NA'}",
        f"- ISS primary median locus GC: {fmt_pct(iss_primary['locus_gc'].median()) if not iss_primary.empty else 'NA'}",
        f"- ISS primary median delta locus GC versus genome: {median_text(iss_primary, 'delta_locus_gc_vs_genome', scale=100, suffix=' percentage points')}",
        f"- Earth primary median delta locus GC versus genome: {median_text(earth_primary, 'delta_locus_gc_vs_genome', scale=100, suffix=' percentage points')}",
        "",
        "## Contig-edge context",
        "",
        f"- ISS primary regions near contig edge (<{EDGE_WARNING_BP:,} bp): {int(iss_primary['near_contig_edge_2kb'].sum()) if not iss_primary.empty else 0}/{len(iss_primary)}",
        f"- ISS primary median minimum distance to contig edge: {median_text(iss_primary, 'min_distance_to_contig_edge_bp', suffix=' bp')}",
        "",
        "## Mobility-gene context",
        "",
        f"- ISS primary regions with >=1 nearby mobility annotation: {int((iss_primary['mobility_gene_count_padded_region'] > 0).sum()) if not iss_primary.empty else 0}/{len(iss_primary)}",
        f"- Earth primary regions with >=1 nearby mobility annotation: {int((earth_primary['mobility_gene_count_padded_region'] > 0).sum()) if not earth_primary.empty else 0}/{len(earth_primary)}",
        f"- Total mobility-like CDS hits in padded target-gene regions: {len(mobility)}",
        "",
        "## ISS flanking conservation",
        "",
        "Most common upstream flank signatures among ISS primary regions:",
    ]
    for sig, count in conserved_up.items():
        lines.append(f"- {count} regions: {sig}")
    lines.append("")
    lines.append("Most common downstream flank signatures among ISS primary regions:")
    for sig, count in conserved_down.items():
        lines.append(f"- {count} regions: {sig}")
    lines.extend(
        [
            "",
            "## Interpretation guide",
            "",
            "Internal contig positions, conserved flanking genes, and absence of nearby mobility annotations support a stable lineage-associated locus in the available assemblies. The lower locus GC relative to genome background indicates compositional divergence of the O-antigen/rfb region, so historical recombination or ancient horizontal acquisition cannot be excluded from these data alone.",
            "",
            f"Files: `{OUT_CONTEXT.relative_to(ROOT)}`, `{OUT_ISS_FLANKS.relative_to(ROOT)}`, `{OUT_MOBILITY.relative_to(ROOT)}`",
            "",
        ]
    )
    OUT_SUMMARY.write_text("\n".join(lines))
    print(OUT_CONTEXT)
    print(OUT_ISS_FLANKS)
    print(OUT_MOBILITY)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
