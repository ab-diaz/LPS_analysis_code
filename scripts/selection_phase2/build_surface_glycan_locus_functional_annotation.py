#!/usr/bin/env python3
from pathlib import Path
import re

import pandas as pd


ROOT = Path(__file__).resolve().parent
OUTDIR = ROOT / "pantoea_route1_strengthening"
GENOME_INDEX = ROOT / "phylogeny_aware_inputs" / "genome_file_index.tsv"
CONTEXT = OUTDIR / "locus_evolutionary_context_by_region.tsv"
VALIDATION = OUTDIR / "oantigen_locus_validation_table.tsv"
OUT_TSV = OUTDIR / "surface_glycan_locus_functional_annotation.tsv"
OUT_SUMMARY = OUTDIR / "surface_glycan_locus_functional_annotation_summary.md"

TARGET_MAP = {
    "rfba": "rfbA",
    "rfbb": "rfbB",
    "rfbc": "rfbC",
    "rfbd": "rfbD",
    "waal": "waaL",
    "wzm": "wzm",
    "wzt": "wzt",
    "rfaz": "rfaZ",
}

MOBILITY_RE = re.compile(
    r"transposase|integrase|recombinase|insertion sequence|\bIS\b|phage|prophage|"
    r"plasmid|conjugative|relaxase|mobilization|mobilisation|resolvase|invertase",
    re.I,
)


def clean(value):
    if pd.isna(value):
        return ""
    value = str(value)
    return "" if value.lower() == "nan" or value == "-" else value


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
            if line.startswith("#") or header is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            rows.append(dict(zip(header, parts)))
    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame()
    df = df.rename(
        columns={
            "Sequence Id": "contig",
            "Type": "type",
            "Start": "start",
            "Stop": "end",
            "Strand": "strand",
            "Locus Tag": "locus_tag",
            "Gene": "bakta_gene",
            "Product": "bakta_product",
            "DbXrefs": "bakta_dbxrefs",
        }
    )
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    return df


def read_emapper(path):
    if not Path(path).exists():
        return {}
    header = None
    rows = {}
    with Path(path).open(errors="replace") as handle:
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
            rows[rec.get("query", "")] = rec
    return rows


def classify_function(gene, product, preferred, description, kegg_ko, kegg_pathway):
    text = " ".join([gene, product, preferred, description, kegg_ko, kegg_pathway]).lower()
    gene_l = gene.lower()

    if gene_l in {"rfba", "rfbb", "rfbc", "rfbd"}:
        return "dTDP-L-rhamnose/nucleotide-sugar biosynthesis"
    if gene_l == "waal":
        return "O-antigen ligation"
    if gene_l in {"wzm", "wzt"} or any(ko in text for ko in ["ko:k09688", "ko:k09690", "ko:k09691", "ko:k09692"]):
        return "O-antigen ABC transport"
    if gene_l.startswith("rfa") or gene_l.startswith("waa") or "lipopolysaccharide" in text:
        return "LPS core/outer-membrane assembly"
    if MOBILITY_RE.search(text):
        return "mobile element/recombination"
    if any(x in text for x in ["glycosyltransferase", "glycosyl transferase", "glycosyl-transferase"]):
        return "glycosyltransferase/surface-glycan assembly"
    if any(
        x in text
        for x in [
            "udp-glucose",
            "utp--glucose",
            "glucose-1-phosphate",
            "thymidylyltransferase",
            "epimerase",
            "dehydratase",
            "reductase",
            "rhamnose",
            "nucleotide sugar",
            "nucleotide-sugar",
        ]
    ):
        return "nucleotide-sugar/polysaccharide precursor metabolism"
    if any(
        x in text
        for x in [
            "amylovoran",
            "colanic acid",
            "polysaccharide",
            "exopolysaccharide",
            "capsule",
            "hyaluronan",
            "cell wall biosynthesis",
        ]
    ):
        return "surface polysaccharide/exopolysaccharide biosynthesis"
    if any(x in text for x in ["abc transporter", "transport permease", "membrane", "permease", "transporter"]):
        return "membrane transport"
    if "hypothetical" in text:
        return "hypothetical/unknown"
    return "other flanking function"


def evidence_tags(category, gene, product, preferred, description, kegg_ko, kegg_pathway):
    tags = []
    text = " ".join([gene, product, preferred, description, kegg_ko, kegg_pathway]).lower()
    if gene.lower() in TARGET_MAP:
        tags.append("focal O-antigen/LPS gene")
    if "ko00521" in kegg_pathway or "ko00523" in kegg_pathway or "ko00525" in kegg_pathway:
        tags.append("KEGG sugar/O-antigen-related pathway")
    if "ko00540" in kegg_pathway:
        tags.append("KEGG LPS pathway")
    if "ko02010" in kegg_pathway:
        tags.append("KEGG ABC transporter pathway")
    if any(x in text for x in ["glycosyl", "polysaccharide", "amylovoran", "colanic", "lipopolysaccharide"]):
        tags.append("surface-glycan product keyword")
    if "hypothetical" in text:
        tags.append("hypothetical product")
    if category == "mobile element/recombination":
        tags.append("mobility keyword")
    return ";".join(dict.fromkeys(tags))


def main():
    index = pd.read_csv(GENOME_INDEX, sep="\t").set_index("genome")
    context = pd.read_csv(CONTEXT, sep="\t")
    validation = pd.read_csv(VALIDATION, sep="\t")
    validation_lookup = {
        (r.isolate, r.bakta_locus_tag): r.gene for r in validation.itertuples(index=False)
    }

    primary = context[context["is_primary_rfb_region"].eq(1)].copy()
    rows = []
    for region in primary.itertuples(index=False):
        isolate = region.isolate
        if isolate not in index.index:
            continue
        meta = index.loc[isolate]
        bakta_tsv = clean(meta["bakta_tsv"])
        eggnog = clean(meta["eggnog_annotations"])
        if not bakta_tsv:
            continue
        features = read_bakta_tsv(bakta_tsv)
        if features.empty:
            continue
        emapper = read_emapper(eggnog)
        cds = features[
            features["type"].str.lower().eq("cds")
            & features["contig"].eq(region.contig)
            & (features["start"].astype(int) <= int(region.padded_region_end))
            & (features["end"].astype(int) >= int(region.padded_region_start))
        ].sort_values(["start", "end"])

        for feature in cds.itertuples(index=False):
            locus = clean(feature.locus_tag)
            ann = emapper.get(locus, {})
            gene = clean(feature.bakta_gene)
            product = clean(feature.bakta_product)
            preferred = clean(ann.get("Preferred_name", ""))
            description = clean(ann.get("Description", ""))
            kegg_ko = clean(ann.get("KEGG_ko", ""))
            kegg_pathway = clean(ann.get("KEGG_Pathway", ""))
            standard_target = validation_lookup.get((isolate, locus), "")
            category = classify_function(gene or preferred or standard_target, product, preferred, description, kegg_ko, kegg_pathway)
            rows.append(
                {
                    "isolate": isolate,
                    "source": clean(region.source),
                    "year": clean(region.year),
                    "species": clean(region.species),
                    "contig": clean(region.contig),
                    "target_cluster_index_on_contig": int(region.target_cluster_index_on_contig),
                    "region_start": int(region.padded_region_start),
                    "region_end": int(region.padded_region_end),
                    "region_target_genes": clean(region.target_genes_in_region),
                    "locus_tag": locus,
                    "start": int(feature.start),
                    "end": int(feature.end),
                    "strand": clean(feature.strand),
                    "bakta_gene": gene,
                    "standard_focal_gene": standard_target,
                    "bakta_product": product,
                    "eggnog_preferred_name": preferred,
                    "eggnog_description": description,
                    "kegg_ko": kegg_ko,
                    "kegg_pathway": kegg_pathway,
                    "functional_category": category,
                    "surface_glycan_support": int(
                        category
                        in {
                            "dTDP-L-rhamnose/nucleotide-sugar biosynthesis",
                            "O-antigen ligation",
                            "O-antigen ABC transport",
                            "LPS core/outer-membrane assembly",
                            "glycosyltransferase/surface-glycan assembly",
                            "nucleotide-sugar/polysaccharide precursor metabolism",
                            "surface polysaccharide/exopolysaccharide biosynthesis",
                            "membrane transport",
                        }
                    ),
                    "evidence_tags": evidence_tags(category, gene or preferred or standard_target, product, preferred, description, kegg_ko, kegg_pathway),
                }
            )

    out = pd.DataFrame(rows)
    out.to_csv(OUT_TSV, sep="\t", index=False)

    iss = out[out["source"].eq("ISS")]
    cat_counts = (
        iss.groupby("functional_category")
        .size()
        .sort_values(ascending=False)
        .reset_index(name="cds_count_in_iss_primary_regions")
    )
    support_by_region = (
        iss.groupby(["isolate", "year", "contig", "target_cluster_index_on_contig"])
        .agg(
            cds_total=("locus_tag", "count"),
            surface_glycan_supporting_cds=("surface_glycan_support", "sum"),
            categories=("functional_category", lambda x: ";".join(sorted(set(x)))),
        )
        .reset_index()
    )

    lines = [
        "# Surface-glycan locus functional annotation",
        "",
        f"- Regions annotated: {primary.shape[0]} primary rfb-like regions",
        f"- CDS annotated in primary regions: {out.shape[0]}",
        f"- ISS CDS annotated in primary regions: {iss.shape[0]}",
        f"- Output table: `{OUT_TSV.relative_to(ROOT)}`",
        "",
        "This table classifies CDSs within the +/-10 kb neighborhoods of primary rfb-like regions using Bakta products, eggNOG preferred names/descriptions, and KEGG annotations.",
        "",
        "## ISS functional-category counts",
        "",
        "| functional category | CDS count |",
        "| --- | ---: |",
    ]
    for r in cat_counts.itertuples(index=False):
        lines.append(f"| {r.functional_category} | {int(r.cds_count_in_iss_primary_regions)} |")

    lines.extend(
        [
            "",
            "## Region-level support in ISS",
            "",
            f"- ISS primary regions: {support_by_region.shape[0]}",
            f"- Median CDS per ISS region: {support_by_region['cds_total'].median():.1f}",
            f"- Median surface-glycan-supporting CDS per ISS region: {support_by_region['surface_glycan_supporting_cds'].median():.1f}",
            "",
            "Interpretation: concentration of rfb genes, nucleotide-sugar enzymes, glycosyltransferases, O-antigen/ABC transport components, LPS-associated products, and polysaccharide-biosynthesis genes supports interpretation of the conserved ISS-associated region as a surface-glycan/O-antigen-associated locus.",
            "",
        ]
    )
    OUT_SUMMARY.write_text("\n".join(lines))
    print(OUT_TSV)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
