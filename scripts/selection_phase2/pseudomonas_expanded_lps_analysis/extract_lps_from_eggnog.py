#!/usr/bin/env python3
import csv
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
META = ROOT / "selected_genomes.tsv"
HITS = ROOT / "expanded_lps_eggnog_hits.tsv"
CATEGORIZED = ROOT / "expanded_lps_hits_categorized.tsv"
COUNTS = ROOT / "expanded_lps_category_counts.tsv"
O_MATRIX = ROOT / "expanded_o_antigen_gene_presence.tsv"

CURATED_LPS_GENES = {
    "diaa",
    "gmha",
    "gmhb",
    "hlde",
    "kdsa",
    "kdsb",
    "kdsc",
    "kdsd",
    "kdta",
    "lpta",
    "lptb",
    "lptc",
    "lptd",
    "lpte",
    "lptf",
    "lptg",
    "lpxa",
    "lpxb",
    "lpxc",
    "lpxd",
    "lpxh",
    "lpxk",
    "lpxl",
    "lpxm",
    "msba",
    "rfah",
    "rfba",
    "rfbb",
    "rfbc",
    "rfbd",
    "rfbf",
    "rfbg",
    "rmlA".lower(),
    "rmlB".lower(),
    "rmlC".lower(),
    "rmlD".lower(),
    "rfac",
    "rfaf",
    "rfag",
    "rfap",
    "rfaz",
    "waaa",
    "waac",
    "waaf",
    "waag",
    "waal",
    "waap",
    "wbb",
    "wbp",
    "wzm",
    "wzt",
    "wzx",
    "wzy",
    "wzz",
}

LPS_KEYWORDS = [
    "lipopolysaccharide",
    "lipid a",
    "o-antigen",
    "o antigen",
    "o-polysaccharide",
    "heptosyltransferase",
    "kdo",
    "dtdp",
    "rhamnose",
    "rhamnosyl",
]

LPS_KOS = {
    "ko:K00973",
    "ko:K01710",
    "ko:K01790",
    "ko:K00067",
    "ko:K09690",
    "ko:K09691",
    "ko:K02847",
    "ko:K02527",
    "ko:K02841",
    "ko:K02843",
    "ko:K02844",
    "ko:K02848",
    "ko:K03271",
    "ko:K03272",
    "ko:K03273",
    "ko:K04744",
    "ko:K03643",
    "ko:K11720",
    "ko:K07091",
}


def classify_gene(gene):
    g = (gene or "").strip().lower()

    if re.match(r"^lpx[a-z0-9]*$", g):
        return "Lipid A biosynthesis" if g != "lpxt" else "LPS modification/transport"
    if re.match(r"^(rfa|waa|kds|kdt|hld|gmh)", g):
        return "Core oligosaccharide assembly"
    if re.match(r"^(rfb|rml|wbb|wzm|wzt|wzx|wzy|wbp|wzz)", g):
        return "O-antigen pathways"
    if re.match(r"^(ept|arn|pmr|pag|lpt|msb)", g):
        return "LPS modification/transport"
    if g in {"rfah", "upps", "diaa"}:
        return "LPS modification/transport"
    return "Other LPS-related"


def norm_gene(name):
    return (name or "").strip().split(",")[0].lower()


def load_metadata():
    with open(META, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_emapper(path, isolate):
    rows = []
    header = None
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").split("\t")
                continue
            if line.startswith("#"):
                continue
            if header is None:
                continue
            values = line.split("\t")
            if len(values) < len(header):
                values += [""] * (len(header) - len(values))
            row = dict(zip(header, values))
            preferred = row.get("Preferred_name", "") or row.get("query", "")
            gene = norm_gene(preferred)
            desc = row.get("Description", "")
            kos = {x for x in re.split(r"[,;]", row.get("KEGG_ko", "")) if x and x != "-"}
            pathways = row.get("KEGG_Pathway", "")
            pathway_hit = any(x in pathways for x in ["ko00540", "map00540"])
            keyword_hit = any(k in desc.lower() for k in LPS_KEYWORDS)
            curated_hit = gene in CURATED_LPS_GENES or bool(kos & LPS_KOS)
            if not (pathway_hit or keyword_hit or curated_hit):
                continue
            rows.append(
                {
                    "isolate": isolate,
                    "query": row.get("query", ""),
                    "preferred": gene,
                    "description": desc,
                    "kegg_ko": row.get("KEGG_ko", ""),
                    "kegg_pathway": pathways,
                    "pathway_hit": int(pathway_hit),
                    "keyword_hit": int(keyword_hit),
                    "curated_hit": int(curated_hit),
                }
            )
    return rows


def main():
    metadata = load_metadata()
    by_iso = {m["isolate_id"]: m for m in metadata}

    hits = []
    for meta in metadata:
        isolate = meta["isolate_id"]
        annot = ROOT / "results" / isolate / "eggnog_out" / f"{isolate}.emapper.annotations"
        if not annot.exists():
            print(f"Missing eggNOG annotation, skipping: {annot}")
            continue
        hits.extend(parse_emapper(annot, isolate))

    with open(HITS, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "isolate",
                "query",
                "preferred",
                "description",
                "kegg_ko",
                "kegg_pathway",
                "pathway_hit",
                "keyword_hit",
                "curated_hit",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(hits)

    categorized = []
    for row in hits:
        meta = by_iso.get(row["isolate"], {})
        gene = row["preferred"]
        categorized.append(
            {
                "isolate": row["isolate"],
                "accession": meta.get("accession", ""),
                "group": meta.get("group", ""),
                "location": meta.get("location", ""),
                "source_note": meta.get("source_note", ""),
                "query": row["query"],
                "gene": gene,
                "category": classify_gene(gene),
                "description": row["description"],
                "kegg_ko": row["kegg_ko"],
                "kegg_pathway": row["kegg_pathway"],
            }
        )

    with open(CATEGORIZED, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "isolate",
                "accession",
                "group",
                "location",
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
        writer.writerows(categorized)

    categories = [
        "Lipid A biosynthesis",
        "Core oligosaccharide assembly",
        "O-antigen pathways",
        "LPS modification/transport",
        "Other LPS-related",
    ]
    genes_by_iso_cat = defaultdict(set)
    for row in categorized:
        genes_by_iso_cat[(row["isolate"], row["category"])].add(row["gene"])

    with open(COUNTS, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "isolate",
                "accession",
                "group",
                "location",
                "lipid_a_genes",
                "core_genes",
                "o_antigen_genes",
                "lps_modification_transport_genes",
                "other_lps_related_genes",
                "total_unique_lps_genes",
            ]
        )
        for meta in metadata:
            isolate = meta["isolate_id"]
            vals = [len(genes_by_iso_cat[(isolate, c)]) for c in categories]
            writer.writerow(
                [
                    isolate,
                    meta["accession"],
                    meta["group"],
                    meta["location"],
                    *vals,
                    sum(vals),
                ]
            )

    o_rows = [r for r in categorized if r["category"] == "O-antigen pathways"]
    o_genes = sorted({r["gene"] for r in o_rows})
    with open(O_MATRIX, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene"] + [m["isolate_id"] for m in metadata])
        for gene in o_genes:
            writer.writerow(
                [gene]
                + [
                    "1"
                    if any(r["isolate"] == m["isolate_id"] and r["gene"] == gene for r in o_rows)
                    else "0"
                    for m in metadata
                ]
            )

    for path in [HITS, CATEGORIZED, COUNTS, O_MATRIX]:
        print(path)


if __name__ == "__main__":
    main()
