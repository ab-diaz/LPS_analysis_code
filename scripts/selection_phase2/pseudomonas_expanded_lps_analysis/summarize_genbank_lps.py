#!/usr/bin/env python3
import csv
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
META = ROOT / "expanded_genbank_metadata.tsv"
HITS = ROOT / "expanded_lps_hits_genbank.tsv"
COUNTS = ROOT / "expanded_lps_category_counts.tsv"
O_MATRIX = ROOT / "expanded_o_antigen_gene_presence.tsv"

LPS_GENE_PATTERNS = [
    r"^lpx[a-z0-9]*$",
    r"^(rfa|waa|kds|kdt|hld|gmh)[a-z0-9]*$",
    r"^(rfb|rml|wbb|wzm|wzt|wzx|wzy|wbp|wzz)[a-z0-9]*$",
    r"^(ept|arn|pmr|pag|lpt|msb)[a-z0-9]*$",
    r"^(rfah|upps|diaa)$",
]

LPS_PRODUCT_TERMS = [
    "lipopolysaccharide",
    "lipid a",
    "kdo",
    "heptosyltransferase",
    "o-antigen",
    "o antigen",
    "o-polysaccharide",
    "o polysaccharide",
    "dtdp",
    "rhamnose",
    "rhamnosyl",
    "wzx",
    "wzy",
    "wzm",
    "wzt",
]

KO_CATEGORY = {
    "K00973": ("rfbA", "O-antigen pathways"),
    "K01710": ("rfbB", "O-antigen pathways"),
    "K01790": ("rfbC", "O-antigen pathways"),
    "K00067": ("rfbD", "O-antigen pathways"),
    "K09690": ("wzm", "O-antigen pathways"),
    "K09691": ("wzt", "O-antigen pathways"),
    "K02847": ("waaL", "Core oligosaccharide assembly"),
    "K02527": ("waaA", "Core oligosaccharide assembly"),
    "K02841": ("rfaC", "Core oligosaccharide assembly"),
    "K02843": ("rfaF", "Core oligosaccharide assembly"),
    "K02844": ("rfaG", "Core oligosaccharide assembly"),
    "K02848": ("rfaP", "Core oligosaccharide assembly"),
    "K03271": ("gmhA", "Core oligosaccharide assembly"),
    "K03272": ("hldE", "Core oligosaccharide assembly"),
    "K03273": ("gmhB", "Core oligosaccharide assembly"),
    "K04744": ("lptD", "LPS modification/transport"),
    "K03643": ("lptE", "LPS modification/transport"),
    "K11720": ("lptG", "LPS modification/transport"),
    "K07091": ("lptF", "LPS modification/transport"),
}


def classify_gene(gene, product="", kos=None):
    kos = kos or []
    for ko in kos:
        if ko in KO_CATEGORY:
            return KO_CATEGORY[ko][1]

    g = (gene or "").strip().lower()
    p = (product or "").strip().lower()

    if re.match(r"^lpx[a-z0-9]*$", g) or "lipid a" in p:
        return "Lipid A biosynthesis"
    if re.match(r"^(rfa|waa|kds|kdt|hld|gmh)", g):
        return "Core oligosaccharide assembly"
    if re.match(r"^(rfb|rml|wbb|wzm|wzt|wzx|wzy|wbp|wzz)", g):
        return "O-antigen pathways"
    if re.match(r"^(ept|arn|pmr|pag|lpt|msb)", g) or g in {"rfah", "upps", "diaa"}:
        return "LPS modification/transport"
    if "o-antigen" in p or "o antigen" in p or "rhamnose" in p or "dtdp" in p:
        return "O-antigen pathways"
    if "lipopolysaccharide" in p:
        return "Other LPS-related"
    return "Other LPS-related"


def preferred_gene(gene, product, kos):
    gene = (gene or "").strip()
    if gene:
        return gene
    for ko in kos:
        if ko in KO_CATEGORY:
            return KO_CATEGORY[ko][0]
    p = (product or "").lower()
    if "o-antigen ligase" in p or "o antigen ligase" in p:
        return "waaL"
    if "abc transporter permease" in p and "polysaccharide" in p:
        return "wzm"
    if "abc transporter atp" in p and "polysaccharide" in p:
        return "wzt"
    if "dtdp-glucose 4,6-dehydratase" in p:
        return "rfbB"
    if "dtdp-4-dehydrorhamnose reductase" in p:
        return "rfbD"
    if "glucose-1-phosphate thymidylyltransferase" in p:
        return "rfbA"
    return ""


def interesting(gene, product, kos):
    g = (gene or "").lower()
    p = (product or "").lower()
    if any(re.match(pattern, g) for pattern in LPS_GENE_PATTERNS):
        return True
    if any(term in p for term in LPS_PRODUCT_TERMS):
        return True
    if any(ko in KO_CATEGORY for ko in kos):
        return True
    return False


def parse_qualifier(line):
    text = line[21:].rstrip()
    if not text.startswith("/"):
        return None, None
    if "=" not in text:
        return text[1:], ""
    key, value = text[1:].split("=", 1)
    value = value.strip()
    if value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    return key, value


def parse_cds_features(path):
    rows = []
    feature = None
    current_key = None

    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("     CDS"):
                if feature:
                    rows.append(feature)
                feature = {"location": line[21:].strip(), "qualifiers": defaultdict(list)}
                current_key = None
                continue
            if not feature:
                continue
            if line.startswith("     ") and not line.startswith("                     "):
                rows.append(feature)
                feature = None
                current_key = None
                continue
            if line.startswith("                     /"):
                key, value = parse_qualifier(line)
                if key:
                    feature["qualifiers"][key].append(value)
                    current_key = key
                continue
            if line.startswith("                     ") and current_key:
                value = line[21:].strip()
                if feature["qualifiers"][current_key]:
                    feature["qualifiers"][current_key][-1] += " " + value.strip('"')
        if feature:
            rows.append(feature)
    return rows


def ko_terms(qualifiers):
    kos = []
    for xref in qualifiers.get("db_xref", []):
        for ko in re.findall(r"K\d{5}", xref):
            kos.append(ko)
    for note in qualifiers.get("note", []):
        for ko in re.findall(r"K\d{5}", note):
            kos.append(ko)
    return sorted(set(kos))


def first(qualifiers, key):
    vals = qualifiers.get(key, [])
    return vals[0] if vals else ""


def load_metadata():
    with open(META, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    metadata = load_metadata()
    all_hits = []

    for meta in metadata:
        gbff = Path(meta["gbff"])
        if not gbff.exists():
            print(f"Missing GenBank file: {gbff}")
            continue
        for feature in parse_cds_features(gbff):
            q = feature["qualifiers"]
            gene = first(q, "gene") or first(q, "locus_tag")
            product = first(q, "product")
            kos = ko_terms(q)
            pref = preferred_gene(first(q, "gene"), product, kos) or gene
            if not interesting(pref, product, kos):
                continue
            category = classify_gene(pref, product, kos)
            all_hits.append(
                {
                    "isolate": meta["isolate_id"],
                    "assembly_accession": meta["assembly_accession"],
                    "organism": meta["organism"],
                    "group": meta["group"],
                    "source_note": meta["source_note"],
                    "locus_tag": first(q, "locus_tag"),
                    "gene": pref,
                    "category": category,
                    "product": product,
                    "kegg_ko": ",".join(kos),
                    "location": feature["location"],
                }
            )

    fields = [
        "isolate",
        "assembly_accession",
        "organism",
        "group",
        "source_note",
        "locus_tag",
        "gene",
        "category",
        "product",
        "kegg_ko",
        "location",
    ]
    with open(HITS, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_hits)

    categories = [
        "Lipid A biosynthesis",
        "Core oligosaccharide assembly",
        "O-antigen pathways",
        "LPS modification/transport",
        "Other LPS-related",
    ]
    by_iso_cat = defaultdict(set)
    by_iso_gene = defaultdict(set)
    for row in all_hits:
        by_iso_cat[(row["isolate"], row["category"])].add(row["gene"])
        if row["category"] == "O-antigen pathways":
            by_iso_gene[row["isolate"]].add(row["gene"])

    with open(COUNTS, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "isolate",
                "assembly_accession",
                "organism",
                "group",
                "lipid_a_genes",
                "core_genes",
                "o_antigen_genes",
                "lps_modification_transport_genes",
                "other_lps_related_genes",
                "total_unique_lps_genes",
            ]
        )
        for meta in metadata:
            iso = meta["isolate_id"]
            vals = [len(by_iso_cat[(iso, c)]) for c in categories]
            writer.writerow(
                [
                    iso,
                    meta["assembly_accession"],
                    meta["organism"],
                    meta["group"],
                    *vals,
                    sum(vals),
                ]
            )

    genes = sorted({gene for geneset in by_iso_gene.values() for gene in geneset})
    with open(O_MATRIX, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene"] + [m["isolate_id"] for m in metadata])
        for gene in genes:
            writer.writerow([gene] + ["1" if gene in by_iso_gene[m["isolate_id"]] else "0" for m in metadata])

    print(HITS)
    print(COUNTS)
    print(O_MATRIX)


if __name__ == "__main__":
    main()
