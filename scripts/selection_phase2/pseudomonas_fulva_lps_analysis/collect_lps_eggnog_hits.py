#!/usr/bin/env python3
"""Collect LPS/O-antigen related eggNOG hits for the P. fulva side analysis."""

from __future__ import annotations

import csv
import os
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
META = Path(os.environ.get("META_FILE", ROOT / "genome_metadata.tsv"))
OUT = Path(os.environ.get("HITS_FILE", ROOT / "lps_eggnog_hits.tsv"))

GENE_RE = re.compile(
    r"^(lpx|kds|kdt|gmh|hld|waa|rfa|rfb|rml|wbb|wbp|wzx|wzy|wzz|wzm|wzt|lpt|msb|ept|arn|pmr|pag)",
    re.I,
)
KEYWORD_RE = re.compile(
    r"lipopolysaccharide|lipid a|o-antigen|o antigen|rhamnose|dtdp|glycosyltransferase|polysaccharide|outer membrane",
    re.I,
)
PATHWAY_RE = re.compile(r"ko00540|lipopolysaccharide", re.I)


def read_metadata() -> list[str]:
    with META.open(newline="", encoding="utf-8") as handle:
        return [r["isolate_id"] for r in csv.DictReader(handle, delimiter="\t")]


def parse_annotations(path: Path):
    header = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                continue
            if not line.strip() or header is None:
                continue
            vals = line.rstrip("\n").split("\t")
            if len(vals) < len(header):
                vals += [""] * (len(header) - len(vals))
            yield dict(zip(header, vals))


def clean(value: str) -> str:
    return "" if value in {"-", "NA", "None"} else value


def main() -> None:
    isolates = read_metadata()
    rows = []
    for isolate in isolates:
        annot = ROOT / "results" / isolate / "eggnog_out" / f"{isolate}.emapper.annotations"
        if not annot.exists():
            print(f"Missing eggNOG annotations, skipping: {isolate}")
            continue
        for row in parse_annotations(annot):
            preferred = clean(row.get("Preferred_name", ""))
            desc = clean(row.get("Description", ""))
            ko = clean(row.get("KEGG_ko", ""))
            pathway = clean(row.get("KEGG_Pathway", ""))
            text = " ".join([preferred, desc, ko, pathway])
            pathway_hit = bool(PATHWAY_RE.search(pathway))
            keyword_hit = bool(KEYWORD_RE.search(text))
            curated_hit = bool(GENE_RE.search(preferred))
            if not (pathway_hit or keyword_hit or curated_hit):
                continue
            rows.append(
                {
                    "isolate": isolate,
                    "query": row.get("query", row.get("#query", "")),
                    "preferred": preferred,
                    "description": desc,
                    "kegg_ko": ko,
                    "kegg_pathway": pathway,
                    "pathway_hit": int(pathway_hit),
                    "keyword_hit": int(keyword_hit),
                    "curated_hit": int(curated_hit),
                }
            )

    with OUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
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
        )
        writer.writeheader()
        writer.writerows(rows)
    print(OUT)
    print(f"LPS/O-antigen related hits: {len(rows)}")


if __name__ == "__main__":
    main()
