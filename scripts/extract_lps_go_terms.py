#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re

PATHWAY_TAGS = {"map00540", "ko00540"}
KEYWORD_RE = re.compile(
    r"\b(lipopolysaccharide|lipid\s*A|O-antigen|LPS|"
    r"lpx\w*|lpt\w*|waa\w*|rfa\w*|rfb\w*|"
    r"wzm\w*|wzt\w*|kds\w*|kdtA)\b",
    re.IGNORECASE,
)
CURATED_GENES = {
    "lpxa", "lpxb", "lpxc", "lpxd", "lpxe", "lpxh", "lpxl", "lpxm",
    "lpta", "lptb", "lptc", "lptd", "lpte", "lptf", "lptg",
    "waaa", "waac", "waaf", "waag", "waal", "waao", "waap",
    "rfae", "rfaa", "rfab", "rfac", "rfad", "rfaz",
    "rfba", "rfbb", "rfbc", "rfbd", "rfbf", "rfbg", "rfbh", "rfbi", "rfbj",
    "wzm", "wzt", "wbbm", "wbbn", "wbbp", "wbbq",
    "kdsa", "kdsb", "kdsc", "kdta",
    "hldd", "hlde", "gmma", "gmha", "gmhb", "rfaq",
    "epta", "eptb", "diaa",
}
CURATED_GENES_LOWER = {g.lower() for g in CURATED_GENES}

GO_SPLIT_RE = re.compile(r"[;,|]")
GO_ID_RE = re.compile(r"GO:\d{7}")


def find_emapper(sample_dir):
    paths = glob.glob(os.path.join(sample_dir, "eggnog_out", "*.emapper.annotations"))
    if not paths:
        return None
    if len(paths) == 1:
        return paths[0]
    sample_name = os.path.basename(sample_dir)
    for p in paths:
        if os.path.basename(p).startswith(sample_name + "."):
            return p
    return paths[0]


def load_go_obo(obo_path):
    go = {}
    if not obo_path or not os.path.exists(obo_path):
        return go
    current = None
    in_term = False
    with open(obo_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                if current and "id" in current:
                    go[current["id"]] = current
                current = {}
                in_term = True
                continue
            if line.startswith("[") and line.endswith("]"):
                if current and "id" in current:
                    go[current["id"]] = current
                current = None
                in_term = False
                continue
            if not in_term or current is None:
                continue
            if line.startswith("id: "):
                current["id"] = line.split("id: ", 1)[1].strip()
            elif line.startswith("name: "):
                current["name"] = line.split("name: ", 1)[1].strip()
            elif line.startswith("namespace: "):
                current["namespace"] = line.split("namespace: ", 1)[1].strip()
            elif line.startswith("def: "):
                definition = line.split("def: ", 1)[1].strip()
                if definition.startswith("\""):
                    definition = definition.split("\"", 2)[1]
                current["def"] = definition
        if current and "id" in current:
            go[current["id"]] = current
    return go


def find_go_column(header):
    candidates = ["GOs", "GO_terms", "GO", "GOterm", "GO_terms"]
    for c in candidates:
        if c in header:
            return header.index(c)
    for i, col in enumerate(header):
        if col.lower().startswith("go"):
            return i
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Build GO-term table for eggNOG LPS-related genes."
    )
    parser.add_argument(
        "--root",
        default="results",
        help="Results root directory with isolate folders.",
    )
    parser.add_argument(
        "--obo",
        default="results/go-basic.obo",
        help="GO OBO file for term names/definitions.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Use curated/pathway-only filtering (stricter).",
    )
    parser.add_argument(
        "--out",
        default="results/lps_analysis/lps_eggnog_go_terms.tsv",
        help="Output TSV file.",
    )
    args = parser.parse_args()

    go_info = load_go_obo(args.obo)

    sample_dirs = sorted(
        [
            os.path.join(args.root, d)
            for d in os.listdir(args.root)
            if os.path.isdir(os.path.join(args.root, d))
        ]
    )

    term_hits = {}

    for sdir in sample_dirs:
        isolate = os.path.basename(sdir)
        anno_path = find_emapper(sdir)
        if not anno_path:
            continue
        with open(anno_path, "r", encoding="utf-8", newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = None
            for row in reader:
                if not row:
                    continue
                if row[0].startswith("##"):
                    continue
                if row[0].startswith("#"):
                    header = [c.lstrip("#") for c in row]
                    continue
                if header is None:
                    continue

                rec = dict(zip(header, row))
                desc = rec.get("Description", "") or ""
                pref = rec.get("Preferred_name", "") or ""
                kegg_ko = rec.get("KEGG_ko", "") or ""
                kegg_path = rec.get("KEGG_Pathway", "") or ""
                query = rec.get("query", "") or ""

                pathway_hit = any(tag in kegg_path for tag in PATHWAY_TAGS)
                keyword_hit = KEYWORD_RE.search(" ".join([desc, pref, kegg_ko, kegg_path])) is not None
                curated_hit = pref.strip().lower() in CURATED_GENES_LOWER

                if args.strict:
                    keep = pathway_hit or curated_hit
                else:
                    keep = pathway_hit or keyword_hit

                if not keep:
                    continue

                go_col_idx = find_go_column(header)
                if go_col_idx is None or go_col_idx >= len(row):
                    continue
                go_field = row[go_col_idx].strip()
                if not go_field or go_field == "-":
                    continue

                terms = []
                for part in GO_SPLIT_RE.split(go_field):
                    part = part.strip()
                    if not part:
                        continue
                    match = GO_ID_RE.search(part)
                    if match:
                        terms.append(match.group(0))

                if not terms:
                    continue

                for go_id in terms:
                    entry = term_hits.setdefault(go_id, {"isolates": set(), "genes": set(), "hits": 0})
                    entry["isolates"].add(isolate)
                    if query:
                        entry["genes"].add(query)
                    entry["hits"] += 1

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "go_id",
            "go_name",
            "go_namespace",
            "go_definition",
            "lps_gene_count",
            "isolates_count",
            "hit_count",
        ])
        for go_id in sorted(term_hits.keys()):
            info = go_info.get(go_id, {})
            writer.writerow([
                go_id,
                info.get("name", ""),
                info.get("namespace", ""),
                info.get("def", ""),
                str(len(term_hits[go_id]["genes"])),
                str(len(term_hits[go_id]["isolates"])),
                str(term_hits[go_id]["hits"]),
            ])


if __name__ == "__main__":
    main()
