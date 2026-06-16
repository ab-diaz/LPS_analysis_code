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


def find_faa(sample_dir):
    paths = glob.glob(os.path.join(sample_dir, "bakta_output", "*.faa"))
    if not paths:
        return None
    if len(paths) == 1:
        return paths[0]
    sample_name = os.path.basename(sample_dir)
    for p in paths:
        if os.path.basename(p).startswith(sample_name + "."):
            return p
    return paths[0]


def find_ffn(sample_dir):
    paths = glob.glob(os.path.join(sample_dir, "bakta_output", "*.ffn"))
    if not paths:
        return None
    if len(paths) == 1:
        return paths[0]
    sample_name = os.path.basename(sample_dir)
    for p in paths:
        if os.path.basename(p).startswith(sample_name + "."):
            return p
    return paths[0]


def load_fasta(path):
    seqs = {}
    current = None
    chunks = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                if current is not None:
                    chunks.append(line)
    if current is not None:
        seqs[current] = "".join(chunks)
    return seqs


def parse_emapper(path, strict):
    hits = []
    with open(path, "r", encoding="utf-8", newline="") as fh:
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

            if not query:
                continue

            pathway_hit = any(tag in kegg_path for tag in PATHWAY_TAGS)
            keyword_hit = KEYWORD_RE.search(" ".join([desc, pref, kegg_ko, kegg_path])) is not None
            pref_lower = pref.strip().lower()
            curated_hit = pref_lower in CURATED_GENES_LOWER

            if strict:
                keep = pathway_hit or curated_hit
            else:
                keep = pathway_hit or keyword_hit

            if not keep:
                continue

            hits.append(
                {
                    "query": query,
                    "preferred": pref,
                    "description": desc,
                    "kegg_ko": kegg_ko,
                    "kegg_pathway": kegg_path,
                    "pathway_hit": "1" if pathway_hit else "0",
                    "keyword_hit": "1" if keyword_hit else "0",
                    "curated_hit": "1" if curated_hit else "0",
                }
            )
    return hits


def main():
    parser = argparse.ArgumentParser(
        description="Extract LPS-related sequences from eggNOG annotations."
    )
    parser.add_argument(
        "--root",
        default="results",
        help="Results root directory with isolate folders.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Use curated/pathway-only filtering (stricter).",
    )
    parser.add_argument(
        "--out-proteins",
        default="results/lps_analysis/lps_eggnog_proteins.faa",
        help="Output FASTA with LPS-related protein sequences.",
    )
    parser.add_argument(
        "--out-nucleotides",
        default="results/lps_analysis/lps_eggnog_nucleotides.ffn",
        help="Output FASTA with LPS-related nucleotide sequences.",
    )
    parser.add_argument(
        "--out-hits",
        default="results/lps_analysis/lps_eggnog_hits.tsv",
        help="Output TSV with LPS-related eggNOG hits.",
    )
    args = parser.parse_args()

    sample_dirs = sorted(
        [
            os.path.join(args.root, d)
            for d in os.listdir(args.root)
            if os.path.isdir(os.path.join(args.root, d))
        ]
    )

    all_hits = []
    protein_records = []
    nucleotide_records = []

    for sdir in sample_dirs:
        isolate = os.path.basename(sdir)
        anno_path = find_emapper(sdir)
        if not anno_path:
            continue
        hits = parse_emapper(anno_path, args.strict)
        if not hits:
            continue
        query_ids = {h["query"] for h in hits}

        for h in hits:
            row = {"isolate": isolate}
            row.update(h)
            all_hits.append(row)

        faa_path = find_faa(sdir)
        if faa_path:
            faa_seqs = load_fasta(faa_path)
            for q in query_ids:
                seq = faa_seqs.get(q)
                if not seq:
                    continue
                header = f">{isolate}|{q}"
                protein_records.append((header, seq))

        ffn_path = find_ffn(sdir)
        if ffn_path:
            ffn_seqs = load_fasta(ffn_path)
            for q in query_ids:
                seq = ffn_seqs.get(q)
                if not seq:
                    continue
                header = f">{isolate}|{q}"
                nucleotide_records.append((header, seq))

    os.makedirs(os.path.dirname(args.out_hits), exist_ok=True)

    with open(args.out_hits, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
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
        for row in all_hits:
            writer.writerow(row)

    with open(args.out_proteins, "w", encoding="utf-8") as fh:
        for header, seq in protein_records:
            fh.write(header + "\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")

    with open(args.out_nucleotides, "w", encoding="utf-8") as fh:
        for header, seq in nucleotide_records:
            fh.write(header + "\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")


if __name__ == "__main__":
    main()
