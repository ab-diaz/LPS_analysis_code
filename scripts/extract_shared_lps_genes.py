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


def parse_emapper_strict(path):
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
            kegg_path = rec.get("KEGG_Pathway", "") or ""
            query = rec.get("query", "") or ""

            if not query or not pref or pref == "-":
                continue

            pathway_hit = any(tag in kegg_path for tag in PATHWAY_TAGS)
            curated_hit = pref.strip().lower() in CURATED_GENES_LOWER
            keep = pathway_hit or curated_hit

            if not keep:
                continue

            hits.append({
                "query": query,
                "preferred": pref.strip(),
                "description": desc,
            })
    return hits


def sanitize_gene(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def main():
    parser = argparse.ArgumentParser(
        description="Extract per-gene FASTAs for LPS genes shared across all isolates (eggNOG strict, Preferred_name)."
    )
    parser.add_argument(
        "--root",
        default="results",
        help="Results root directory with isolate folders.",
    )
    parser.add_argument(
        "--outdir",
        default="results/lps_analysis/lps_shared_genes",
        help="Output directory for per-gene FASTA files.",
    )
    args = parser.parse_args()

    sample_dirs = sorted(
        [
            os.path.join(args.root, d)
            for d in os.listdir(args.root)
            if os.path.isdir(os.path.join(args.root, d))
        ]
    )

    isolates = []
    per_isolate_hits = {}

    for sdir in sample_dirs:
        isolate = os.path.basename(sdir)
        anno_path = find_emapper(sdir)
        if not anno_path:
            continue
        hits = parse_emapper_strict(anno_path)
        if not hits:
            continue
        isolates.append(isolate)
        per_isolate_hits[isolate] = hits

    if not isolates:
        return

    # Genes present in every isolate
    gene_sets = []
    for isolate in isolates:
        gene_sets.append({h["preferred"] for h in per_isolate_hits[isolate]})
    shared_genes = set.intersection(*gene_sets) if gene_sets else set()

    os.makedirs(args.outdir, exist_ok=True)

    # Collect sequences per gene with one FASTA load per isolate
    gene_prot = {g: [] for g in shared_genes}
    gene_nucl = {g: [] for g in shared_genes}

    for isolate in isolates:
        hits = [h for h in per_isolate_hits[isolate] if h["preferred"] in shared_genes]
        if not hits:
            continue

        faa_path = find_faa(os.path.join(args.root, isolate))
        ffn_path = find_ffn(os.path.join(args.root, isolate))
        if not faa_path or not ffn_path:
            continue

        faa_seqs = load_fasta(faa_path)
        ffn_seqs = load_fasta(ffn_path)

        for h in hits:
            gene = h["preferred"]
            q = h["query"]
            prot = faa_seqs.get(q)
            nucl = ffn_seqs.get(q)
            if prot:
                header = f">{isolate}|{gene}|{q}"
                gene_prot[gene].append((header, prot))
            if nucl:
                header = f">{isolate}|{gene}|{q}"
                gene_nucl[gene].append((header, nucl))

    # Write per-gene FASTAs
    for gene in sorted(shared_genes):
        gene_safe = sanitize_gene(gene)
        out_faa = os.path.join(args.outdir, f"{gene_safe}.faa")
        out_ffn = os.path.join(args.outdir, f"{gene_safe}.ffn")

        with open(out_faa, "w", encoding="utf-8") as f_faa:
            for header, seq in gene_prot.get(gene, []):
                f_faa.write(header + "\n")
                for i in range(0, len(seq), 80):
                    f_faa.write(seq[i:i+80] + "\n")

        with open(out_ffn, "w", encoding="utf-8") as f_ffn:
            for header, seq in gene_nucl.get(gene, []):
                f_ffn.write(header + "\n")
                for i in range(0, len(seq), 80):
                    f_ffn.write(seq[i:i+80] + "\n")


if __name__ == "__main__":
    main()
