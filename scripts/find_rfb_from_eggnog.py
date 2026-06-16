#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re

TARGET = {"rfba", "rfbb", "rfbc", "rfbd"}
GENE_RE = re.compile(r"\brfb[abcd]\b", re.IGNORECASE)


def parse_gff3_attributes(attr_str):
    attrs = {}
    for part in attr_str.split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
        else:
            attrs[part] = ""
    return attrs


def load_gff3_index(gff3_path):
    idx = {}
    with open(gff3_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            contig, source, feature, start, end, score, strand, phase, attrs = parts
            if feature != "CDS":
                continue
            attr_map = parse_gff3_attributes(attrs)
            gid = attr_map.get("ID", "")
            locus = attr_map.get("locus_tag", "")
            if gid:
                idx[gid] = (contig, int(start), int(end), strand, attr_map)
            if locus:
                idx[locus] = (contig, int(start), int(end), strand, attr_map)
    return idx


def find_gff3(sample_dir):
    gff3_paths = glob.glob(os.path.join(sample_dir, "bakta_output", "*.gff3"))
    if not gff3_paths:
        return None
    if len(gff3_paths) == 1:
        return gff3_paths[0]
    sample_name = os.path.basename(sample_dir)
    for p in gff3_paths:
        if os.path.basename(p).startswith(sample_name + "."):
            return p
    return gff3_paths[0]


def find_fna(sample_dir):
    fna_paths = glob.glob(os.path.join(sample_dir, "bakta_output", "*.fna"))
    if not fna_paths:
        return None
    if len(fna_paths) == 1:
        return fna_paths[0]
    sample_name = os.path.basename(sample_dir)
    for p in fna_paths:
        if os.path.basename(p).startswith(sample_name + "."):
            return p
    return fna_paths[0]


def find_faa(sample_dir):
    faa_paths = glob.glob(os.path.join(sample_dir, "bakta_output", "*.faa"))
    if not faa_paths:
        return None
    if len(faa_paths) == 1:
        return faa_paths[0]
    sample_name = os.path.basename(sample_dir)
    for p in faa_paths:
        if os.path.basename(p).startswith(sample_name + "."):
            return p
    return faa_paths[0]


def load_fna(fna_path):
    seqs = {}
    current = None
    chunks = []
    with open(fna_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(chunks)
                header = line[1:].split()[0]
                current = header
                chunks = []
            else:
                if current is not None:
                    chunks.append(line)
    if current is not None:
        seqs[current] = "".join(chunks)
    return seqs


def load_faa(faa_path):
    seqs = {}
    current = None
    chunks = []
    with open(faa_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(chunks)
                header = line[1:].split()[0]
                current = header
                chunks = []
            else:
                if current is not None:
                    chunks.append(line)
    if current is not None:
        seqs[current] = "".join(chunks)
    return seqs


def reverse_complement(seq):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def parse_emapper(path):
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
            pref = (rec.get("Preferred_name") or "").strip()
            desc = (rec.get("Description") or "").strip()
            text = " ".join([pref, desc])
            match = GENE_RE.search(text)
            if not match:
                continue
            gene = match.group(0).lower()
            if gene not in TARGET:
                continue
            hits.append({
                "query": rec.get("query", ""),
                "gene": gene,
                "preferred": pref,
                "description": desc,
            })
    return hits


def main():
    parser = argparse.ArgumentParser(
        description="Map eggNOG rfbA-D hits to Bakta coordinates and extract operon sequences."
    )
    parser.add_argument(
        "--root",
        default="results",
        help="Results root directory with isolate folders.",
    )
    parser.add_argument(
        "--gap",
        type=int,
        default=500,
        help="Max intergenic gap (bp) to call operon-like.",
    )
    parser.add_argument(
        "--out-hits",
        default="results/lps_analysis/rfb_eggnog_locations.tsv",
        help="Output TSV with mapped locations.",
    )
    parser.add_argument(
        "--out-operons",
        default="results/lps_analysis/rfb_operon_sequences.fna",
        help="Output FASTA with operon span sequences.",
    )
    parser.add_argument(
        "--out-proteins",
        default="results/lps_analysis/rfb_proteins.faa",
        help="Output FASTA with rfbA-D protein sequences.",
    )
    args = parser.parse_args()

    sample_dirs = sorted(
        [
            os.path.join(args.root, d)
            for d in os.listdir(args.root)
            if os.path.isdir(os.path.join(args.root, d))
        ]
    )

    all_rows = []
    operon_records = []
    protein_records = []

    for sdir in sample_dirs:
        isolate = os.path.basename(sdir)
        annos = glob.glob(os.path.join(sdir, "eggnog_out", "*.emapper.annotations"))
        if not annos:
            continue
        anno_path = annos[0]
        hits = parse_emapper(anno_path)
        if not hits:
            continue

        gff3_path = find_gff3(sdir)
        if not gff3_path:
            continue
        idx = load_gff3_index(gff3_path)

        mapped = []
        for h in hits:
            query = h["query"]
            if not query:
                continue
            if query not in idx:
                continue
            contig, start, end, strand, attr_map = idx[query]
            row = {
                "isolate": isolate,
                "gene": h["gene"],
                "query": query,
                "contig": contig,
                "start": start,
                "end": end,
                "strand": strand,
                "locus_tag": attr_map.get("locus_tag", ""),
                "product": attr_map.get("product", ""),
                "preferred": h["preferred"],
                "description": h["description"],
            }
            mapped.append(row)
            all_rows.append(row)

        # Collect protein sequences for mapped hits
        faa_path = find_faa(sdir)
        if faa_path:
            faa_seqs = load_faa(faa_path)
            for row in mapped:
                seq = faa_seqs.get(row["query"])
                if not seq:
                    continue
                header = (
                    f">{isolate}|{row['gene']}|{row['query']}|"
                    f"{row['contig']}:{row['start']}-{row['end']}({row['strand']})"
                )
                protein_records.append((header, seq))

        # Operon detection for mapped hits
        if not mapped:
            continue
        groups = {}
        for row in mapped:
            key = (row["contig"], row["strand"])
            groups.setdefault(key, []).append(row)

        for (contig, strand), rows in groups.items():
            genes_found = {r["gene"] for r in rows}
            if not TARGET.issubset(genes_found):
                continue
            ordered = sorted(rows, key=lambda r: r["start"])
            gaps = []
            for a, b in zip(ordered, ordered[1:]):
                gaps.append(b["start"] - a["end"] - 1)
            if gaps and max(gaps) > args.gap:
                continue

            # Extract operon span sequence
            fna_path = find_fna(sdir)
            if not fna_path:
                continue
            seqs = load_fna(fna_path)
            if contig not in seqs:
                continue
            span_start = min(r["start"] for r in ordered)
            span_end = max(r["end"] for r in ordered)
            seq = seqs[contig][span_start - 1:span_end]
            if strand == "-":
                seq = reverse_complement(seq)
            header = f">{isolate}|{contig}|{strand}|{span_start}-{span_end}|rfbA-rfbD"
            operon_records.append((header, seq))

    os.makedirs(os.path.dirname(args.out_hits), exist_ok=True)

    with open(args.out_hits, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "isolate",
                "gene",
                "query",
                "contig",
                "start",
                "end",
                "strand",
                "locus_tag",
                "product",
                "preferred",
                "description",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    with open(args.out_operons, "w", encoding="utf-8") as fh:
        for header, seq in operon_records:
            fh.write(header + "\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")

    with open(args.out_proteins, "w", encoding="utf-8") as fh:
        for header, seq in protein_records:
            fh.write(header + "\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")


if __name__ == "__main__":
    main()
