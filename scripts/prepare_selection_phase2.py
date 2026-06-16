#!/usr/bin/env python3
import argparse
import csv
import os
import re
import shutil
import subprocess
from collections import defaultdict


TARGETS = ("rfbA", "rfbB", "rfbC", "rfbD", "waaL")


def read_fasta(path):
    records = {}
    header = None
    seq_parts = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            records[header] = "".join(seq_parts)
    return records


def write_fasta(records, path, width=80):
    with open(path, "w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i : i + width] + "\n")


def parse_header(header):
    parts = header.split("|")
    if len(parts) < 4:
        raise ValueError(f"Unexpected curated FASTA header: {header}")
    return {
        "isolate": parts[0],
        "group": parts[1],
        "gene": parts[2],
        "locus": parts[3],
    }


def safe_id(header):
    meta = parse_header(header)
    raw = f"{meta['gene']}_{meta['group']}_{meta['isolate']}_{meta['locus']}"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", raw)


def find_diamond(explicit=None):
    candidates = [
        explicit,
        shutil.which("diamond"),
        "diamond",
        "diamond",
        "diamond",
    ]
    for path in candidates:
        if path and os.path.exists(path):
            return path
    return None


def load_curated(indir):
    nuc = {}
    prot = {}
    for gene in TARGETS:
        nuc[gene] = read_fasta(os.path.join(indir, f"{gene}.ffn"))
        prot[gene] = read_fasta(os.path.join(indir, f"{gene}.faa"))
    return nuc, prot


def choose_reference(prot_records, requested):
    complete = defaultdict(set)
    headers_by_gene_iso = {}
    for gene, records in prot_records.items():
        for header in records:
            meta = parse_header(header)
            complete[meta["isolate"]].add(gene)
            headers_by_gene_iso[(gene, meta["isolate"])] = header

    if requested:
        if complete.get(requested, set()) >= set(TARGETS):
            ref = requested
        else:
            raise ValueError(f"Requested reference isolate does not have all target genes: {requested}")
    else:
        valid = sorted(iso for iso, genes in complete.items() if genes >= set(TARGETS))
        if not valid:
            raise ValueError("No isolate has all target genes after curation.")
        ref = valid[0]

    reference = []
    for gene in TARGETS:
        header = headers_by_gene_iso[(gene, ref)]
        reference.append((f"REF|{gene}|{ref}|{parse_header(header)['locus']}", prot_records[gene][header]))
    return ref, reference


def run_diamond(diamond, reference_faa, query_faa, out_tsv, outdir):
    db_prefix = os.path.join(outdir, "reference_targets.dmnd")
    subprocess.run([diamond, "makedb", "--in", reference_faa, "-d", db_prefix], check=True)
    subprocess.run(
        [
            diamond,
            "blastp",
            "-d",
            db_prefix,
            "-q",
            query_faa,
            "-o",
            out_tsv,
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "qlen",
            "slen",
            "evalue",
            "bitscore",
            "--max-target-seqs",
            "1",
            "--quiet",
        ],
        check=True,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Prepare curated rfbABCD/waaL FASTAs for codon alignment and HyPhy RELAX."
    )
    parser.add_argument("--indir", default="lps_analysis/selection_curation")
    parser.add_argument("--outdir", default="lps_analysis/selection_phase2")
    parser.add_argument("--reference-isolate", default="Pantoea_piersonii_F10_5S-D1_P5")
    parser.add_argument("--diamond", default=None)
    parser.add_argument("--min-pident", type=float, default=40.0)
    parser.add_argument("--min-qcov", type=float, default=70.0)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    nuc_records, prot_records = load_curated(args.indir)
    ref_isolate, reference = choose_reference(prot_records, args.reference_isolate)

    reference_faa = os.path.join(args.outdir, "reference_targets.faa")
    query_faa = os.path.join(args.outdir, "curated_all_targets.faa")
    diamond_tsv = os.path.join(args.outdir, "diamond_best_hits.tsv")

    write_fasta(reference, reference_faa)
    all_queries = []
    for gene in TARGETS:
        for header, seq in prot_records[gene].items():
            all_queries.append((header, seq))
    write_fasta(all_queries, query_faa)

    diamond = find_diamond(args.diamond)
    if diamond:
        run_diamond(diamond, reference_faa, query_faa, diamond_tsv, args.outdir)
    else:
        raise RuntimeError("DIAMOND was not found. Use --diamond /path/to/diamond.")

    best = {}
    with open(diamond_tsv) as handle:
        for line in handle:
            qseqid, sseqid, pident, length, qlen, slen, evalue, bitscore = line.rstrip("\n").split("\t")
            subject_gene = sseqid.split("|")[1]
            qcov = 100.0 * float(length) / float(qlen)
            best[qseqid] = {
                "subject": sseqid,
                "subject_gene": subject_gene,
                "pident": float(pident),
                "qcov": qcov,
                "evalue": evalue,
                "bitscore": float(bitscore),
            }

    decisions = []
    kept_headers = defaultdict(set)
    for gene in TARGETS:
        for header, prot_seq in prot_records[gene].items():
            hit = best.get(header)
            keep = False
            reason = "missing_diamond_hit"
            if hit:
                if hit["subject_gene"] != gene:
                    reason = "best_hit_wrong_gene"
                elif hit["pident"] < args.min_pident:
                    reason = "low_identity"
                elif hit["qcov"] < args.min_qcov:
                    reason = "low_query_coverage"
                else:
                    keep = True
                    reason = "diamond_best_hit_pass"
                    kept_headers[gene].add(header)
            meta = parse_header(header)
            decisions.append(
                {
                    "query": header,
                    "isolate": meta["isolate"],
                    "group": meta["group"],
                    "gene": gene,
                    "keep": int(keep),
                    "reason": reason,
                    "subject": hit["subject"] if hit else "",
                    "subject_gene": hit["subject_gene"] if hit else "",
                    "pident": f"{hit['pident']:.3f}" if hit else "",
                    "qcov": f"{hit['qcov']:.3f}" if hit else "",
                    "evalue": hit["evalue"] if hit else "",
                    "bitscore": f"{hit['bitscore']:.1f}" if hit else "",
                }
            )

    for gene in TARGETS:
        checked_nuc = [(h, s) for h, s in nuc_records[gene].items() if h in kept_headers[gene]]
        checked_prot = [(h, s) for h, s in prot_records[gene].items() if h in kept_headers[gene]]
        write_fasta(checked_nuc, os.path.join(args.outdir, f"{gene}.diamond_checked.ffn"))
        write_fasta(checked_prot, os.path.join(args.outdir, f"{gene}.diamond_checked.faa"))
        write_fasta(
            [(safe_id(h), s) for h, s in checked_nuc],
            os.path.join(args.outdir, f"{gene}.relax_input.ffn"),
        )
        write_fasta(
            [(safe_id(h), s) for h, s in checked_prot],
            os.path.join(args.outdir, f"{gene}.relax_input.faa"),
        )

    with open(os.path.join(args.outdir, "diamond_curation_decisions.tsv"), "w", newline="") as handle:
        fieldnames = [
            "query",
            "isolate",
            "group",
            "gene",
            "keep",
            "reason",
            "subject",
            "subject_gene",
            "pident",
            "qcov",
            "evalue",
            "bitscore",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(decisions)

    with open(os.path.join(args.outdir, "sequence_id_map.tsv"), "w", newline="") as handle:
        fieldnames = ["safe_id", "original_id", "gene", "group", "isolate", "locus"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for gene in TARGETS:
            for header in sorted(kept_headers[gene]):
                meta = parse_header(header)
                writer.writerow(
                    {
                        "safe_id": safe_id(header),
                        "original_id": header,
                        "gene": gene,
                        "group": meta["group"],
                        "isolate": meta["isolate"],
                        "locus": meta["locus"],
                    }
                )

    with open(os.path.join(args.outdir, "diamond_curation_summary.tsv"), "w", newline="") as handle:
        fieldnames = ["gene", "kept_total", "kept_ISS", "kept_Earth", "rejected_total"]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for gene in TARGETS:
            gene_rows = [row for row in decisions if row["gene"] == gene]
            kept = [row for row in gene_rows if row["keep"]]
            writer.writerow(
                {
                    "gene": gene,
                    "kept_total": len(kept),
                    "kept_ISS": sum(1 for row in kept if row["group"] == "ISS"),
                    "kept_Earth": sum(1 for row in kept if row["group"] == "Earth"),
                    "rejected_total": len(gene_rows) - len(kept),
                }
            )

    script_path = os.path.join(args.outdir, "run_alignment_relax_template.sh")
    with open(script_path, "w") as handle:
        handle.write(
            """#!/usr/bin/env bash
set -euo pipefail

# Required external tools: mafft, pal2nal.pl, iqtree2 or iqtree, hyphy.
# This template starts from the DIAMOND-checked FASTAs in this directory.

for gene in rfbA rfbB rfbC rfbD waaL; do
  mafft --auto "${gene}.relax_input.faa" > "${gene}.aa.aln.faa"
  pal2nal.pl "${gene}.aa.aln.faa" "${gene}.relax_input.ffn" -output fasta > "${gene}.codon.aln.fasta"

  if command -v iqtree2 >/dev/null 2>&1; then
    iqtree2 -s "${gene}.codon.aln.fasta" -m MFP -B 1000 -T AUTO
    tree="${gene}.codon.aln.fasta.treefile"
  else
    iqtree -s "${gene}.codon.aln.fasta" -m MFP -bb 1000 -nt AUTO
    tree="${gene}.codon.aln.fasta.treefile"
  fi

  python3 label_relax_tree.py --tree "${tree}" --map sequence_id_map.tsv --gene "${gene}" --out "${gene}.relax_labeled.tree"
  hyphy relax --alignment "${gene}.codon.aln.fasta" --tree "${gene}.relax_labeled.tree" --test ISS --reference Earth
done
"""
        )
    os.chmod(script_path, 0o755)

    labeler_path = os.path.join(args.outdir, "label_relax_tree.py")
    with open(labeler_path, "w") as handle:
        handle.write(
            """#!/usr/bin/env python3
import argparse
import csv
import re


def main():
    parser = argparse.ArgumentParser(description="Add HyPhy RELAX branch labels to terminal branches by ISS/Earth group.")
    parser.add_argument("--tree", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    groups = {}
    with open(args.map, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\\t")
        for row in reader:
            if row["gene"] == args.gene:
                groups[row["safe_id"]] = row["group"]

    with open(args.tree) as handle:
        tree = handle.read().strip()

    for leaf, group in sorted(groups.items(), key=lambda item: len(item[0]), reverse=True):
        tree = re.sub(rf"(?<![A-Za-z0-9_.-])({re.escape(leaf)})(?=[:),])", rf"\\1{{{group}}}", tree)

    with open(args.out, "w") as handle:
        handle.write(tree + "\\n")


if __name__ == "__main__":
    main()
"""
        )
    os.chmod(labeler_path, 0o755)

    print(f"Reference isolate: {ref_isolate}")
    print(f"Wrote Phase 2 files to {args.outdir}")


if __name__ == "__main__":
    main()
