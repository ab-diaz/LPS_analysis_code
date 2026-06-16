#!/usr/bin/env python3
import csv
from collections import Counter
from pathlib import Path

from Bio import SeqIO


ROOT = Path("results/lps_analysis/selection_phase2")
OUT = ROOT / "relax_sensitivity_20260530"


def consensus(seqs):
    if not seqs:
        return ""
    length = len(seqs[0])
    out = []
    for i in range(length):
        chars = [s[i].upper() for s in seqs if s[i].upper() not in {"-", "N", "?"}]
        if not chars:
            out.append("N")
        else:
            out.append(Counter(chars).most_common(1)[0][0])
    return "".join(out)


def p_distance(a, b):
    total = 0
    diff = 0
    for x, y in zip(a.upper(), b.upper()):
        if x in {"-", "N", "?"} or y in {"-", "N", "?"}:
            continue
        total += 1
        if x != y:
            diff += 1
    return diff / total if total else 0


def main():
    divergent = list(csv.DictReader(open(OUT / "divergent_iss_tips.tsv"), delimiter="\t"))
    rows = []
    for row in divergent:
        gene = row["gene"]
        div_id = row["divergent_safe_id"]
        records = list(SeqIO.parse(ROOT / f"{gene}.codon.aln.fasta", "fasta"))
        seq_map = {r.id: str(r.seq).upper() for r in records}
        iss_records = [r for r in records if r.id.startswith(f"{gene}_ISS_")]
        iss_seqs = [str(r.seq).upper() for r in iss_records if r.id != div_id]
        cons = consensus(iss_seqs)
        seq = seq_map[div_id]
        gap = seq.count("-")
        amb = sum(seq.count(c) for c in ["N", "?"])
        rows.append(
            {
                "gene": gene,
                "divergent_safe_id": div_id,
                "isolate": row["isolate"],
                "alignment_nt_length": len(seq),
                "gap_count": gap,
                "ambiguous_count": amb,
                "gap_fraction": gap / len(seq) if seq else 0,
                "ambiguous_fraction": amb / len(seq) if seq else 0,
                "p_distance_to_other_ISS_consensus": p_distance(seq, cons),
            }
        )

    out = OUT / "divergent_sequence_qc.tsv"
    with open(out, "w", newline="") as handle:
        fields = list(rows[0])
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(out)


if __name__ == "__main__":
    main()
