#!/usr/bin/env python3
import argparse
import csv
import os
from collections import defaultdict
from statistics import median
from urllib.parse import unquote


TARGETS = ("rfbA", "rfbB", "rfbC", "rfbD", "waaL")
RFB_TARGETS = {"rfbA", "rfbB", "rfbC", "rfbD"}
GENE_CANON = {g.lower(): g for g in TARGETS}


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


def load_groups(path):
    groups = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            groups[row["Folder"]] = row["Location"]
    return groups


def attrs_from_gff(field):
    attrs = {}
    for item in field.split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = unquote(value)
    return attrs


def load_gff_coords(root):
    coords = {}
    for isolate in os.listdir(root):
        gff_path = os.path.join(root, isolate, "bakta_output", f"{isolate}.gff3")
        if not os.path.exists(gff_path):
            continue
        with open(gff_path) as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9 or parts[2] != "CDS":
                    continue
                attrs = attrs_from_gff(parts[8])
                locus = attrs.get("locus_tag") or attrs.get("ID", "").replace("cds-", "")
                if not locus:
                    continue
                coords[(isolate, locus)] = {
                    "contig": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "strand": parts[6],
                    "product": attrs.get("product", ""),
                }
    return coords


def fasta_key(isolate, locus):
    return f"{isolate}|{locus}"


def add_candidate(candidates, isolate, gene, locus, source, groups, coords=None, extra=None):
    gene = GENE_CANON.get(gene.lower())
    if not gene:
        return
    group = groups.get(isolate, "UNKNOWN")
    item = {
        "isolate": isolate,
        "group": group,
        "gene": gene,
        "locus": locus,
        "source": source,
        "contig": "",
        "start": "",
        "end": "",
        "strand": "",
        "product": "",
        "description": "",
    }
    if coords:
        item.update(
            {
                "contig": coords.get("contig", ""),
                "start": coords.get("start", ""),
                "end": coords.get("end", ""),
                "strand": coords.get("strand", ""),
                "product": coords.get("product", ""),
            }
        )
    if extra:
        item.update({k: v for k, v in extra.items() if k in item})
    candidates.append(item)


def load_candidates(args, groups):
    candidates = []
    gff_coords = load_gff_coords(args.root)

    with open(args.rfb_locations, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            isolate = row["isolate"]
            if isolate not in groups:
                continue
            gene = GENE_CANON.get(row["gene"].lower())
            if gene not in RFB_TARGETS:
                continue
            add_candidate(
                candidates,
                isolate,
                gene,
                row["query"],
                "rfb_eggnog_locations",
                groups,
                coords={
                    "contig": row["contig"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "strand": row["strand"],
                    "product": row["product"],
                },
                extra={"description": row.get("description", "")},
            )

    with open(args.lps_hits, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            isolate = row["isolate"]
            if isolate not in groups:
                continue
            gene = GENE_CANON.get(row["preferred_name"].lower())
            if gene != "waaL":
                continue
            locus = row["query"]
            add_candidate(
                candidates,
                isolate,
                gene,
                locus,
                "lps_strict_hits",
                groups,
                coords=gff_coords.get((isolate, locus), {}),
                extra={"description": row.get("description", "")},
            )

    seen = set()
    unique = []
    for item in candidates:
        key = (item["isolate"], item["gene"], item["locus"])
        if key in seen:
            continue
        seen.add(key)
        unique.append(item)
    return unique


def distance(a, b):
    if not a.get("contig") or a.get("contig") != b.get("contig"):
        return None
    if not a.get("start") or not b.get("start"):
        return None
    a_start, a_end = int(a["start"]), int(a["end"])
    b_start, b_end = int(b["start"]), int(b["end"])
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0


def cluster_score(candidate, isolate_candidates, max_gap):
    distances = []
    for other in isolate_candidates:
        if other is candidate or other["gene"] == candidate["gene"]:
            continue
        if other["gene"] not in RFB_TARGETS:
            continue
        dist = distance(candidate, other)
        if dist is not None and dist <= max_gap:
            distances.append(dist)
    if not distances:
        return 0, 10**12
    return len(distances), sum(distances)


def curate(candidates, nuc_records, prot_records, groups, min_frac, max_frac, max_gap):
    by_gene = defaultdict(list)
    for item in candidates:
        key = fasta_key(item["isolate"], item["locus"])
        item["nuc_seq"] = nuc_records.get(key, "")
        item["prot_seq"] = prot_records.get(key, "")
        item["nt_len"] = len(item["nuc_seq"])
        item["aa_len"] = len(item["prot_seq"].rstrip("*"))
        item["has_sequence"] = bool(item["nuc_seq"] and item["prot_seq"])
        item["nt_multiple_of_3"] = item["nt_len"] % 3 == 0 if item["nuc_seq"] else False
        by_gene[item["gene"]].append(item)

    medians = {}
    for gene, items in by_gene.items():
        iss_lengths = [
            item["aa_len"]
            for item in items
            if item["group"] == "ISS" and item["has_sequence"] and item["aa_len"] > 0
        ]
        all_lengths = [
            item["aa_len"] for item in items if item["has_sequence"] and item["aa_len"] > 0
        ]
        medians[gene] = median(iss_lengths or all_lengths) if (iss_lengths or all_lengths) else 0

    for item in candidates:
        med = medians.get(item["gene"], 0)
        item["median_aa_len"] = med
        item["length_min"] = int(med * min_frac) if med else 0
        item["length_max"] = int(med * max_frac) if med else 0
        item["length_pass"] = (
            item["has_sequence"]
            and item["nt_multiple_of_3"]
            and item["aa_len"] >= item["length_min"]
            and item["aa_len"] <= item["length_max"]
        )
        item["decision"] = "candidate"
        item["reason"] = ""
        if not item["has_sequence"]:
            item["decision"] = "reject"
            item["reason"] = "missing_sequence"
        elif not item["nt_multiple_of_3"]:
            item["decision"] = "reject"
            item["reason"] = "cds_length_not_multiple_of_3"
        elif not item["length_pass"]:
            item["decision"] = "reject"
            item["reason"] = "length_outlier"

    accepted = []
    by_isolate_gene = defaultdict(list)
    by_isolate = defaultdict(list)
    for item in candidates:
        by_isolate[item["isolate"]].append(item)
        if item["length_pass"]:
            by_isolate_gene[(item["isolate"], item["gene"])].append(item)

    for (isolate, gene), items in by_isolate_gene.items():
        if len(items) == 1:
            chosen = items[0]
            chosen["decision"] = "keep"
            chosen["reason"] = "single_copy_after_length_filter"
            accepted.append(chosen)
            continue

        if gene in RFB_TARGETS:
            scored = []
            for item in items:
                neighbor_count, total_distance = cluster_score(item, by_isolate[isolate], max_gap)
                item["cluster_neighbors"] = neighbor_count
                item["cluster_distance"] = total_distance
                length_delta = abs(item["aa_len"] - item["median_aa_len"])
                scored.append((-neighbor_count, total_distance, length_delta, item["locus"], item))
            chosen = sorted(scored)[0][-1]
            reason = "duplicate_resolved_by_operon_proximity"
        else:
            chosen = sorted(
                items,
                key=lambda item: (
                    abs(item["aa_len"] - item["median_aa_len"]),
                    item["locus"],
                ),
            )[0]
            reason = "duplicate_resolved_by_length_closest_to_median"

        for item in items:
            if item is chosen:
                item["decision"] = "keep"
                item["reason"] = reason
                accepted.append(item)
            else:
                item["decision"] = "reject"
                item["reason"] = "duplicate_nonselected"

    return accepted, medians


def sort_key(item):
    return (item["gene"], item["group"], item["isolate"], item["locus"])


def main():
    parser = argparse.ArgumentParser(
        description="Curate single-copy rfbABCD and waaL orthologs for codon selection analysis."
    )
    parser.add_argument("--root", default=".")
    parser.add_argument("--groups", default="old/pantoea_folders_location.tsv")
    parser.add_argument("--nucleotides", default="lps_analysis/lps_eggnog_nucleotides.ffn")
    parser.add_argument("--proteins", default="lps_analysis/lps_eggnog_proteins.faa")
    parser.add_argument("--rfb-locations", default="lps_analysis/rfb_eggnog_locations.tsv")
    parser.add_argument("--lps-hits", default="lps_analysis/lps_strict_hits.tsv")
    parser.add_argument("--outdir", default="lps_analysis/selection_curation")
    parser.add_argument("--min-frac", type=float, default=0.80)
    parser.add_argument("--max-frac", type=float, default=1.20)
    parser.add_argument("--max-operon-gap", type=int, default=20000)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    groups = load_groups(args.groups)
    nuc_records = read_fasta(args.nucleotides)
    prot_records = read_fasta(args.proteins)
    candidates = load_candidates(args, groups)
    accepted, medians = curate(
        candidates,
        nuc_records,
        prot_records,
        groups,
        args.min_frac,
        args.max_frac,
        args.max_operon_gap,
    )

    for gene in TARGETS:
        gene_items = sorted([item for item in accepted if item["gene"] == gene], key=sort_key)
        nuc_out = [
            (
                f"{item['isolate']}|{item['group']}|{gene}|{item['locus']}",
                item["nuc_seq"],
            )
            for item in gene_items
        ]
        prot_out = [
            (
                f"{item['isolate']}|{item['group']}|{gene}|{item['locus']}",
                item["prot_seq"],
            )
            for item in gene_items
        ]
        write_fasta(nuc_out, os.path.join(args.outdir, f"{gene}.ffn"))
        write_fasta(prot_out, os.path.join(args.outdir, f"{gene}.faa"))

    decision_fields = [
        "isolate",
        "group",
        "gene",
        "locus",
        "decision",
        "reason",
        "aa_len",
        "nt_len",
        "median_aa_len",
        "length_min",
        "length_max",
        "contig",
        "start",
        "end",
        "strand",
        "source",
        "product",
        "description",
    ]
    with open(os.path.join(args.outdir, "curation_decisions.tsv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=decision_fields)
        writer.writeheader()
        for item in sorted(candidates, key=sort_key):
            writer.writerow({field: item.get(field, "") for field in decision_fields})

    count_fields = [
        "gene",
        "median_aa_len",
        "candidate_total",
        "kept_total",
        "kept_ISS",
        "kept_Earth",
        "rejected_missing_sequence",
        "rejected_length_outlier",
        "rejected_duplicate_nonselected",
        "rejected_other",
    ]
    with open(os.path.join(args.outdir, "curation_summary.tsv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=count_fields)
        writer.writeheader()
        for gene in TARGETS:
            gene_candidates = [item for item in candidates if item["gene"] == gene]
            kept = [item for item in gene_candidates if item["decision"] == "keep"]
            reasons = defaultdict(int)
            for item in gene_candidates:
                if item["decision"] != "keep":
                    reasons[item["reason"] or "other"] += 1
            writer.writerow(
                {
                    "gene": gene,
                    "median_aa_len": medians.get(gene, ""),
                    "candidate_total": len(gene_candidates),
                    "kept_total": len(kept),
                    "kept_ISS": sum(1 for item in kept if item["group"] == "ISS"),
                    "kept_Earth": sum(1 for item in kept if item["group"] == "Earth"),
                    "rejected_missing_sequence": reasons["missing_sequence"],
                    "rejected_length_outlier": reasons["length_outlier"],
                    "rejected_duplicate_nonselected": reasons["duplicate_nonselected"],
                    "rejected_other": sum(
                        count
                        for reason, count in reasons.items()
                        if reason
                        not in {
                            "missing_sequence",
                            "length_outlier",
                            "duplicate_nonselected",
                        }
                    ),
                }
            )

    print(f"Wrote curated FASTA and QC tables to {args.outdir}")


if __name__ == "__main__":
    main()
