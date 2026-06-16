#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re

TARGET_GENES = {"rfba", "rfbb", "rfbc", "rfbd"}
GENE_RE = re.compile(r"\brfb[abcd]\b", re.IGNORECASE)


def parse_attributes(attr_str):
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


def collect_hits(gff3_path, isolate):
    hits = []
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
            attr_map = parse_attributes(attrs)
            gene = attr_map.get("gene", "")
            name = attr_map.get("Name", "")
            product = attr_map.get("product", "")

            text = " ".join([gene, name, product])
            match = GENE_RE.search(text)
            if not match:
                continue

            gene_name = match.group(0).lower()
            if gene_name not in TARGET_GENES:
                continue

            locus_tag = attr_map.get("locus_tag", attr_map.get("ID", ""))
            hits.append(
                {
                    "isolate": isolate,
                    "gene": gene_name,
                    "contig": contig,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "locus_tag": locus_tag,
                    "name": name,
                    "product": product,
                }
            )
    return hits


def group_hits(hits):
    groups = {}
    for h in hits:
        key = (h["contig"], h["strand"])
        groups.setdefault(key, []).append(h)
    return groups


def order_and_gaps(group_hits):
    ordered = sorted(group_hits, key=lambda x: x["start"])
    gaps = []
    for a, b in zip(ordered, ordered[1:]):
        gap = b["start"] - a["end"] - 1
        gaps.append(gap)
    return ordered, gaps


def summarize_isolate(isolate, hits, gap_threshold):
    summary_rows = []
    if not hits:
        summary_rows.append(
            {
                "isolate": isolate,
                "contig": ".",
                "strand": ".",
                "genes_found": "",
                "has_all_four": "0",
                "order": "",
                "order_matches": "0",
                "max_gap_bp": "",
                "span_bp": "",
                "gaps_bp": "",
                "notes": "no_hits",
            }
        )
        return summary_rows

    groups = group_hits(hits)
    any_full = False
    for (contig, strand), ghits in groups.items():
        found = sorted({h["gene"] for h in ghits})
        has_all = TARGET_GENES.issubset(found)
        if not has_all:
            continue
        any_full = True
        ordered, gaps = order_and_gaps(ghits)
        order = [h["gene"] for h in ordered]
        order_str = ">".join(order)
        max_gap = max(gaps) if gaps else 0
        span = max(h["end"] for h in ordered) - min(h["start"] for h in ordered) + 1
        order_matches = order == ["rfba", "rfbb", "rfbc", "rfbd"]
        operon_like = (max_gap <= gap_threshold) and order_matches

        summary_rows.append(
            {
                "isolate": isolate,
                "contig": contig,
                "strand": strand,
                "genes_found": ",".join(found),
                "has_all_four": "1",
                "order": order_str,
                "order_matches": "1" if order_matches else "0",
                "max_gap_bp": str(max_gap),
                "span_bp": str(span),
                "gaps_bp": ",".join(str(g) for g in gaps),
                "notes": "operon_like" if operon_like else "",
            }
        )

    if not any_full:
        found = sorted({h["gene"] for h in hits})
        summary_rows.append(
            {
                "isolate": isolate,
                "contig": ".",
                "strand": ".",
                "genes_found": ",".join(found),
                "has_all_four": "0",
                "order": "",
                "order_matches": "0",
                "max_gap_bp": "",
                "span_bp": "",
                "gaps_bp": "",
                "notes": "partial_set",
            }
        )

    return summary_rows


def main():
    parser = argparse.ArgumentParser(
        description="Find rfbA-D co-localization from Bakta GFF3 files."
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
        "--out-summary",
        default="results/lps_analysis/rfb_operon_summary.tsv",
        help="Output TSV summary.",
    )
    parser.add_argument(
        "--out-hits",
        default="results/lps_analysis/rfb_operon_hits.tsv",
        help="Output TSV with all rfbA-D hits.",
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
    all_summary = []

    for sdir in sample_dirs:
        isolate = os.path.basename(sdir)
        gff3_path = find_gff3(sdir)
        if not gff3_path:
            continue
        hits = collect_hits(gff3_path, isolate)
        if hits:
            all_hits.extend(hits)
        all_summary.extend(summarize_isolate(isolate, hits, args.gap))

    os.makedirs(os.path.dirname(args.out_summary), exist_ok=True)

    with open(args.out_hits, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "isolate",
                "gene",
                "contig",
                "start",
                "end",
                "strand",
                "locus_tag",
                "name",
                "product",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in all_hits:
            writer.writerow(row)

    with open(args.out_summary, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "isolate",
                "contig",
                "strand",
                "genes_found",
                "has_all_four",
                "order",
                "order_matches",
                "max_gap_bp",
                "span_bp",
                "gaps_bp",
                "notes",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in all_summary:
            writer.writerow(row)


if __name__ == "__main__":
    main()
