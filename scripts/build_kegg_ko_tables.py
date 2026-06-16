#!/usr/bin/env python3
import argparse
import csv
import glob
import os


def load_group_map(path, isolate_col="isolate", group_col="group"):
    group_map = {}
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return group_map
        isolate_idx = header.index(isolate_col)
        group_idx = header.index(group_col)
        for row in reader:
            if len(row) <= max(isolate_idx, group_idx):
                continue
            isolate = row[isolate_idx].strip()
            group = row[group_idx].strip()
            if isolate and group:
                group_map[isolate] = group
    return group_map


def normalize_pathway_id(token):
    token = token.strip()
    if token.startswith("map"):
        return "ko" + token[3:]
    return token


def find_annotation_files(root):
    pattern = os.path.join(root, "*", "eggnog_out", "*.emapper.annotations")
    return glob.glob(pattern)


def parse_emapper_for_pathways(path, target_pathways):
    pathway_idx = None
    ko_idx = None
    data = {}
    with open(path, "r", newline="") as handle:
        for line in handle:
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                try:
                    pathway_idx = header.index("KEGG_Pathway")
                    ko_idx = header.index("KEGG_ko")
                except ValueError:
                    return data
                continue
            if line.startswith("#"):
                continue
            if pathway_idx is None or ko_idx is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(pathway_idx, ko_idx):
                continue
            pathways = parts[pathway_idx]
            if not pathways or pathways == "-":
                continue
            pathway_tokens = [normalize_pathway_id(p) for p in pathways.split(",")]
            pathway_hits = [p for p in pathway_tokens if p in target_pathways]
            if not pathway_hits:
                continue
            ko_field = parts[ko_idx]
            if not ko_field or ko_field == "-":
                continue
            kos = [ko.strip() for ko in ko_field.split(",") if ko.strip()]
            if not kos:
                continue
            for pathway in pathway_hits:
                data.setdefault(pathway, set()).update(kos)
    return data


def main():
    parser = argparse.ArgumentParser(
        description="Build KO tables per KEGG pathway with Earth/ISS counts and fractions."
    )
    parser.add_argument("--annotations-root", required=True, help="Root directory with genome folders")
    parser.add_argument("--group-map", required=True, help="TSV with isolate and group columns")
    parser.add_argument("--out-dir", required=True, help="Output directory for KO tables")
    parser.add_argument(
        "--pathways",
        nargs="+",
        default=["ko00540", "ko00541", "ko00520"],
        help="Target KEGG pathway IDs",
    )
    parser.add_argument(
        "--groups",
        nargs=2,
        default=["Earth", "ISS"],
        help="Two group names to compare (order matters)",
    )
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    group_map = load_group_map(args.group_map)

    group_totals = {group: 0 for group in args.groups}
    for isolate, group in group_map.items():
        if group in group_totals:
            group_totals[group] += 1

    counts = {pathway: {} for pathway in args.pathways}
    files = find_annotation_files(args.annotations_root)
    for path in files:
        isolate = os.path.basename(path).split(".emapper.annotations")[0]
        group = group_map.get(isolate)
        if group not in args.groups:
            continue
        per_pathway = parse_emapper_for_pathways(path, args.pathways)
        for pathway, kos in per_pathway.items():
            for ko in kos:
                counts[pathway].setdefault(ko, {g: 0 for g in args.groups})
                counts[pathway][ko][group] += 1

    for pathway, ko_map in counts.items():
        out_path = os.path.join(args.out_dir, f"{pathway}_ko.tsv")
        with open(out_path, "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(
                [
                    "ko",
                    f"{args.groups[0]}_count",
                    f"{args.groups[1]}_count",
                    f"{args.groups[0]}_frac",
                    f"{args.groups[1]}_frac",
                ]
            )
            for ko in sorted(ko_map.keys()):
                row = ko_map[ko]
                c1 = row.get(args.groups[0], 0)
                c2 = row.get(args.groups[1], 0)
                f1 = c1 / group_totals[args.groups[0]] if group_totals[args.groups[0]] else 0.0
                f2 = c2 / group_totals[args.groups[1]] if group_totals[args.groups[1]] else 0.0
                writer.writerow([ko.replace("ko:", ""), c1, c2, f1, f2])


if __name__ == "__main__":
    main()
