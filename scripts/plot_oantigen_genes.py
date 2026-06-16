#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re

import matplotlib.pyplot as plt
import pandas as pd


KEYWORDS = [
    "o-antigen",
    "o antigen",
    "o-antigen ligase",
    "o-antigen flippase",
    "o-antigen polymerase",
    "chain length determinant",
]

GENE_PREFIXES = [
    "rfa",
    "rfb",
]

GENE_NAMES = [
    "wzx",
    "wzy",
    "wzz",
    "wzt",
    "wzm",
    "waal",
    "tolA",
    "tolB",
    "tolQ",
    "tolR",
]


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


def find_annotation_files(root):
    pattern = os.path.join(root, "*", "eggnog_out", "*.emapper.annotations")
    return glob.glob(pattern)


def build_gene_regex():
    parts = []
    parts.extend([re.escape(k) for k in GENE_NAMES])
    parts.extend([f"{p}[a-z0-9_]+" for p in GENE_PREFIXES])
    return re.compile(r"\b(" + "|".join(parts) + r")\b", re.IGNORECASE)


def extract_gene_label(preferred, description, gene_re):
    if preferred and preferred != "-":
        pref = preferred.strip()
        if gene_re.search(pref):
            return pref
    text = " ".join([preferred or "", description or ""]).lower()
    match = gene_re.search(text)
    if match:
        return match.group(1)
    return None


def matches_keyword(description):
    if not description:
        return False
    lower = description.lower()
    return any(k in lower for k in KEYWORDS)


def collect_ko_presence(annotations_root, group_map):
    gene_re = build_gene_regex()
    presence = {}
    ko_name_counts = {}
    files = find_annotation_files(annotations_root)
    for path in files:
        isolate = os.path.basename(path).split(".emapper.annotations")[0]
        if isolate not in group_map:
            continue
        preferred_idx = None
        desc_idx = None
        ko_idx = None
        with open(path, "r", newline="") as handle:
            for line in handle:
                if line.startswith("#query"):
                    header = line.lstrip("#").rstrip("\n").split("\t")
                    if "Preferred_name" in header:
                        preferred_idx = header.index("Preferred_name")
                    if "Description" in header:
                        desc_idx = header.index("Description")
                    if "KEGG_ko" in header:
                        ko_idx = header.index("KEGG_ko")
                    break
        if preferred_idx is None or desc_idx is None or ko_idx is None:
            continue
        with open(path, "r", newline="") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(preferred_idx, desc_idx, ko_idx):
                    continue
                preferred = parts[preferred_idx]
                description = parts[desc_idx]
                ko_field = parts[ko_idx]
                if not ko_field or ko_field == "-":
                    continue
                label = extract_gene_label(preferred, description, gene_re)
                if label is None and not matches_keyword(description):
                    continue
                for ko in [k.strip() for k in ko_field.split(",") if k.strip()]:
                    presence.setdefault(ko, set()).add(isolate)
                    if preferred and preferred != "-":
                        ko_name_counts.setdefault(ko, {})
                        ko_name_counts[ko][preferred] = ko_name_counts[ko].get(preferred, 0) + 1
    ko_names = {}
    for ko, counts in ko_name_counts.items():
        ko_names[ko] = max(counts.items(), key=lambda item: (item[1], item[0]))[0]
    return presence, ko_names
    return presence


def build_fraction_table(presence, group_map, groups, ko_names):
    totals = {g: 0 for g in groups}
    for isolate, group in group_map.items():
        if group in totals:
            totals[group] += 1
    rows = []
    for ko in sorted(presence.keys()):
        per_group = {g: 0 for g in groups}
        for isolate in presence[ko]:
            group = group_map.get(isolate)
            if group in per_group:
                per_group[group] += 1
        label = ko.replace("ko:", "")
        if ko in ko_names:
            label = f"{label} {ko_names[ko]}"
        row = {
            "gene": label,
            f"{groups[0]}_count": per_group[groups[0]],
            f"{groups[1]}_count": per_group[groups[1]],
            f"{groups[0]}_frac": per_group[groups[0]] / totals[groups[0]] if totals[groups[0]] else 0.0,
            f"{groups[1]}_frac": per_group[groups[1]] / totals[groups[1]] if totals[groups[1]] else 0.0,
        }
        rows.append(row)
    return pd.DataFrame(rows)


def plot_split_bar(df, groups, output_path, top_n=None):
    df = df.sort_values(by=[f"{groups[0]}_frac", f"{groups[1]}_frac"], ascending=False)
    if top_n:
        df = df.head(top_n)
    y = range(len(df))
    left_vals = -df[f"{groups[0]}_frac"].values
    right_vals = df[f"{groups[1]}_frac"].values
    fig, ax = plt.subplots(figsize=(9, max(6, len(df) * 0.25)))
    ax.barh(y, left_vals, color="#4F81BD", alpha=0.85, label=groups[0])
    ax.barh(y, right_vals, color="#C0504D", alpha=0.85, label=groups[1])
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"])
    ax.set_xlabel("Fraction of genomes")
    ax.set_xlim(-1, 1)
    ax.axvline(0, color="#333333", linewidth=0.8)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=len(groups),
        frameon=False,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Plot split bar chart of O-antigen gene presence by group."
    )
    parser.add_argument("--annotations-root", required=True, help="Root directory with genome folders")
    parser.add_argument("--group-map", required=True, help="TSV with isolate and group columns")
    parser.add_argument("--out-tsv", required=True, help="Output TSV with gene fractions")
    parser.add_argument("--out-plot", required=True, help="Output plot path (png/pdf)")
    parser.add_argument("--groups", nargs=2, default=["Earth", "ISS"])
    parser.add_argument("--top-n", type=int, default=None)
    args = parser.parse_args()

    group_map = load_group_map(args.group_map)
    presence, ko_names = collect_ko_presence(args.annotations_root, group_map)
    df = build_fraction_table(presence, group_map, args.groups, ko_names)
    df.to_csv(args.out_tsv, sep="\t", index=False)
    plot_split_bar(df, args.groups, args.out_plot, top_n=args.top_n)


if __name__ == "__main__":
    main()
