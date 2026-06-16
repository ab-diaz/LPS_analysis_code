#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd


def read_lps_terms(path):
    terms = {}
    with open(path, "r", newline="") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            go_id = parts[0].strip() if parts else ""
            if not go_id or not go_id.startswith("GO:"):
                continue
            name = parts[1].strip() if len(parts) > 1 else ""
            terms[go_id] = name
    return terms


def load_lps_hits(path):
    isolates = {}
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return isolates
        try:
            isolate_idx = header.index("isolate")
            query_idx = header.index("query")
        except ValueError:
            raise ValueError("lps_hits.tsv must include 'isolate' and 'query' columns.")
        for row in reader:
            if len(row) <= max(isolate_idx, query_idx):
                continue
            isolate = row[isolate_idx].strip()
            query = row[query_idx].strip()
            if not isolate or not query:
                continue
            isolates.setdefault(isolate, set()).add(query)
    return isolates


def find_annotation_file(root, isolate):
    pattern = os.path.join(root, isolate, "eggnog_out", "*.emapper.annotations")
    matches = glob.glob(pattern)
    if not matches:
        return None
    return matches[0]


def parse_emapper_annotations(path, target_queries, target_go_terms):
    total_genes = 0
    query_to_gos = {}
    gos_idx = None
    query_idx = None
    with open(path, "r", newline="") as handle:
        for line in handle:
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                try:
                    query_idx = header.index("query")
                    gos_idx = header.index("GOs")
                except ValueError:
                    raise ValueError(f"Missing query/GOs columns in {path}")
                continue
            if line.startswith("#"):
                continue
            if gos_idx is None or query_idx is None:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(query_idx, gos_idx):
                continue
            total_genes += 1
            query = parts[query_idx]
            if query not in target_queries:
                continue
            gos_field = parts[gos_idx]
            if not gos_field or gos_field == "-":
                continue
            for go in gos_field.split(","):
                go = go.strip()
                if go in target_go_terms:
                    query_to_gos.setdefault(query, []).append(go)
    return total_genes, query_to_gos


def build_counts(isolates, root, target_go_terms, count_mode):
    counts = {}
    totals = {}
    lps_totals = {}
    for isolate, queries in isolates.items():
        ann_path = find_annotation_file(root, isolate)
        if not ann_path:
            print(f"Warning: no annotations for {isolate}", file=sys.stderr)
            continue
        total_genes, query_to_gos = parse_emapper_annotations(
            ann_path, queries, target_go_terms
        )
        totals[isolate] = total_genes
        lps_totals[isolate] = len(queries)
        if total_genes == 0:
            continue
        if count_mode == "gene":
            for gos in query_to_gos.values():
                for go in set(gos):
                    counts.setdefault(isolate, {}).setdefault(go, 0)
                    counts[isolate][go] += 1
        else:
            for gos in query_to_gos.values():
                for go in gos:
                    counts.setdefault(isolate, {}).setdefault(go, 0)
                    counts[isolate][go] += 1
    return counts, totals, lps_totals


def build_dataframe(counts, totals):
    df = pd.DataFrame.from_dict(counts, orient="index").fillna(0).astype(float)
    for isolate, total in totals.items():
        if total > 0 and isolate in df.index:
            df.loc[isolate] = df.loc[isolate] / total
    df = df.sort_index()
    return df


def parse_rename_rules(rules):
    renames = {}
    for rule in rules or []:
        if "=" not in rule:
            raise ValueError(f"Rename rule must use OLD=NEW format: {rule}")
        old, new = rule.split("=", 1)
        old = old.strip()
        new = new.strip()
        if not old or not new:
            raise ValueError(f"Rename rule must include both OLD and NEW: {rule}")
        renames[old] = new
    return renames


def apply_renames(df, renames):
    if not renames:
        return df
    collisions = [
        new for old, new in renames.items() if old in df.index and new in df.index and old != new
    ]
    if collisions:
        raise ValueError(
            "Rename would create duplicate isolate labels: " + ", ".join(sorted(collisions))
        )
    return df.rename(index=renames)


def sniff_delimiter(path):
    with open(path, "r", newline="") as handle:
        sample = handle.read(4096)
    try:
        return csv.Sniffer().sniff(sample, delimiters="\t,").delimiter
    except csv.Error:
        return "\t"


def load_metadata(path, renames=None):
    delimiter = sniff_delimiter(path)
    metadata = {}
    with open(path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            return metadata
        fieldnames = {name.strip(): name for name in reader.fieldnames}
        sample_col = fieldnames.get("Sample ID")
        if not sample_col:
            raise ValueError("Metadata file must include a 'Sample ID' column.")
        for row in reader:
            sample = row.get(sample_col, "").strip()
            if not sample:
                continue
            normalized = {key.strip(): value.strip() for key, value in row.items() if key}
            metadata[sample] = normalized
    for old, new in (renames or {}).items():
        if new in metadata:
            continue
        if old in metadata:
            metadata[new] = dict(metadata[old])
            metadata[new]["Sample ID"] = new
    return metadata


def parse_flight_rank(isolate):
    if isolate.startswith("Pantoea_piersonii_IIIF"):
        return (3, isolate)
    match = re.search(r"_F(\d+)", isolate)
    if match:
        return (int(match.group(1)), isolate)
    return (999, isolate)


def metadata_group(isolate, metadata):
    meta = metadata.get(isolate, {})
    category = meta.get("Category", "")
    location = meta.get("Location", "")
    is_piersonii = "piersonii" in isolate
    if not is_piersonii:
        return "Non-piersonii"
    if location == "ISS":
        return "ISS"
    if category == "Clinical":
        return "Earth clinical piersonii"
    return "Earth non-clinical piersonii"


def metadata_sort_key(isolate, metadata, iss_order):
    meta = metadata.get(isolate, {})
    category = meta.get("Category", "")
    location = meta.get("Location", "")
    source = meta.get("Source", "")
    is_piersonii = "piersonii" in isolate

    if not is_piersonii:
        return (0, category, isolate)
    if location != "ISS" and category != "Clinical":
        return (1, category, isolate)
    if location != "ISS" and category == "Clinical":
        return (2, isolate)
    if location == "ISS":
        if iss_order == "flight":
            return (3, *parse_flight_rank(isolate))
        return (3, source, isolate)
    return (4, isolate)


def order_by_metadata(df, metadata, iss_order):
    missing = [isolate for isolate in df.index if isolate not in metadata]
    if missing:
        print(
            "Warning: no metadata for " + ", ".join(sorted(missing)),
            file=sys.stderr,
        )
    ordered = sorted(df.index, key=lambda isolate: metadata_sort_key(isolate, metadata, iss_order))
    return df.loc[ordered]


def build_metadata_group_map(isolates, metadata):
    return {isolate: metadata_group(isolate, metadata) for isolate in isolates}


def load_column_order(path):
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
    if not header:
        return []
    return [column for column in header[1:] if column]


def apply_column_order(df, columns, keep_only=False):
    if not columns:
        return df
    ordered = [column for column in columns if column in df.columns]
    if keep_only:
        return df.loc[:, ordered]
    extras = [column for column in df.columns if column not in ordered]
    return df.loc[:, ordered + extras]


def plot_stacked_bar(
    df, output_path, labels, top_n=None, group_map=None, ylabel=None, sort_columns=True
):
    if top_n is not None and top_n > 0:
        totals = df.sum(axis=0).sort_values(ascending=False)
        keep = totals.head(top_n).index
        df = df[keep]
    if sort_columns:
        df = df.loc[:, df.sum(axis=0).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(22, 7))
    df.plot(kind="bar", stacked=True, ax=ax, width=0.9)
    ax.set_ylabel(ylabel or "Normalized count", labelpad=10)
    ax.set_xlabel("Genome")
    ax.tick_params(axis="x", labelsize=6)
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right")
    legend_labels = [
        f"{go} {labels.get(go, '')}".strip() for go in df.columns.tolist()
    ]
    legend = ax.legend(
        legend_labels,
        title="GO term",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )
    extra_artists = [legend]
    if group_map:
        handles = color_ticklabels_by_group(ax, group_map)
        if handles:
            group_legend = ax.legend(
                handles=handles,
                title="Group",
                bbox_to_anchor=(1.02, -0.2),
                loc="upper left",
                ncol=1,
                frameon=False,
            )
            ax.add_artist(legend)
            extra_artists.append(group_legend)
    fig.subplots_adjust(left=0.12, bottom=0.35, right=0.6)
    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
        bbox_extra_artists=extra_artists,
        pad_inches=0.5,
    )
    plt.close(fig)


def load_group_map(path, isolate_col="isolate", group_col="group"):
    group_map = {}
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if not header:
            return group_map
        try:
            isolate_idx = header.index(isolate_col)
            group_idx = header.index(group_col)
        except ValueError as exc:
            raise ValueError("Group file must include isolate and group columns.") from exc
        for row in reader:
            if len(row) <= max(isolate_idx, group_idx):
                continue
            isolate = row[isolate_idx].strip()
            group = row[group_idx].strip()
            if isolate and group:
                group_map[isolate] = group
    return group_map


def order_by_group(df, group_map):
    return df.loc[sorted(df.index, key=lambda iso: (group_map.get(iso, "ZZZ"), iso))]


def color_ticklabels_by_group(ax, group_map):
    groups = sorted(
        {
            group_map.get(label.get_text())
            for label in ax.get_xticklabels()
            if group_map.get(label.get_text())
        }
    )
    if not groups:
        return None
    cmap = plt.get_cmap("Dark2")
    group_colors = {group: cmap(i % cmap.N) for i, group in enumerate(groups)}
    for label in ax.get_xticklabels():
        group = group_map.get(label.get_text())
        if group:
            label.set_color(group_colors.get(group, "black"))
    handles = [Patch(color=group_colors[g], label=g) for g in groups]
    return handles


def main():
    parser = argparse.ArgumentParser(
        description="Plot stacked bar of LPS_terms GO IDs for LPS hits normalized by total genes."
    )
    parser.add_argument("--lps-hits", required=True, help="Path to lps_hits.tsv")
    parser.add_argument("--annotations-root", required=True, help="Root directory with genome folders")
    parser.add_argument("--lps-terms", required=True, help="TSV with GO IDs (e.g., LPS_terms.tsv)")
    parser.add_argument("--out-tsv", required=True, help="Output TSV for normalized counts")
    parser.add_argument("--out-plot", required=True, help="Output plot path (png/pdf)")
    parser.add_argument("--top-n", type=int, default=None, help="Keep only top N GO terms")
    parser.add_argument(
        "--normalize-by",
        choices=["total", "lps"],
        default="total",
        help="Normalize by total genes or LPS gene count",
    )
    parser.add_argument(
        "--count-mode",
        choices=["go", "gene"],
        default="go",
        help="Count GO term per GO occurrence or per gene presence",
    )
    parser.add_argument(
        "--group-map",
        default=None,
        help="TSV with isolate and group columns to order/label genomes",
    )
    parser.add_argument(
        "--metadata",
        default=None,
        help="Metadata CSV/TSV with Sample ID, Location, Source, and Category columns.",
    )
    parser.add_argument(
        "--iss-order",
        choices=["location", "flight"],
        default="location",
        help="Order ISS isolates by metadata Source/location or parsed flight label.",
    )
    parser.add_argument(
        "--rename-isolate",
        action="append",
        default=[],
        metavar="OLD=NEW",
        help="Rename an isolate label after counts are computed; may be supplied multiple times.",
    )
    parser.add_argument(
        "--column-order-from",
        default=None,
        help="Existing TSV whose header should be used to order GO-term columns.",
    )
    parser.add_argument(
        "--only-columns-from",
        default=None,
        help="Existing TSV whose header should be used as the complete GO-term column set.",
    )
    args = parser.parse_args()

    renames = parse_rename_rules(args.rename_isolate)
    terms = read_lps_terms(args.lps_terms)
    if not terms:
        raise ValueError("LPS_terms file appears empty or has no GO IDs.")
    isolates = load_lps_hits(args.lps_hits)
    counts, totals, lps_totals = build_counts(
        isolates, args.annotations_root, set(terms.keys()), args.count_mode
    )
    denom = totals if args.normalize_by == "total" else lps_totals
    df = build_dataframe(counts, denom)
    df = apply_renames(df, renames)
    group_map = None
    if args.metadata:
        metadata = load_metadata(args.metadata, renames)
        df = order_by_metadata(df, metadata, args.iss_order)
        group_map = build_metadata_group_map(df.index, metadata)
    elif args.group_map:
        group_map = load_group_map(args.group_map)
        if group_map:
            df = order_by_group(df, group_map)
    if args.only_columns_from:
        df = apply_column_order(df, load_column_order(args.only_columns_from), keep_only=True)
    elif args.column_order_from:
        df = apply_column_order(df, load_column_order(args.column_order_from))
    df.to_csv(args.out_tsv, sep="\t", index=True)
    ylabel = (
        "Normalized count per total genes"
        if args.normalize_by == "total"
        else "Normalized count per LPS genes"
    )
    plot_stacked_bar(
        df,
        args.out_plot,
        labels=terms,
        top_n=args.top_n,
        group_map=group_map,
        ylabel=ylabel,
        sort_columns=not (args.column_order_from or args.only_columns_from),
    )


if __name__ == "__main__":
    main()
