#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import sys
import xml.etree.ElementTree as ET

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


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


def build_ko_name_map(annotations_root):
    name_counts = {}
    files = find_annotation_files(annotations_root)
    for path in files:
        ko_idx = None
        name_idx = None
        with open(path, "r", newline="") as handle:
            for line in handle:
                if line.startswith("#query"):
                    header = line.lstrip("#").rstrip("\n").split("\t")
                    try:
                        ko_idx = header.index("KEGG_ko")
                        name_idx = header.index("Preferred_name")
                    except ValueError:
                        break
                    continue
                if line.startswith("#"):
                    continue
                if ko_idx is None or name_idx is None:
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(ko_idx, name_idx):
                    continue
                ko_field = parts[ko_idx]
                if not ko_field or ko_field == "-":
                    continue
                pref = parts[name_idx].strip()
                if not pref or pref == "-":
                    continue
                for ko in ko_field.split(","):
                    ko = ko.strip()
                    if not ko:
                        continue
                    name_counts.setdefault(ko, {})
                    name_counts[ko][pref] = name_counts[ko].get(pref, 0) + 1
    ko_name = {}
    for ko, counts in name_counts.items():
        best = max(counts.items(), key=lambda item: (item[1], item[0]))
        ko_name[ko] = best[0]
    return ko_name


def find_annotation_files(root):
    pattern = os.path.join(root, "*", "eggnog_out", "*.emapper.annotations")
    return glob.glob(pattern)


def parse_emapper_for_kos(path, target_pathways):
    kos = set()
    pathway_idx = None
    ko_idx = None
    with open(path, "r", newline="") as handle:
        for line in handle:
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                try:
                    pathway_idx = header.index("KEGG_Pathway")
                    ko_idx = header.index("KEGG_ko")
                except ValueError:
                    return kos
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
            if not any(p in pathways.split(",") for p in target_pathways):
                continue
            ko_field = parts[ko_idx]
            if not ko_field or ko_field == "-":
                continue
            for ko in ko_field.split(","):
                ko = ko.strip()
                if ko:
                    kos.add(ko)
    return kos


def build_ko_group_counts(annotations_root, group_map, target_pathways):
    ko_counts = {}
    group_totals = {}
    files = find_annotation_files(annotations_root)
    for path in files:
        isolate = os.path.basename(path).split(".emapper.annotations")[0]
        group = group_map.get(isolate)
        if not group:
            continue
        group_totals[group] = group_totals.get(group, 0) + 1
        kos = parse_emapper_for_kos(path, target_pathways)
        for ko in kos:
            ko_counts.setdefault(ko, {}).setdefault(group, 0)
            ko_counts[ko][group] += 1
    return ko_counts, group_totals


def parse_kgml(path):
    tree = ET.parse(path)
    root = tree.getroot()
    entries = {}
    for entry in root.findall("entry"):
        entry_id = entry.get("id")
        entry_type = entry.get("type")
        name = entry.get("name", "")
        graphics = entry.find("graphics")
        if graphics is None:
            continue
        coords = graphics.get("coords")
        if coords:
            try:
                pts = [float(v) for v in coords.split(",")]
            except ValueError:
                continue
            xs = pts[0::2]
            ys = pts[1::2]
            if not xs or not ys:
                continue
            min_x, max_x = min(xs), max(xs)
            min_y, max_y = min(ys), max(ys)
            x = (min_x + max_x) / 2
            y = (min_y + max_y) / 2
            w = max(10.0, max_x - min_x)
            h = max(10.0, max_y - min_y)
        else:
            try:
                x = float(graphics.get("x"))
                y = float(graphics.get("y"))
                w = float(graphics.get("width"))
                h = float(graphics.get("height"))
            except (TypeError, ValueError):
                continue
        entries[entry_id] = {
            "type": entry_type,
            "name": name,
            "x": x,
            "y": y,
            "w": w,
            "h": h,
        }
    relations = []
    for rel in root.findall("relation"):
        e1 = rel.get("entry1")
        e2 = rel.get("entry2")
        if e1 and e2:
            relations.append((e1, e2))
    return entries, relations


def node_kos(entry_name):
    kos = []
    for token in entry_name.split():
        token = token.strip()
        if token.startswith("ko:"):
            kos.append(token)
    return kos


def plot_pathway(kgml_path, ko_counts, group_totals, groups, output_path, ko_name, background_png=None):
    entries, relations = parse_kgml(kgml_path)
    gene_entries = {
        eid: data
        for eid, data in entries.items()
        if data.get("type") in {"gene", "ortholog"}
    }
    if not gene_entries:
        print(f"Warning: no gene entries in {kgml_path}", file=sys.stderr)
        return

    max_x = max(data["x"] + data["w"] for data in gene_entries.values())
    max_y = max(data["y"] + data["h"] for data in gene_entries.values())
    fig_w = max(10, max_x / 50)
    fig_h = max(6, max_y / 50)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    if background_png:
        img = mpimg.imread(background_png)
        ax.imshow(img, extent=[0, max_x, max_y, 0], zorder=0)

    # Precompute max counts for color scaling.
    max_counts = {group: 1.0 for group in groups}
    for data in gene_entries.values():
        kos = node_kos(data["name"])
        for group in groups:
            total = 0
            for ko in kos:
                total += ko_counts.get(ko, {}).get(group, 0)
            denom = group_totals.get(group, 1)
            normalized = total / denom if denom else 0.0
            if normalized > max_counts[group]:
                max_counts[group] = normalized

    cmap = {
        groups[0]: cm.get_cmap("Blues"),
        groups[1]: cm.get_cmap("Reds"),
    }
    norm = colors.Normalize(vmin=0.0, vmax=1.0)

    if not background_png:
        for e1, e2 in relations:
            if e1 not in gene_entries or e2 not in gene_entries:
                continue
            x1, y1 = gene_entries[e1]["x"], gene_entries[e1]["y"]
            x2, y2 = gene_entries[e2]["x"], gene_entries[e2]["y"]
            ax.plot([x1, x2], [y1, y2], color="#888888", linewidth=0.6, zorder=1)

    for eid, data in gene_entries.items():
        x = data["x"] - data["w"] / 2
        y = data["y"] - data["h"] / 2
        w = data["w"]
        h = data["h"]
        kos = node_kos(data["name"])
        if not kos:
            continue
        counts = {}
        for group in groups:
            counts[group] = sum(ko_counts.get(ko, {}).get(group, 0) for ko in kos)
        normalized = {
            group: (counts[group] / group_totals.get(group, 1) if group_totals.get(group, 1) else 0.0)
            for group in groups
        }

        left_value = normalized[groups[0]]
        right_value = normalized[groups[1]]
        if left_value == 0.0 and right_value == 0.0:
            ax.add_patch(
                Rectangle(
                    (x, y),
                    w,
                    h,
                    facecolor="white",
                    edgecolor="#BBBBBB",
                    linewidth=0.5,
                    zorder=2,
                )
            )
        elif abs(left_value - right_value) < 1e-6:
            ax.add_patch(
                Rectangle(
                    (x, y),
                    w,
                    h,
                    facecolor="#B0B0B0",
                    edgecolor="black",
                    linewidth=0.5,
                    zorder=2,
                )
            )
        else:
            left_color = cmap[groups[0]](norm(left_value))
            right_color = cmap[groups[1]](norm(right_value))
            ax.add_patch(
                Rectangle(
                    (x, y),
                    w / 2,
                    h,
                    facecolor=left_color,
                    edgecolor="black",
                    linewidth=0.5,
                    alpha=0.8,
                    zorder=2,
                )
            )
            ax.add_patch(
                Rectangle(
                    (x + w / 2, y),
                    w / 2,
                    h,
                    facecolor=right_color,
                    edgecolor="black",
                    linewidth=0.5,
                    alpha=0.8,
                    zorder=2,
                )
            )

        if not background_png:
            ko = kos[0]
            label = ko.replace("ko:", "")
            if ko in ko_name:
                label = ko_name[ko]
            if len(kos) > 1:
                label = f"{label}(+{len(kos)-1})"
            ax.text(
                data["x"],
                data["y"] - h * 0.1,
                label,
                ha="center",
                va="center",
                fontsize=6,
                zorder=3,
            )
            ax.text(
                data["x"],
                data["y"] + h * 0.25,
                f"{groups[0][0]}:{normalized[groups[0]]:.2f} {groups[1][0]}:{normalized[groups[1]]:.2f}",
                ha="center",
                va="center",
                fontsize=5,
                zorder=3,
            )

    ax.set_xlim(0, max_x + 20)
    ax.set_ylim(max_y + 20, 0)
    ax.axis("off")
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Plot KEGG pathway networks with Earth/ISS KO presence counts."
    )
    parser.add_argument("--annotations-root", required=True, help="Root directory with genome folders")
    parser.add_argument("--group-map", required=True, help="TSV with isolate and group columns")
    parser.add_argument("--kgml", nargs="+", required=True, help="KGML pathway files")
    parser.add_argument("--out-dir", required=True, help="Output directory for pathway plots")
    parser.add_argument(
        "--png-dir",
        default=None,
        help="Directory containing KEGG pathway PNGs (e.g., map00540.png)",
    )
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
    ko_counts, group_totals = build_ko_group_counts(args.annotations_root, group_map, args.pathways)
    ko_name = build_ko_name_map(args.annotations_root)

    for kgml_path in args.kgml:
        base = os.path.splitext(os.path.basename(kgml_path))[0]
        out_path = os.path.join(args.out_dir, f"{base}_network.png")
        background = None
        if args.png_dir:
            png_name = base.replace("ko", "map") + ".png"
            candidate = os.path.join(args.png_dir, png_name)
            if os.path.exists(candidate):
                background = candidate
        plot_pathway(
            kgml_path,
            ko_counts,
            group_totals,
            args.groups,
            out_path,
            ko_name,
            background_png=background,
        )


if __name__ == "__main__":
    main()
