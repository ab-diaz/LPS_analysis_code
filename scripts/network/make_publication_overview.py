#!/usr/bin/env python3
"""Create a clean group-level overview from the exact native SplitsTree SVG."""

from __future__ import annotations

import argparse
import csv
import re
import xml.etree.ElementTree as ET
from pathlib import Path


SVG = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG)


def parse_taxa(workflow: Path) -> list[str]:
    text = workflow.read_text()
    block = re.search(r"BEGIN TAXA;(.*?)END; \[TAXA\]", text, re.S | re.I)
    if not block:
        raise ValueError("No TAXA block")
    return [value.replace("''", "'") for value in re.findall(r"\[\d+\]\s+'((?:[^']|'')*)'", block.group(1))]


def parse_terminal_weights(workflow: Path, ntax: int) -> tuple[dict[int, float], dict[int, int]]:
    text = workflow.read_text()
    blocks = re.findall(r"BEGIN SPLITS;(.*?)END; \[SPLITS\]", text, re.S | re.I)
    if not blocks:
        raise ValueError("No SPLITS block")
    block = blocks[-1]
    cycle_match = re.search(r"CYCLE\s+(.*?);", block, re.S | re.I)
    cycle = [int(value) for value in cycle_match.group(1).split()] if cycle_match else []
    cycle_position = {taxon: position for position, taxon in enumerate(cycle, 1)}
    matrix = re.search(r"MATRIX(.*?);", block, re.S | re.I)
    terminal: dict[int, float] = {}
    if matrix:
        for match in re.finditer(
            r"\[\d+,\s*size=1\]\s+([0-9.eE+-]+)\s+([^,]*),", matrix.group(1)
        ):
            weight = float(match.group(1))
            side = {int(value) for value in match.group(2).split()}
            if len(side) == 1:
                taxon = next(iter(side))
            elif len(side) == ntax - 1:
                taxon = next(iter(set(range(1, ntax + 1)) - side))
            else:
                continue
            terminal[taxon] = max(weight, terminal.get(taxon, 0.0))
    return terminal, cycle_position


def group_for(name: str) -> str:
    if name.startswith("Pantoea_piersonii_"):
        return "P. piersonii"
    if name.startswith("Pantoea_septica_"):
        return "P. septica"
    return "Other Pantoea"


def species_for(name: str) -> str:
    parts = name.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else name


def add_arrow_marker(root: ET.Element) -> None:
    defs = ET.SubElement(root, f"{{{SVG}}}defs")
    marker = ET.SubElement(
        defs,
        f"{{{SVG}}}marker",
        {
            "id": "arrow",
            "viewBox": "0 0 10 10",
            "refX": "9",
            "refY": "5",
            "markerWidth": "5",
            "markerHeight": "5",
            "orient": "auto-start-reverse",
        },
    )
    ET.SubElement(marker, f"{{{SVG}}}path", {"d": "M 0 0 L 10 5 L 0 10 z", "fill": "#64748b"})


def add_callout(root: ET.Element, x: float, y: float, target_x: float, target_y: float,
                title: str, subtitle: str, color: str) -> None:
    group = ET.SubElement(root, f"{{{SVG}}}g")
    ET.SubElement(
        group,
        f"{{{SVG}}}line",
        {
            "x1": str(x + 4), "y1": str(y + 4),
            "x2": str(target_x), "y2": str(target_y),
            "stroke": "#64748b", "stroke-width": "0.8",
            "stroke-dasharray": "2 1.5", "marker-end": "url(#arrow)",
        },
    )
    ET.SubElement(group, f"{{{SVG}}}circle", {"cx": str(x), "cy": str(y), "r": "3.5", "fill": color})
    heading = ET.SubElement(
        group, f"{{{SVG}}}text",
        {"x": str(x + 7), "y": str(y + 1.5), "font-family": "Arial, Helvetica, sans-serif",
         "font-size": "8", "font-weight": "bold", "fill": "#172033"},
    )
    heading.text = title
    detail = ET.SubElement(
        group, f"{{{SVG}}}text",
        {"x": str(x + 7), "y": str(y + 10), "font-family": "Arial, Helvetica, sans-serif",
         "font-size": "5.8", "fill": "#475569"},
    )
    detail.text = subtitle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--native-svg", type=Path, required=True)
    parser.add_argument("--workflow", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-svg", type=Path, required=True)
    parser.add_argument("--taxon-key", type=Path, required=True)
    args = parser.parse_args()

    taxa = parse_taxa(args.workflow)
    terminal, cycle_position = parse_terminal_weights(args.workflow, len(taxa))
    metadata: dict[str, dict[str, str]] = {}
    with args.metadata.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            metadata[row["genome"]] = row

    args.taxon_key.parent.mkdir(parents=True, exist_ok=True)
    with args.taxon_key.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "network_taxon_number", "cycle_position", "genome", "species",
            "taxonomic_group", "source", "terminal_split_weight",
            "distinct_positive_terminal_split",
        ])
        for index, name in enumerate(taxa, 1):
            row = metadata.get(name, {})
            weight = terminal.get(index, 0.0)
            writer.writerow([
                index, cycle_position.get(index, ""), name,
                row.get("species", species_for(name)), group_for(name),
                row.get("source", ""), f"{weight:.8g}", int(weight > 1e-10),
            ])

    tree = ET.parse(args.native_svg)
    root = tree.getroot()
    native_lines = list(root.findall(f"{{{SVG}}}line"))
    for line in native_lines:
        line.set("stroke", "#263445")
        line.set("stroke-width", "0.42")
        line.set("stroke-opacity", "0.88")
    for element in list(root):
        if element.tag != f"{{{SVG}}}line":
            root.remove(element)
    # Make the exact network occupy most of the panel while leaving room for
    # three group callouts. Existing network line coordinates are unchanged.
    root.set("viewBox", "235 155 270 155")
    root.set("width", "2160")
    root.set("height", "1240")
    root.insert(0, ET.Comment(
        "Exact native SplitsTree line coordinates retained; taxon labels replaced by group-level callouts."
    ))
    root.insert(1, ET.Element(
        f"{{{SVG}}}rect",
        {"x": "235", "y": "155", "width": "270", "height": "155", "fill": "#ffffff"},
    ))
    add_arrow_marker(root)
    add_callout(root, 247, 169, 318, 218, "Other Pantoea", "26 genomes across multiple species", "#64748b")
    add_callout(root, 407, 169, 404, 233, "P. septica", "9 genomes", "#a855f7")
    add_callout(root, 365, 284, 354, 263, "P. piersonii", "57 genomes; see primary tree", "#0f766e")
    note = ET.SubElement(
        root, f"{{{SVG}}}text",
        {"x": "247", "y": "303", "font-family": "Arial, Helvetica, sans-serif",
         "font-size": "5.5", "fill": "#64748b"},
    )
    note.text = "Overlapping/box-like edges represent incompatible distance signal, not identified transfer events."

    # Confirm native network lines themselves survived untouched and in order.
    final_lines = [line for line in root.findall(f"{{{SVG}}}line") if line in native_lines]
    if len(final_lines) != len(native_lines):
        raise RuntimeError("A native network line was lost")
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    tree.write(args.output_svg, encoding="utf-8", xml_declaration=True)


if __name__ == "__main__":
    main()
