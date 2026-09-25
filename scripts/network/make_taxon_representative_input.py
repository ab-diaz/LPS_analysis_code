#!/usr/bin/env python3
"""Create a readable taxon-representative distance matrix from the 92-genome matrix.

For taxa represented by multiple genomes, select the medoid: the genome with the
smallest mean JTT distance to the other genomes assigned to the same taxon.
Singleton taxa are retained unchanged.
"""

from __future__ import annotations

import csv
import re
import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np


def parse_distances(path: Path) -> tuple[list[str], np.ndarray]:
    text = path.read_text()
    block = re.search(r"BEGIN\s+DISTANCES\s*;(.*?)END\s*;", text, re.I | re.S)
    if not block:
        raise ValueError("No DISTANCES block found")
    matrix = re.search(r"\bMATRIX\b(.*?);", block.group(1), re.I | re.S)
    if not matrix:
        raise ValueError("No MATRIX found")
    row_re = re.compile(r"^\[(\d+)\]\s+'((?:[^']|'')*)'\s+(.*)$")
    names, rows = [], []
    for raw in matrix.group(1).splitlines():
        match = row_re.match(raw.strip())
        if match:
            names.append(match.group(2).replace("''", "'"))
            rows.append([float(value) for value in match.group(3).split()])
    values = np.asarray(rows)
    if values.shape != (len(names), len(names)):
        raise ValueError(f"Unexpected distance-matrix shape: {values.shape}")
    return names, values


def quote(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def main() -> None:
    parser = argparse.ArgumentParser(description="Select within-taxon JTT-distance medoids")
    parser.add_argument("--distances", type=Path, required=True)
    parser.add_argument("--key", type=Path, required=True)
    parser.add_argument("--out-nexus", type=Path, required=True)
    parser.add_argument("--out-table", type=Path, required=True)
    args = parser.parse_args()
    distances_path = args.distances
    key_path = args.key
    out_nexus = args.out_nexus
    out_table = args.out_table
    names, distances = parse_distances(distances_path)
    name_to_index = {name: index for index, name in enumerate(names)}
    with key_path.open() as handle:
        metadata = {row["genome"]: row for row in csv.DictReader(handle, delimiter="\t")}
    missing = sorted(set(names) - set(metadata))
    if missing:
        raise ValueError(f"Taxa absent from metadata: {missing}")

    by_species: dict[str, list[int]] = defaultdict(list)
    for index, name in enumerate(names):
        by_species[metadata[name]["species"]].append(index)

    chosen: list[tuple[str, int, float]] = []
    for species, indices in sorted(by_species.items()):
        if len(indices) == 1:
            selected, mean_distance = indices[0], 0.0
        else:
            within = distances[np.ix_(indices, indices)]
            means = within.sum(axis=1) / (len(indices) - 1)
            local = min(range(len(indices)), key=lambda i: (means[i], names[indices[i]]))
            selected, mean_distance = indices[local], float(means[local])
        chosen.append((species, selected, mean_distance))

    # Alphabetic species order makes the file deterministic; NeighborNet determines
    # its own circular order and is unaffected by this input ordering.
    selected_indices = [index for _, index, _ in chosen]
    selected_names = [names[index] for index in selected_indices]
    subset = distances[np.ix_(selected_indices, selected_indices)]

    out_nexus.parent.mkdir(parents=True, exist_ok=True)
    with out_nexus.open("w") as out:
        out.write("#NEXUS\n\nBEGIN DISTANCES;\n")
        out.write(f"DIMENSIONS ntax={len(selected_names)};\n")
        out.write("FORMAT labels=left diagonal triangle=Both;\nMATRIX\n")
        for row, (name, values) in enumerate(zip(selected_names, subset), 1):
            payload = " ".join(f"{value:.8f}" for value in values)
            out.write(f"[{row}] {quote(name)} {payload}\n")
        out.write(";\nEND;\n")

    out_table.parent.mkdir(parents=True, exist_ok=True)
    with out_table.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "species", "representative_genome", "genomes_represented",
            "selection_rule", "mean_within_taxon_JTT_distance"
        ])
        for species, selected, mean_distance in chosen:
            count = len(by_species[species])
            rule = "singleton" if count == 1 else "within-taxon medoid"
            writer.writerow([species, names[selected], count, rule, f"{mean_distance:.8g}"])

    print(f"Wrote {out_nexus} ({len(selected_names)} representatives)")
    print(f"Wrote {out_table}")


if __name__ == "__main__":
    main()
