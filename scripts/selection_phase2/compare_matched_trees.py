#!/usr/bin/env python3
"""Check the focal ISS-rich split and unrooted topology in matched trees."""

import csv
import json
import sys
from pathlib import Path

from Bio import Phylo


manifest, before_path, after_path, output_path = map(Path, sys.argv[1:5])
with manifest.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
taxa = {row["genome"] for row in rows}
iss = {row["genome"] for row in rows if row["source"] == "ISS"}
earth = {row["genome"] for row in rows if row["source"] == "Earth"}
focal = iss | {"Pantoea_piersonii_Y3_I4"}
if len(rows) != 57 or len(taxa) != 57 or len(iss) != 40 or len(earth) != 17:
    raise ValueError("Unexpected manifest size or source counts")


def read_tree(path):
    tree = Phylo.read(path, "newick")
    names = [tip.name for tip in tree.get_terminals()]
    if len(names) != 57 or len(set(names)) != 57 or set(names) != taxa:
        raise ValueError(f"Tree labels do not match the manifest: {path}")
    return tree


def tree_splits(tree):
    splits = {}
    for clade in tree.get_nonterminals():
        desc = frozenset(tip.name for tip in clade.get_terminals())
        other = frozenset(taxa - desc)
        if len(desc) < 2 or len(other) < 2:
            continue
        key = desc if (len(desc), sorted(desc)) <= (len(other), sorted(other)) else other
        splits[key] = {"support": clade.name, "branch_length": clade.branch_length}
    return splits


def focal_support(splits):
    key = frozenset(focal)
    other = frozenset(taxa - focal)
    match = splits.get(key if (len(key), sorted(key)) <= (len(other), sorted(other)) else other)
    return match


def strong_splits(splits):
    result = set()
    for key, value in splits.items():
        try:
            sh_alrt, ufboot = (float(part) for part in value["support"].split("/"))
        except (AttributeError, TypeError, ValueError):
            continue
        if sh_alrt >= 80 and ufboot >= 95:
            result.add(key)
    return result


before = tree_splits(read_tree(before_path))
after = tree_splits(read_tree(after_path))
strong_before = strong_splits(before)
strong_after = strong_splits(after)
result = {
    "taxa": len(taxa),
    "iss": len(iss),
    "earth": len(earth),
    "focal_split_size": len(focal),
    "focal_earth_taxon": "Pantoea_piersonii_Y3_I4",
    "focal_split_before": focal_support(before),
    "focal_split_after": focal_support(after),
    "unrooted_nontrivial_splits_before": len(before),
    "unrooted_nontrivial_splits_after": len(after),
    "shared_splits": len(set(before) & set(after)),
    "unrooted_robinson_foulds_distance": len(set(before) ^ set(after)),
    "strong_splits_before": len(strong_before),
    "strong_splits_after": len(strong_after),
    "strong_splits_shared": len(strong_before & strong_after),
    "strong_splits_only_before": len(strong_before - set(after)),
    "strong_splits_only_after": len(strong_after - set(before)),
    "strong_discordant_details": [
        {"tree": label, "smaller_split_side": sorted(key), "support": splits[key]["support"],
         "branch_length": splits[key]["branch_length"]}
        for label, splits, keys in (
            ("before", before, strong_before - set(after)),
            ("after", after, strong_after - set(before)),
        )
        for key in sorted(keys, key=lambda item: (len(item), sorted(item)))
    ],
}
output_path.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
