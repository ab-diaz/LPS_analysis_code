#!/usr/bin/env python3
import argparse
import csv
import re


def main():
    parser = argparse.ArgumentParser(description="Add HyPhy RELAX branch labels to terminal branches by ISS/Earth group.")
    parser.add_argument("--tree", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    groups = {}
    with open(args.map, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["gene"] == args.gene:
                groups[row["safe_id"]] = row["group"]

    with open(args.tree) as handle:
        tree = handle.read().strip()

    for leaf, group in sorted(groups.items(), key=lambda item: len(item[0]), reverse=True):
        tree = re.sub(rf"(?<![A-Za-z0-9_.-])({re.escape(leaf)})(?=[:),])", rf"\1{{{group}}}", tree)

    with open(args.out, "w") as handle:
        handle.write(tree + "\n")


if __name__ == "__main__":
    main()
