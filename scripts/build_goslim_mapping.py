#!/usr/bin/env python3
import argparse
import sys


def parse_obo_terms(path):
    terms = {}
    current = None
    with open(path, "r", newline="") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line == "[Term]":
                if current and "id" in current:
                    terms[current["id"]] = current
                current = {"is_a": [], "part_of": [], "obsolete": False}
                continue
            if line.startswith("["):
                if current and "id" in current:
                    terms[current["id"]] = current
                current = None
                continue
            if current is None:
                continue
            if line.startswith("id:"):
                current["id"] = line.split("id:", 1)[1].strip()
            elif line.startswith("name:"):
                current["name"] = line.split("name:", 1)[1].strip()
            elif line.startswith("is_a:"):
                parent = line.split("is_a:", 1)[1].split("!")[0].strip()
                if parent:
                    current["is_a"].append(parent)
            elif line.startswith("relationship:"):
                parts = line.split()
                if len(parts) >= 3 and parts[1] == "part_of":
                    parent = parts[2].strip()
                    if parent:
                        current["part_of"].append(parent)
            elif line.startswith("is_obsolete:"):
                value = line.split("is_obsolete:", 1)[1].strip().lower()
                current["obsolete"] = value == "true"
    if current and "id" in current:
        terms[current["id"]] = current
    return terms


def build_parent_map(terms):
    parents = {}
    for term_id, data in terms.items():
        rels = []
        rels.extend(data.get("is_a", []))
        rels.extend(data.get("part_of", []))
        parents[term_id] = rels
    return parents


def compute_slim_ancestors(term_id, parents, slim_ids, memo, visiting):
    if term_id in memo:
        return memo[term_id]
    if term_id in visiting:
        return set()
    visiting.add(term_id)
    slim_hits = set()
    if term_id in slim_ids:
        slim_hits.add(term_id)
    for parent in parents.get(term_id, []):
        slim_hits.update(compute_slim_ancestors(parent, parents, slim_ids, memo, visiting))
    visiting.remove(term_id)
    memo[term_id] = slim_hits
    return slim_hits


def main():
    parser = argparse.ArgumentParser(
        description="Generate GO term to GO-slim mapping TSV using OBO files."
    )
    parser.add_argument("--go-obo", required=True, help="Path to go-basic.obo")
    parser.add_argument("--slim-obo", required=True, help="Path to GO-slim OBO (prokaryote)")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    go_terms = parse_obo_terms(args.go_obo)
    slim_terms = parse_obo_terms(args.slim_obo)
    slim_ids = {term_id for term_id, data in slim_terms.items() if not data.get("obsolete")}
    slim_names = {
        term_id: data.get("name", term_id) for term_id, data in slim_terms.items()
    }

    parents = build_parent_map(go_terms)
    memo = {}

    with open(args.out, "w", newline="") as out:
        for term_id, data in go_terms.items():
            if data.get("obsolete"):
                continue
            slim_anc = compute_slim_ancestors(term_id, parents, slim_ids, memo, set())
            if not slim_anc:
                continue
            for slim_id in sorted(slim_anc):
                slim_label = f"{slim_id} {slim_names.get(slim_id, slim_id)}"
                out.write(f"{term_id}\t{slim_label}\n")


if __name__ == "__main__":
    main()
