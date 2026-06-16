#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

from Bio import Phylo, SeqIO


GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL"]


def load_groups(map_path, gene):
    groups = {}
    rows = {}
    with open(map_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["gene"] != gene:
                continue
            groups[row["safe_id"]] = row["group"]
            rows[row["safe_id"]] = row
    return groups, rows


def strip_labels(tree):
    for clade in tree.find_clades():
        if clade.name:
            clade.name = clade.name.replace("{ISS}", "").replace("{Earth}", "")
        clade.comment = None
    return tree


def label_tree(tree, groups, mode, single_test=None):
    tree = strip_labels(tree)

    for tip in tree.get_terminals():
        group = groups.get(tip.name)
        if mode == "single_test":
            if tip.name == single_test:
                tip.name = f"{tip.name}{{ISS}}"
            elif group == "Earth":
                tip.name = f"{tip.name}{{Earth}}"
            elif group == "ISS":
                # Other ISS branches are left unclassified in this diagnostic test.
                pass
        elif group in {"ISS", "Earth"}:
            tip.name = f"{tip.name}{{{group}}}"

    if mode == "internal_all_descendants":
        for clade in tree.get_nonterminals(order="postorder"):
            desc_groups = {groups.get(t.name) for t in clade.get_terminals()}
            if desc_groups == {"ISS"}:
                clade.name = "{ISS}"
            elif desc_groups == {"Earth"}:
                clade.name = "{Earth}"

    return tree


def terminal_branch_lengths(tree, groups):
    rows = []
    for tip in tree.get_terminals():
        rows.append(
            {
                "safe_id": tip.name,
                "group": groups.get(tip.name, ""),
                "terminal_branch_length": tip.branch_length if tip.branch_length is not None else 0.0,
            }
        )
    return rows


def choose_divergent_tip(tree, groups):
    iss = [r for r in terminal_branch_lengths(tree, groups) if r["group"] == "ISS"]
    if not iss:
        return None
    return max(iss, key=lambda r: (r["terminal_branch_length"], r["safe_id"]))["safe_id"]


def write_pruned_alignment(input_fasta, output_fasta, remove_id):
    kept = []
    for rec in SeqIO.parse(input_fasta, "fasta"):
        if rec.id == remove_id:
            continue
        kept.append(rec)
    SeqIO.write(kept, output_fasta, "fasta")
    return len(kept)


def main():
    parser = argparse.ArgumentParser(description="Prepare RELAX sensitivity inputs.")
    parser.add_argument("--root", default=".")
    parser.add_argument("--outdir", default="relax_sensitivity_20260530")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    outdir = root / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    map_path = root / "sequence_id_map.tsv"

    divergent_rows = []
    branch_rows = []
    run_rows = []

    for gene in GENES:
        gene_dir = outdir / gene
        gene_dir.mkdir(exist_ok=True)
        groups, meta = load_groups(map_path, gene)
        tree_path = root / f"{gene}.codon.aln.fasta.treefile"
        aln_path = root / f"{gene}.codon.aln.fasta"
        tree = Phylo.read(tree_path, "newick")
        tree = strip_labels(tree)

        for row in terminal_branch_lengths(tree, groups):
            row["gene"] = gene
            branch_rows.append(row)

        divergent = choose_divergent_tip(tree, groups)
        divergent_rows.append(
            {
                "gene": gene,
                "divergent_safe_id": divergent or "",
                "group": groups.get(divergent, "") if divergent else "",
                "isolate": meta.get(divergent, {}).get("isolate", "") if divergent else "",
                "locus": meta.get(divergent, {}).get("locus", "") if divergent else "",
            }
        )

        internal_tree = label_tree(Phylo.read(tree_path, "newick"), groups, "internal_all_descendants")
        internal_tree_path = gene_dir / f"{gene}.relax_internal_all_descendants.tree"
        Phylo.write(internal_tree, internal_tree_path, "newick")
        run_rows.append(
            {
                "gene": gene,
                "analysis": "internal_all_descendants",
                "alignment": str(aln_path),
                "tree": str(internal_tree_path),
                "output": str(gene_dir / f"{gene}.RELAX.internal_all_descendants.json"),
                "test": "ISS",
                "reference": "Earth",
                "divergent_safe_id": "",
            }
        )

        if divergent:
            single_tree = label_tree(Phylo.read(tree_path, "newick"), groups, "single_test", single_test=divergent)
            single_tree_path = gene_dir / f"{gene}.relax_single_divergent.tree"
            Phylo.write(single_tree, single_tree_path, "newick")
            run_rows.append(
                {
                    "gene": gene,
                    "analysis": "single_divergent_iss_tip",
                    "alignment": str(aln_path),
                    "tree": str(single_tree_path),
                    "output": str(gene_dir / f"{gene}.RELAX.single_divergent_iss_tip.json"),
                    "test": "ISS",
                    "reference": "Earth",
                    "divergent_safe_id": divergent,
                }
            )

            pruned_aln = gene_dir / f"{gene}.codon.aln.minus_divergent.fasta"
            n_kept = write_pruned_alignment(aln_path, pruned_aln, divergent)
            pruned_tree = Phylo.read(tree_path, "newick")
            pruned_tree = strip_labels(pruned_tree)
            pruned_tree.prune(divergent)
            pruned_labeled = label_tree(pruned_tree, groups, "internal_all_descendants")
            pruned_tree_path = gene_dir / f"{gene}.relax_minus_divergent_internal.tree"
            Phylo.write(pruned_labeled, pruned_tree_path, "newick")
            run_rows.append(
                {
                    "gene": gene,
                    "analysis": "minus_divergent_internal_all_descendants",
                    "alignment": str(pruned_aln),
                    "tree": str(pruned_tree_path),
                    "output": str(gene_dir / f"{gene}.RELAX.minus_divergent_internal_all_descendants.json"),
                    "test": "ISS",
                    "reference": "Earth",
                    "divergent_safe_id": divergent,
                    "n_sequences": str(n_kept),
                }
            )

    with open(outdir / "terminal_branch_lengths.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["gene", "safe_id", "group", "terminal_branch_length"])
        writer.writeheader()
        writer.writerows(branch_rows)

    with open(outdir / "divergent_iss_tips.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["gene", "divergent_safe_id", "group", "isolate", "locus"])
        writer.writeheader()
        writer.writerows(divergent_rows)

    all_fields = sorted({k for row in run_rows for k in row})
    preferred = ["gene", "analysis", "alignment", "tree", "output", "test", "reference", "divergent_safe_id", "n_sequences"]
    fields = preferred + [f for f in all_fields if f not in preferred]
    with open(outdir / "relax_sensitivity_runs.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(run_rows)

    print(outdir)


if __name__ == "__main__":
    main()
