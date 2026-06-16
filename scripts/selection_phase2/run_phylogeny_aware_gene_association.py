#!/usr/bin/env python3
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.stats import fisher_exact


BASE = Path(__file__).resolve().parent
MATRIX = BASE / "phylogeny_aware_inputs" / "focal_gene_matrix.tsv"
TREE = (
    BASE
    / "phylogenomic_context_phylophlan"
    / "output"
    / "pantoea_92_phylophlan_marker_tree"
    / "input_proteomes.tre.treefile"
)
OUTDIR = BASE / "phylogeny_aware_association"
GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL", "wzm", "wzt", "rfaZ"]
N_PERMUTATIONS = 10000
RANDOM_SEED = 20260531


def bh_qvalues(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def association_stat(gene_values, labels):
    iss = labels == "ISS"
    earth = labels == "Earth"
    iss_prev = gene_values[iss].mean()
    earth_prev = gene_values[earth].mean()
    return iss_prev - earth_prev


def permute_within_groups(labels, groups, rng):
    out = labels.copy()
    for group in np.unique(groups):
        idx = np.where(groups == group)[0]
        if len(idx) > 1:
            out[idx] = rng.permutation(out[idx])
    return out


def empirical_p(observed, null_values):
    # Two-sided empirical p-value with +1 correction.
    return (np.sum(np.abs(null_values) >= abs(observed)) + 1) / (len(null_values) + 1)


def directional_p(observed, null_values):
    # Direction follows the observed sign; useful for interpretation.
    if observed >= 0:
        count = np.sum(null_values >= observed)
    else:
        count = np.sum(null_values <= observed)
    return (count + 1) / (len(null_values) + 1)


def main():
    OUTDIR.mkdir(exist_ok=True)
    df = pd.read_csv(MATRIX, sep="\t")
    tree = Phylo.read(TREE, "newick")
    tips = [tip.name for tip in tree.get_terminals()]

    missing_in_tree = sorted(set(df["genome"]) - set(tips))
    extra_in_tree = sorted(set(tips) - set(df["genome"]))
    if missing_in_tree or extra_in_tree:
        raise SystemExit(
            "Tree/matrix mismatch:\n"
            f"missing_in_tree={missing_in_tree[:10]}\n"
            f"extra_in_tree={extra_in_tree[:10]}"
        )

    df = df.set_index("genome").loc[tips].reset_index()
    labels = df["source"].to_numpy()
    species = df["species"].to_numpy()
    rng = np.random.default_rng(RANDOM_SEED)

    global_label_perms = []
    species_label_perms = []
    for _ in range(N_PERMUTATIONS):
        global_label_perms.append(rng.permutation(labels))
        species_label_perms.append(permute_within_groups(labels, species, rng))
    global_label_perms = np.asarray(global_label_perms)
    species_label_perms = np.asarray(species_label_perms)

    rows = []
    null_rows = []
    for gene in GENES:
        values = df[gene].to_numpy(dtype=int)
        observed = association_stat(values, labels)

        iss_present = int(values[labels == "ISS"].sum())
        iss_absent = int((labels == "ISS").sum() - iss_present)
        earth_present = int(values[labels == "Earth"].sum())
        earth_absent = int((labels == "Earth").sum() - earth_present)
        odds_ratio, fisher_p = fisher_exact(
            [[iss_present, iss_absent], [earth_present, earth_absent]],
            alternative="two-sided",
        )

        null_global = np.array(
            [association_stat(values, perm_labels) for perm_labels in global_label_perms]
        )
        null_species = np.array(
            [association_stat(values, perm_labels) for perm_labels in species_label_perms]
        )

        for null_name, null_values in [
            ("global", null_global),
            ("species_restricted", null_species),
        ]:
            null_rows.extend(
                {
                    "gene": gene,
                    "null_model": null_name,
                    "permutation": i,
                    "statistic": stat,
                }
                for i, stat in enumerate(null_values)
            )

        rows.append(
            {
                "gene": gene,
                "iss_present": iss_present,
                "iss_total": int((labels == "ISS").sum()),
                "earth_present": earth_present,
                "earth_total": int((labels == "Earth").sum()),
                "iss_prevalence": iss_present / int((labels == "ISS").sum()),
                "earth_prevalence": earth_present / int((labels == "Earth").sum()),
                "observed_prevalence_difference": observed,
                "fisher_odds_ratio": odds_ratio,
                "fisher_p": fisher_p,
                "global_empirical_p": empirical_p(observed, null_global),
                "global_directional_p": directional_p(observed, null_global),
                "species_restricted_empirical_p": empirical_p(observed, null_species),
                "species_restricted_directional_p": directional_p(observed, null_species),
            }
        )

    results = pd.DataFrame(rows)
    for col in [
        "fisher_p",
        "global_empirical_p",
        "global_directional_p",
        "species_restricted_empirical_p",
        "species_restricted_directional_p",
    ]:
        results[col.replace("_p", "_q")] = bh_qvalues(results[col])

    def interpret(row):
        direction = "ISS-enriched" if row["observed_prevalence_difference"] > 0 else "Earth-enriched"
        if row["species_restricted_empirical_q"] <= 0.05:
            return f"{direction}; remains unusual after species-restricted permutation"
        if row["global_empirical_q"] <= 0.05:
            return f"{direction}; significant globally but not after species restriction"
        return f"{direction}; not significant after permutation"

    results["interpretation"] = results.apply(interpret, axis=1)
    results.to_csv(OUTDIR / "phylogeny_aware_gene_association.tsv", sep="\t", index=False)
    pd.DataFrame(null_rows).to_csv(
        OUTDIR / "phylogeny_aware_gene_association_null_distributions.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    summary_table = results[
        [
            "gene",
            "iss_present",
            "iss_total",
            "earth_present",
            "earth_total",
            "observed_prevalence_difference",
            "global_empirical_q",
            "species_restricted_empirical_q",
            "interpretation",
        ]
    ].copy()
    summary_table = summary_table.round(
        {
            "observed_prevalence_difference": 4,
            "global_empirical_q": 4,
            "species_restricted_empirical_q": 4,
        }
    )
    markdown_lines = [
        "| " + " | ".join(summary_table.columns) + " |",
        "| " + " | ".join(["---"] * len(summary_table.columns)) + " |",
    ]
    for _, row in summary_table.iterrows():
        markdown_lines.append("| " + " | ".join(str(row[col]) for col in summary_table.columns) + " |")

    summary = [
        "# Phylogeny-aware focal gene association",
        "",
        f"- Matrix: `{MATRIX.relative_to(BASE)}`",
        f"- Tree: `{TREE.relative_to(BASE)}`",
        f"- Genomes: {len(df)}",
        f"- ISS genomes: {(labels == 'ISS').sum()}",
        f"- Earth genomes: {(labels == 'Earth').sum()}",
        f"- Permutations per null model: {N_PERMUTATIONS}",
        f"- Random seed: {RANDOM_SEED}",
        "",
        "The test statistic is ISS prevalence minus Earth prevalence.",
        "Two null models were used: unrestricted global label shuffling, and species-restricted shuffling that preserves the ISS/Earth composition within each species.",
        "",
        "## Results",
        "",
        "\n".join(markdown_lines),
        "",
    ]
    (OUTDIR / "phylogeny_aware_gene_association_summary.md").write_text(
        "\n".join(summary) + "\n"
    )

    print(results.to_string(index=False))
    print(f"\nWrote: {OUTDIR / 'phylogeny_aware_gene_association.tsv'}")


if __name__ == "__main__":
    main()
