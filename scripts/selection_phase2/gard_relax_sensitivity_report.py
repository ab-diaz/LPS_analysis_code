#!/usr/bin/env python3
import csv
import json
from pathlib import Path


ROOT = Path("results/lps_analysis/selection_phase2")
OUT = ROOT / "relax_sensitivity_20260530"


def gard_summary():
    path = OUT / "gard_rfbC_default" / "rfbC.GARD.default_clean.json"
    d = json.load(open(path))
    return {
        "gene": "rfbC",
        "run": str(path),
        "sequences": d["input"].get("number of sequences", ""),
        "codon_sites": d["input"].get("number of sites", ""),
        "partition_count": d["input"].get("partition count", ""),
        "best_model_aicc": d.get("bestModelAICc", ""),
        "single_tree_aicc": d.get("singleTreeAICc", ""),
        "improvements": json.dumps(d.get("improvements", {}), sort_keys=True),
        "interpretation": "no retained recombination breakpoint",
        "caveat": "JSON completed; HyPhy process did not return cleanly after writing result",
    }


def main():
    gard = gard_summary()
    with open(OUT / "GARD_rfbC_clean_summary.tsv", "w", newline="") as handle:
        fields = list(gard)
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerow(gard)

    relax_rows = list(csv.DictReader(open(OUT / "RELAX_sensitivity_summary.tsv"), delimiter="\t"))
    key = []
    for gene in ["rfbA", "rfbB", "rfbC", "rfbD", "waaL"]:
        subset = [r for r in relax_rows if r["gene"] == gene]
        original = next(r for r in subset if r["analysis"] == "original_terminal_only")
        minus = next(r for r in subset if r["analysis"] == "minus_divergent_internal_all_descendants")
        single = next(r for r in subset if r["analysis"] == "single_divergent_iss_tip")
        internal = next(r for r in subset if r["analysis"] == "internal_all_descendants")
        key.append(
            {
                "gene": gene,
                "original_k": original["k"],
                "original_p": original["p_value"],
                "internal_k": internal["k"],
                "internal_p": internal["p_value"],
                "minus_divergent_k": minus["k"],
                "minus_divergent_p": minus["p_value"],
                "single_divergent_k": single["k"],
                "single_divergent_p": single["p_value"],
                "interpretation": (
                    "original_not_significant"
                    if float(original["p_value"]) >= 0.05
                    else (
                        "robust_to_minus_divergent"
                        if float(minus["p_value"]) < 0.05
                        else "original_signal_lost_after_minus_divergent"
                    )
                ),
            }
        )
    with open(OUT / "RELAX_sensitivity_key_comparison.tsv", "w", newline="") as handle:
        fields = list(key[0])
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(key)

    print(OUT / "GARD_rfbC_clean_summary.tsv")
    print(OUT / "RELAX_sensitivity_key_comparison.tsv")


if __name__ == "__main__":
    main()
