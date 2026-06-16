#!/usr/bin/env python3
import csv
import json
from pathlib import Path


GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL"]


def parse_relax(path):
    d = json.load(open(path))
    tr = d["test results"]
    return {
        "k": tr.get("relaxation or intensification parameter", ""),
        "LRT": tr.get("LRT", ""),
        "p_value": tr.get("p-value", ""),
        "sequences": d["input"].get("number of sequences", ""),
        "sites": d["input"].get("number of sites", ""),
    }


def main():
    root = Path("results/lps_analysis/selection_phase2")
    outdir = root / "relax_sensitivity_20260530"
    rows = []

    current = {
        "rfbA": root / "rfbA.codon.aln.fasta.RELAX.json",
        "rfbB": root / "rfbB.codon.aln.fasta.RELAX.json",
        "rfbC": root / "rfbC.codon.aln.fasta.RELAX.json",
        "rfbD": root / "rfbD.codon.aln.fasta.RELAX.json",
        "waaL": root / "waaL.codon.aln.fasta.RELAX.json",
    }
    for gene, path in current.items():
        rec = parse_relax(path)
        rec.update({"gene": gene, "analysis": "original_terminal_only", "divergent_safe_id": ""})
        rows.append(rec)

    run_rows = list(csv.DictReader(open(outdir / "relax_sensitivity_runs.tsv"), delimiter="\t"))
    for run in run_rows:
        path = Path(run["output"])
        rec = parse_relax(path)
        rec.update(
            {
                "gene": run["gene"],
                "analysis": run["analysis"],
                "divergent_safe_id": run.get("divergent_safe_id", ""),
            }
        )
        rows.append(rec)

    for row in rows:
        try:
            k = float(row["k"])
            p = float(row["p_value"])
            row["direction"] = "intensified" if k > 1 else "relaxed" if k < 1 else "unchanged"
            row["significant_0.05"] = "1" if p < 0.05 else "0"
        except Exception:
            row["direction"] = ""
            row["significant_0.05"] = ""

    fields = [
        "gene",
        "analysis",
        "sequences",
        "sites",
        "k",
        "LRT",
        "p_value",
        "direction",
        "significant_0.05",
        "divergent_safe_id",
    ]
    out = outdir / "RELAX_sensitivity_summary.tsv"
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(out)


if __name__ == "__main__":
    main()
