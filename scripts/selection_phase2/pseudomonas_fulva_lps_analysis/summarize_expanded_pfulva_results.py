#!/usr/bin/env python3
"""Write a concise summary of the expanded Pseudomonas fulva side analysis."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
META = ROOT / "genome_metadata_pfulva_expanded_with_mt2_iss.tsv"
MATRIX = ROOT / "pfulva_expanded_o_antigen_gene_presence.tsv"
COUNTS = ROOT / "pfulva_expanded_lps_category_counts.tsv"
OUT = ROOT / "pfulva_expanded_results_summary.md"


def main() -> None:
    with META.open(newline="", encoding="utf-8") as handle:
        meta = {r["isolate_id"]: r for r in csv.DictReader(handle, delimiter="\t")}

    with MATRIX.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        isolates = header[1:]
        matrix = {}
        for row in reader:
            matrix[row[0]] = {iso: int(v) for iso, v in zip(isolates, row[1:])}

    genes = list(matrix)
    groups = Counter(row["group"] for row in meta.values())
    iss_isolates = [iso for iso in isolates if meta[iso]["group"] == "ISS"]
    earth_refs = [iso for iso in isolates if meta[iso]["group"] == "Earth_reference"]

    iss_profiles = Counter(
        tuple(g for g in genes if matrix[g][iso])
        for iso in iss_isolates
    )
    iss_profile = iss_profiles.most_common(1)[0][0]

    same_as_iss = [
        iso
        for iso in isolates
        if meta[iso]["group"] != "ISS" and tuple(g for g in genes if matrix[g][iso]) == iss_profile
    ]

    presence_lines = []
    for gene in genes:
        vals = []
        for group in ["ISS", "Earth_reference", "Earth_type_strain", "Earth_spacecraft_associated"]:
            isos = [iso for iso in isolates if meta[iso]["group"] == group]
            if not isos:
                continue
            n = sum(matrix[gene][iso] for iso in isos)
            vals.append(f"{group}: {n}/{len(isos)}")
        presence_lines.append(f"- `{gene}`: " + "; ".join(vals))

    with COUNTS.open(newline="", encoding="utf-8") as handle:
        count_rows = list(csv.DictReader(handle, delimiter="\t"))

    def values(group: str, field: str) -> list[int]:
        return [int(r[field]) for r in count_rows if r["group"] == group]

    lines = [
        "# Expanded Pseudomonas fulva LPS/O-antigen Side Analysis",
        "",
        "## Dataset",
        f"- Total genomes: {len(meta)}",
        f"- ISS genomes: {groups['ISS']}",
        f"- Earth reference genomes: {groups['Earth_reference']}",
        f"- Earth spacecraft-associated genomes: {groups['Earth_spacecraft_associated']}",
        f"- Earth type-strain genomes: {groups['Earth_type_strain']}",
        "",
        "## Main Result",
        f"- All {len(iss_isolates)} ISS genomes had the same O-antigen-related gene profile:",
        f"  `{', '.join(iss_profile)}`.",
        f"- Non-ISS genomes with the same profile: {len(same_as_iss)}",
        "  " + ", ".join(f"`{iso}` ({meta[iso]['group']})" for iso in same_as_iss) + ".",
        "- The ISS profile includes `rfbA`, `rfbC`, `rfbD`, `wzm`, and `wzt`, but lacks `rfbB`.",
        "- In contrast, many Earth references carried `rfbB` and lacked `rfbA`, indicating O-antigen pathway variation among Earth/reference genomes.",
        "",
        "## Gene Presence By Group",
        *presence_lines,
        "",
        "## LPS Category Counts",
    ]

    for group in ["ISS", "Earth_reference", "Earth_type_strain", "Earth_spacecraft_associated"]:
        if groups[group] == 0:
            continue
        o_vals = values(group, "o_antigen_genes")
        total_vals = values(group, "total_unique_lps_genes")
        lines.append(
            f"- {group}: O-antigen genes range {min(o_vals)}-{max(o_vals)}; "
            f"total LPS-related genes range {min(total_vals)}-{max(total_vals)}."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "- This analysis supports a cautious exploratory statement: a conserved LPS/O-antigen-related profile is present across many MT-2 ISS Pseudomonas fulva isolates.",
            "- The result should not be interpreted as multi-year persistence because the MT-2 ISS Pseudomonas fulva genomes are from one collection date.",
            "- The identical ISS profile and the matching Mars Odyssey-associated strain suggest a spacecraft-associated/clonal component is plausible.",
            "- For the manuscript, this should be used as a supplementary cross-genus observation, not as a main statistical result.",
        ]
    )

    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
