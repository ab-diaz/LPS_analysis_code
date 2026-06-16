#!/usr/bin/env python3
from pathlib import Path
import os
import shutil

import pandas as pd


OUT = Path("supplementary_tables_final/Supplementary_Table_4_association_founder_phylogeny.tsv")
BACKUP = OUT.with_suffix(".before_publication_summary.tsv")

KEY = Path("pantoea_route1_strengthening/route1_key_results.tsv")
YEAR = Path("pantoea_route1_strengthening/iss_year_persistence_target_genes.tsv")
PHYLO = Path("phylogeny_aware_association/phylogeny_aware_gene_association.tsv")


def fmt_float(x, digits=3):
    try:
        val = float(x)
    except Exception:
        return str(x)
    if abs(val) < 0.001 and val != 0:
        return f"{val:.2e}"
    return f"{val:.{digits}g}"


def year_summary(year_df: pd.DataFrame, gene: str) -> str:
    sub = year_df[year_df["gene_id"] == gene].sort_values("year")
    parts = []
    for row in sub.itertuples(index=False):
        parts.append(f"{int(row.year)}:{int(row.present)}/{int(row.iss_isolates)}")
    return "; ".join(parts)


def direction_from_diff(diff: float) -> str:
    if diff > 0:
        return "ISS-enriched"
    if diff < 0:
        return "Earth-enriched"
    return "No directional difference"


def overall_interpretation(row) -> str:
    gene = row["gene"]
    direction = row["direction"]
    if gene in {"rfbA", "rfbB", "rfbC", "rfbD"}:
        return (
            "Part of the ISS-retained rfb/rhamnose-biosynthesis module; signal is strong at isolate level "
            "and persists across ISS years, but cluster-collapsed support is weaker because ISS isolates are closely related."
        )
    if gene == "waaL":
        return (
            "ISS-retained O-antigen ligase; strongest founder-aware support among focal ISS-enriched genes "
            "and present in all ISS isolates."
        )
    if gene == "wzm":
        return (
            "Frequently retained in ISS isolates and supported by phylogeny-aware testing, but less universally present than rfbABCD/waaL."
        )
    if gene in {"wzt", "rfaZ"} and direction == "Earth-enriched":
        return (
            "More frequent in Earth-derived genomes and absent from ISS sampling years in this dataset; useful as a contrasting Earth-associated marker."
        )
    return row["phylogeny_aware_interpretation"]


def rebuild_table4() -> None:
    key = pd.read_csv(KEY, sep="\t")
    year = pd.read_csv(YEAR, sep="\t")
    phylo = pd.read_csv(PHYLO, sep="\t")
    merged = key.merge(phylo, on="gene", how="left", suffixes=("_key", "_phylo"))

    rows = []
    for _, r in merged.iterrows():
        diff = float(r["observed_prevalence_difference"])
        iss_clusters_present = int(r["strict_profile_cluster_iss_present"])
        earth_clusters_present = int(r["strict_profile_cluster_earth_present"])
        iss_clusters_total = 5
        earth_clusters_total = 30
        out = {
            "supplementary_table": "Supplementary Table 4",
            "gene": r["gene"],
            "direction": direction_from_diff(diff),
            "iss_present_isolates": f"{int(r['iss_present'])}/{int(r['iss_total'])}",
            "earth_present_genomes": f"{int(r['earth_present'])}/{int(r['earth_total'])}",
            "iss_prevalence": fmt_float(r["iss_prevalence"]),
            "earth_prevalence": fmt_float(r["earth_prevalence"]),
            "prevalence_difference_iss_minus_earth": fmt_float(diff),
            "isolate_level_fisher_q": fmt_float(r["fisher_q"]),
            "strict_lps_profile_clusters_iss_present": f"{iss_clusters_present}/{iss_clusters_total}",
            "strict_lps_profile_clusters_earth_present": f"{earth_clusters_present}/{earth_clusters_total}",
            "cluster_collapsed_fisher_q": fmt_float(r["strict_profile_cluster_q"]),
            "minimum_iss_year_fraction_present": fmt_float(r["min_iss_year_fraction"]),
            "iss_year_persistence": year_summary(year, r["gene"]),
            "leave_one_iss_year_out_max_q": fmt_float(r["leave_one_year_max_q"]),
            "species_restricted_phylogeny_aware_empirical_q": fmt_float(r["species_restricted_empirical_q"]),
            "species_restricted_phylogeny_aware_directional_q": fmt_float(r["species_restricted_directional_q"]),
            "phylogeny_aware_interpretation": r["interpretation"],
        }
        out["overall_interpretation"] = overall_interpretation(out)
        rows.append(out)

    final = pd.DataFrame(rows)
    if OUT.exists() and not BACKUP.exists():
        shutil.copy2(OUT, BACKUP)
    tmp = OUT.with_suffix(".tmp")
    final.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, OUT)


def main() -> None:
    rebuild_table4()
    print(OUT)
    print(BACKUP)


if __name__ == "__main__":
    main()
