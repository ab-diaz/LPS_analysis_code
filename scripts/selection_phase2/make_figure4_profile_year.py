#!/usr/bin/env python3
"""Render the submitted two-panel Figure 4 from frozen source tables.

No permutation or phylogenetic-adjustment statistic is calculated here.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL", "wzm", "wzt", "rfaZ"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True, type=Path,
                    help="cluster_collapsed_strict_lps_fisher.tsv")
    ap.add_argument("--years", required=True, type=Path,
                    help="iss_year_persistence_target_genes.tsv")
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    cl = pd.read_csv(args.clusters, sep="\t")
    years = pd.read_csv(args.years, sep="\t")
    assert set(GENES).issubset(set(cl.gene_id))
    ct = cl.set_index("gene_id").reindex(GENES)
    x = np.arange(len(ct))
    iss = ct.iss_profile_clusters_present / (ct.iss_profile_clusters_present + ct.iss_profile_clusters_absent)
    earth = ct.earth_profile_clusters_present / (ct.earth_profile_clusters_present + ct.earth_profile_clusters_absent)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10.0, 4.6),
                                   gridspec_kw={"width_ratios": [1.55, 1.05]})
    fig.subplots_adjust(left=.08, right=.96, bottom=.20, top=.89, wspace=.38)
    ax0.bar(x-.18, iss, width=.36, color="#b23a48", label="ISS profile clusters")
    ax0.bar(x+.18, earth, width=.36, color="#3d6f9f", label="Earth profile clusters")
    ax0.set_xticks(x, GENES, rotation=45, ha="right")
    ax0.set_ylabel("Fraction of exact LPS-profile clusters")
    ax0.set_ylim(0, 1.08)
    ax0.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0, 1.02))
    ax0.text(0, 1.035, "A", transform=ax0.transAxes, weight="bold", fontsize=11)
    ax0.spines[["top", "right"]].set_visible(False)

    selected = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL", "wzm"]
    subset = years[years.gene_id.isin(selected)]
    pivot = subset.pivot_table(index="gene_id", columns="year", values="fraction_present").reindex(selected).sort_index(axis=1)
    im = ax1.imshow(pivot.fillna(0).to_numpy(), aspect="auto", vmin=0, vmax=1, cmap="Reds")
    ax1.set_yticks(np.arange(len(pivot)), pivot.index)
    ax1.set_xticks(np.arange(len(pivot.columns)), [str(int(c)) for c in pivot.columns])
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.iloc[i, j]
            ax1.text(j, i, "" if pd.isna(val) else f"{val:.2f}", ha="center", va="center",
                     fontsize=8, color="white", weight="semibold")
    ax1.text(0, 1.035, "B", transform=ax1.transAxes, weight="bold", fontsize=11)
    cb = fig.colorbar(im, ax=ax1, shrink=.76, pad=.03)
    cb.set_label("Fraction present", fontsize=8)
    args.outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.outdir / "Figure_4.pdf", bbox_inches="tight")
    fig.savefig(args.outdir / "Figure_4.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
