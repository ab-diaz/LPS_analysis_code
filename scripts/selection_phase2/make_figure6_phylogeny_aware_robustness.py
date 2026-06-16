#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
ROUTE = ROOT / "pantoea_route1_strengthening"
PHYLO = ROOT / "phylogeny_aware_association"
TARGET_GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL", "wzm", "wzt", "rfaZ"]
ISS_COLOR = "#b23a48"
EARTH_COLOR = "#3d6f9f"


def add_panel_label(ax, label, title):
    ax.text(0.0, 1.03, label, transform=ax.transAxes, ha="left", va="bottom", weight="bold", fontsize=10)


def main():
    FIG.mkdir(exist_ok=True)
    clusters = pd.read_csv(ROUTE / "cluster_collapsed_strict_lps_fisher.tsv", sep="\t")
    years = pd.read_csv(ROUTE / "iss_year_persistence_target_genes.tsv", sep="\t")
    phylo = pd.read_csv(PHYLO / "phylogeny_aware_gene_association.tsv", sep="\t")

    genes = [g for g in TARGET_GENES if g in set(clusters["gene_id"])]
    fig = plt.figure(figsize=(13.4, 4.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.42, 1.05, 1.25], wspace=0.48)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])

    ct = clusters.set_index("gene_id").reindex(genes).reset_index()
    x = np.arange(len(ct))
    iss_frac = ct["iss_profile_clusters_present"] / (
        ct["iss_profile_clusters_present"] + ct["iss_profile_clusters_absent"]
    )
    earth_frac = ct["earth_profile_clusters_present"] / (
        ct["earth_profile_clusters_present"] + ct["earth_profile_clusters_absent"]
    )
    ax0.bar(x - 0.18, iss_frac, width=0.36, label="ISS profile clusters", color=ISS_COLOR)
    ax0.bar(x + 0.18, earth_frac, width=0.36, label="Earth profile clusters", color=EARTH_COLOR)
    ax0.set_xticks(x)
    ax0.set_xticklabels(genes, rotation=45, ha="right")
    ax0.set_ylabel("Fraction of exact LPS-profile clusters")
    ax0.set_ylim(0, 1.08)
    add_panel_label(ax0, "A", "Founder-aware gene presence")
    ax0.legend(
        frameon=False,
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(0.0, 1.02),
        ncol=1,
        borderaxespad=0,
    )
    ax0.spines[["top", "right"]].set_visible(False)

    keep_year = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL", "wzm"]
    subset = years[years["gene_id"].isin(keep_year)].copy()
    pivot = (
        subset.pivot_table(index="gene_id", columns="year", values="fraction_present")
        .reindex(keep_year)
        .sort_index(axis=1)
    )
    im = ax1.imshow(pivot.fillna(0).to_numpy(), aspect="auto", vmin=0, vmax=1, cmap="Reds")
    ax1.set_yticks(np.arange(len(pivot.index)))
    ax1.set_yticklabels(pivot.index)
    ax1.set_xticks(np.arange(len(pivot.columns)))
    ax1.set_xticklabels([str(int(c)) for c in pivot.columns], rotation=0)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.iloc[i, j]
            ax1.text(j, i, "" if pd.isna(v) else f"{v:.2f}", ha="center", va="center", fontsize=8)
    add_panel_label(ax1, "B", "ISS year-level persistence")
    cbar = fig.colorbar(im, ax=ax1, shrink=0.78, pad=0.02)
    cbar.set_label("fraction present", fontsize=8)
    cbar.ax.tick_params(labelsize=8)

    pa = phylo.set_index("gene").reindex(TARGET_GENES).reset_index()
    y = np.arange(len(pa))
    q = pa["species_restricted_empirical_q"].clip(lower=1e-4)
    score = -np.log10(q)
    colors = np.where(pa["observed_prevalence_difference"] >= 0, ISS_COLOR, EARTH_COLOR)
    ax2.barh(y, score, color=colors, alpha=0.9)
    ax2.axvline(-np.log10(0.05), color="#333333", linewidth=1, linestyle="--")
    ax2.text(
        -np.log10(0.05) + 0.03,
        -0.55,
        "q=0.05",
        ha="left",
        va="bottom",
        fontsize=8,
        color="#333333",
    )
    ax2.set_yticks(y)
    ax2.set_yticklabels(pa["gene"])
    ax2.invert_yaxis()
    ax2.set_xlabel("-log10 species-restricted q")
    ax2.set_xlim(0, max(4.2, score.max() + 0.35))
    add_panel_label(ax2, "C", "Phylogeny-aware association")
    ax2.spines[["top", "right"]].set_visible(False)
    for yi, row in pa.iterrows():
        direction = "ISS" if row["observed_prevalence_difference"] >= 0 else "Earth"
        label = f"{direction}; {row['iss_present']}/{row['iss_total']} vs {row['earth_present']}/{row['earth_total']}"
        ax2.text(score.iloc[yi] + 0.06, yi, label, va="center", ha="left", fontsize=7.5)

    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"figure6_founder_year_persistence.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    print(FIG / "figure6_founder_year_persistence.png")


if __name__ == "__main__":
    main()
