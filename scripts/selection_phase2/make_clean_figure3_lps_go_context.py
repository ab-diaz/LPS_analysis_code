#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


BASE = Path("results/lps_analysis")
OUT = Path("figures")
OUT.mkdir(exist_ok=True)


def short_go_label(label: str) -> str:
    if " " not in label:
        return label
    go, term = label.split(" ", 1)
    term = term.replace("carbohydrate derivative metabolic process", "carbohydrate-derivative metabolism")
    term = term.replace("carbohydrate metabolic process", "carbohydrate metabolism")
    term = term.replace("small molecule biosynthetic process", "small-molecule biosynthesis")
    term = term.replace("protein-containing complex assembly", "protein-complex assembly")
    term = term.replace("cell wall organization or biogenesis", "cell-wall organization/biogenesis")
    term = term.replace("nucleobase-containing small molecule metabolic process", "nucleobase-small-molecule metabolism")
    return f"{term}\n{go}"


def main() -> None:
    counts = pd.read_csv(BASE / "lps_go_category_counts.tsv", sep="\t")
    goslim = pd.read_csv(BASE / "lps_goslim_norm_metadata_ordered.tsv", sep="\t")
    goslim = goslim.rename(columns={"Unnamed: 0": "isolate"})

    # One row was produced with an older sample-name normalization; the genome is
    # the same isolate used elsewhere in the manuscript.
    goslim["isolate"] = goslim["isolate"].replace(
        {"Pantoea_piersonii_F9_6S_D2_EB4_spades": "Pantoea_piersonii_F8_6S_D2_EB4_spades"}
    )

    meta = counts[["isolate", "group"]].drop_duplicates()
    goslim = goslim.merge(meta, on="isolate", how="left")
    missing = goslim["group"].isna().sum()
    if missing:
        raise RuntimeError(f"Missing group labels for {missing} GO-slim rows")

    go_cols = [c for c in goslim.columns if c.startswith("GO:")]
    top_cols = goslim[go_cols].mean().sort_values(ascending=False).head(12).index.tolist()

    long = goslim.melt(
        id_vars=["isolate", "group"],
        value_vars=top_cols,
        var_name="category",
        value_name="normalized_count",
    )
    summary = (
        long.groupby(["category", "group"], as_index=False)
        .agg(mean=("normalized_count", "mean"), sem=("normalized_count", "sem"))
    )

    order = (
        long.groupby("category")["normalized_count"]
        .mean()
        .sort_values(ascending=True)
        .index.tolist()
    )

    palette = {"Earth": "#4C78A8", "ISS": "#7B61A9"}
    sns.set_theme(style="whitegrid", context="paper", font_scale=1.0)
    fig = plt.figure(figsize=(9.2, 5.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1.0])

    ax = fig.add_subplot(gs[0, 0])
    ypos = np.arange(len(order))
    offsets = {"Earth": -0.12, "ISS": 0.12}
    for group in ["Earth", "ISS"]:
        sdf = summary[summary["group"] == group].set_index("category").loc[order]
        ax.errorbar(
            sdf["mean"],
            ypos + offsets[group],
            xerr=sdf["sem"].fillna(0),
            fmt="o",
            ms=4.5,
            lw=1.2,
            capsize=2,
            color=palette[group],
            label=f"{group} mean +/- SEM",
        )
    for cat_i, cat in enumerate(order):
        vals = summary[summary["category"] == cat].set_index("group")["mean"]
        if {"Earth", "ISS"}.issubset(vals.index):
            ax.plot([vals["Earth"], vals["ISS"]], [cat_i, cat_i], color="#D6D6D6", lw=1, zorder=0)
    ax.set_yticks(ypos)
    ax.set_yticklabels([short_go_label(x) for x in order], fontsize=7)
    ax.set_xlabel("Normalized count per LPS-related gene")
    ax.set_ylabel("")
    ax.set_title("A  Broad LPS-associated GO-slim profiles", loc="left", fontweight="bold")
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.grid(axis="y", visible=False)

    ax2 = fig.add_subplot(gs[0, 1])
    sns.boxplot(
        data=counts,
        x="group",
        y="total_unique_genes",
        order=["Earth", "ISS"],
        palette=palette,
        width=0.55,
        fliersize=0,
        ax=ax2,
    )
    sns.stripplot(
        data=counts,
        x="group",
        y="total_unique_genes",
        order=["Earth", "ISS"],
        color="#222222",
        size=2.4,
        alpha=0.55,
        jitter=0.18,
        ax=ax2,
    )
    n_by_group = counts.groupby("group")["isolate"].nunique().to_dict()
    ax2.set_xticklabels([f"Earth\nn={n_by_group.get('Earth', 0)}", f"ISS\nn={n_by_group.get('ISS', 0)}"])
    ax2.set_xlabel("")
    ax2.set_ylabel("Unique LPS-related GO-category genes per genome")
    ax2.set_title("B  Per-genome LPS GO-category content", loc="left", fontweight="bold")
    ax2.grid(axis="x", visible=False)

    fig.suptitle("LPS-related GO-category context across Pantoea genomes", y=1.02, fontsize=12, fontweight="bold")
    for ext in ["png", "pdf"]:
        fig.savefig(OUT / f"figure3_lps_go_category_context.{ext}", dpi=400, bbox_inches="tight")
    plt.close(fig)

    print(OUT / "figure3_lps_go_category_context.png")
    print(OUT / "figure3_lps_go_category_context.pdf")


if __name__ == "__main__":
    main()
