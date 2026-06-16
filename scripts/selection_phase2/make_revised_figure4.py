#!/usr/bin/env python3
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


LPS = Path("results/lps_analysis")
OUTDIR = LPS / "selection_phase2" / "figures"

PANEL_GENES = [
    "rfbC",
    "rfbA",
    "rfbD",
    "waaL",
    "rfbB",
    "wzm",
    "wzt",
    "ko:K02847",
    "rfaZ",
    "rfbX",
    "rfbI",
    "rfaL",
    "wzx",
]


def load_pantoea_heatmap():
    pres = pd.read_csv(LPS / "lps_strict_presence.tsv", sep="\t")
    long = pd.read_csv(LPS / "lps_strict_presence_long.tsv", sep="\t")
    group_map = long[["isolate", "group"]].drop_duplicates().set_index("isolate")["group"].to_dict()

    isolate_cols = [c for c in pres.columns if c not in {"gene_id", "kegg_ko", "description"}]
    wanted = [g for g in PANEL_GENES if g in set(pres["gene_id"])]
    mat = pres[pres["gene_id"].isin(wanted)].set_index("gene_id")[isolate_cols]
    mat = mat.loc[wanted]

    earth = [c for c in isolate_cols if group_map.get(c) == "Earth"]
    iss = [c for c in isolate_cols if group_map.get(c) == "ISS"]
    earth = sorted(earth)
    iss = sorted(iss)
    ordered = earth + iss
    return mat[ordered].astype(int), len(earth), len(iss)


def plot_heatmap(ax, mat, n_earth):
    cmap = ListedColormap(["#f1f1f1", "#1f6f8b"])
    ax.imshow(mat.values, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(mat.index, fontsize=8)
    ax.set_xticks([])
    ax.text(0.0, 1.06, "A", transform=ax.transAxes, ha="left", va="bottom", fontsize=11, weight="bold")
    ax.set_ylabel("Gene")
    ax.axvline(n_earth - 0.5, color="#222222", linewidth=1.0)
    ax.text((n_earth / 2) / mat.shape[1], 1.004, "Earth", transform=ax.transAxes, ha="center", va="bottom", fontsize=9)
    ax.text((n_earth + (mat.shape[1] - n_earth) / 2) / mat.shape[1], 1.004, "ISS", transform=ax.transAxes, ha="center", va="bottom", fontsize=9)
    ax.legend(
        handles=[Patch(facecolor="#f1f1f1", edgecolor="#999999", label="Absent"), Patch(facecolor="#1f6f8b", label="Present")],
        frameon=False,
        loc="lower right",
        bbox_to_anchor=(1.0, -0.22),
        ncol=2,
        fontsize=8,
    )


def load_volcano():
    df = pd.read_csv(LPS / "lps_strict_gene_fisher.tsv", sep="\t")
    df["odds_ratio_num"] = pd.to_numeric(df["odds_ratio"].replace({"inf": np.inf}), errors="coerce")
    finite = df.loc[np.isfinite(df["odds_ratio_num"]) & (df["odds_ratio_num"] > 0), "odds_ratio_num"]
    finite_logs = np.log2(finite)
    max_abs = max(4, math.ceil(np.nanmax(np.abs(finite_logs))) + 1)
    min_log = -max_abs
    max_log = max_abs

    def log_or(x):
        if x == 0:
            return min_log
        if np.isinf(x):
            return max_log
        if x > 0:
            return np.log2(x)
        return np.nan

    df["log2_or"] = df["odds_ratio_num"].map(log_or)
    q = pd.to_numeric(df["q_value"], errors="coerce")
    min_pos = q[q > 0].min()
    q_plot = q.mask(q <= 0, min_pos / 10 if pd.notna(min_pos) else 1e-6)
    df["neglog10_q"] = -np.log10(q_plot)
    df["direction"] = "Not significant"
    df.loc[(q <= 0.05) & (df["log2_or"] > 0), "direction"] = "ISS-enriched"
    df.loc[(q <= 0.05) & (df["log2_or"] < 0), "direction"] = "Earth-enriched"
    return df, min_log, max_log


def plot_volcano(ax, df, min_log, max_log):
    colors = {
        "ISS-enriched": "#1f6f8b",
        "Earth-enriched": "#c85a3a",
        "Not significant": "#b8b8b8",
    }
    for direction, sub in df.groupby("direction"):
        ax.scatter(sub["log2_or"], sub["neglog10_q"], s=28, alpha=0.85, c=colors[direction], edgecolors="none", label=direction)

    ax.axhline(-np.log10(0.05), color="#555555", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_xlim(min_log - 0.5, max_log + 0.5)
    ax.set_xlabel("log2 odds ratio (ISS / Earth)")
    ax.set_ylabel("-log10(q value)")
    ax.text(0.0, 1.02, "B", transform=ax.transAxes, ha="left", va="bottom", fontsize=11, weight="bold")
    ax.legend(
        frameon=False,
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=3,
        borderaxespad=0,
        handletextpad=0.4,
        columnspacing=0.9,
    )

    label_offsets = {
        "rfbA": (-28, -18),
        "rfbB": (6, 6),
        "rfbC": (8, -18),
        "rfbD": (-6, -32),
        "waaL": (8, -32),
        "wzm": (6, 4),
        "wzt": (6, 4),
        "rfaZ": (6, -2),
        "ko:K02847": (6, 8),
    }
    labels = list(label_offsets)
    for label in labels:
        hit = df[df["gene_id"] == label]
        if hit.empty:
            continue
        row = hit.iloc[0]
        x, y = row["log2_or"], row["neglog10_q"]
        ax.annotate(
            label,
            xy=(x, y),
            xytext=label_offsets[label],
            textcoords="offset points",
            fontsize=7,
            ha="left" if label_offsets[label][0] >= 0 else "right",
            va="bottom",
            arrowprops={"arrowstyle": "-", "lw": 0.4, "color": "#555555"},
        )

    ax.text(max_log, 0.08, "inf", fontsize=7, ha="center", va="bottom")
    ax.text(min_log, 0.08, "0", fontsize=7, ha="center", va="bottom")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    mat, n_earth, _ = load_pantoea_heatmap()
    volcano, min_log, max_log = load_volcano()

    fig = plt.figure(figsize=(12.8, 7.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.32)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    plot_heatmap(ax_a, mat, n_earth)
    plot_volcano(ax_b, volcano, min_log, max_log)

    fig.savefig(OUTDIR / "figure4_lps_oantigen_gene_enrichment.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "figure4_lps_oantigen_gene_enrichment.pdf", bbox_inches="tight")
    print(OUTDIR / "figure4_lps_oantigen_gene_enrichment.png")
    print(OUTDIR / "figure4_lps_oantigen_gene_enrichment.pdf")


if __name__ == "__main__":
    main()
