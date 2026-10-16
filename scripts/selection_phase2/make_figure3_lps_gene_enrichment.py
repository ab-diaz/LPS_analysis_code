#!/usr/bin/env python3
"""Regenerate Figure 3 after collapsing equivalent O-antigen ligase labels.

waaL, rfaL and KO:K02847 are aliases for the same O-antigen ligase feature in
the annotation tables.  This script merges them by genome-level logical OR,
recalculates Fisher tests and BH q-values, and makes a publication-ready
heatmap/volcano figure without modifying the original figure or tables.
"""
import argparse
from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.stats import fisher_exact

LPS = Path(".")
FIGDIR = Path(".")
OUT_TSV = Path("figure3_waaL_merged_fisher.tsv")

ALIASES = ["waaL", "rfaL", "ko:K02847"]
MERGED_ID = "waaL/rfaL (KO:K02847)"
PANEL_GENES = [
    "rfbC", "rfbA", "rfbD", MERGED_ID, "rfbB", "wzm", "wzt",
    "rfaZ", "rfbX", "rfbI", "wzx",
]


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    m = len(p)
    order = np.argsort(p)
    ranked = p[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m, dtype=float)
    out[order] = np.minimum(ranked, 1.0)
    return out


def load_tables():
    pres = pd.read_csv(LPS / "lps_strict_presence.tsv", sep="\t")
    long = pd.read_csv(LPS / "lps_strict_presence_long.tsv", sep="\t")
    group_map = long[["isolate", "group"]].drop_duplicates().set_index("isolate")["group"].to_dict()
    isolate_cols = [c for c in pres.columns if c not in {"gene_id", "kegg_ko", "description"}]

    # Build one logical presence row for the three equivalent annotations.
    alias = pres[pres["gene_id"].isin(ALIASES)].set_index("gene_id")[isolate_cols].astype(int)
    merged = (alias.loc[ALIASES].max(axis=0) > 0).astype(int)
    merged_row = pd.Series(merged, name=MERGED_ID)

    # Preserve the original panel order, replacing the three aliases by one row.
    rows = {}
    for gid in PANEL_GENES:
        if gid == MERGED_ID:
            rows[gid] = merged_row
        elif gid in set(pres["gene_id"]):
            rows[gid] = pres[pres["gene_id"] == gid].iloc[0][isolate_cols].astype(int)
    mat = pd.DataFrame(rows).T[isolate_cols]
    earth = sorted([c for c in isolate_cols if group_map.get(c) == "Earth"])
    iss = sorted([c for c in isolate_cols if group_map.get(c) == "ISS"])
    ordered = earth + iss

    # Recalculate all feature-level tests after replacing the three aliases.
    feature_rows = []
    for _, row in pres.iterrows():
        gid = row["gene_id"]
        if gid in ALIASES:
            continue
        vals = row[isolate_cols].astype(int)
        feature_rows.append({
            "gene_id": gid,
            "kegg_ko": row.get("kegg_ko", ""),
            "description": row.get("description", ""),
            "values": vals,
        })
    feature_rows.append({
        "gene_id": MERGED_ID,
        "kegg_ko": "ko:K02847",
        "description": "O-antigen ligase (waaL/rfaL; KO:K02847)",
        "values": merged,
    })
    stats = []
    for item in feature_rows:
        vals = item["values"]
        ip = int(vals[iss].sum()); ia = int(len(iss) - ip)
        ep = int(vals[earth].sum()); ea = int(len(earth) - ep)
        odds, pval = fisher_exact([[ip, ia], [ep, ea]], alternative="two-sided")
        stats.append({
            "gene_id": item["gene_id"], "kegg_ko": item["kegg_ko"],
            "description": item["description"], "iss_present": ip,
            "iss_absent": ia, "earth_present": ep, "earth_absent": ea,
            "odds_ratio": odds, "p_value": pval,
        })
    stats = pd.DataFrame(stats)
    stats["q_value"] = bh_adjust(stats["p_value"].to_numpy())
    # Stable, readable output; represent infinite odds ratios as inf.
    stats.to_csv(OUT_TSV, sep="\t", index=False, float_format="%.9g")
    return mat[ordered], len(earth), stats


def plot_heatmap(ax, mat, n_earth):
    cmap = ListedColormap(["#f1f1f1", "#1f6f8b"])
    ax.imshow(mat.values, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=1)
    labels = ["waaL/rfaL\n(KO:K02847)" if x == MERGED_ID else x for x in mat.index]
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xticks([])
    ax.text(0.0, 1.06, "A", transform=ax.transAxes, ha="left", va="bottom", fontsize=11, weight="bold")
    ax.set_ylabel("Feature")
    ax.axvline(n_earth - 0.5, color="#222222", linewidth=1.0)
    ax.text((n_earth / mat.shape[1]) / 2, 1.004, "Earth", transform=ax.transAxes, ha="center", va="bottom", fontsize=9)
    ax.text((n_earth + (mat.shape[1] - n_earth) / 2) / mat.shape[1], 1.004, "ISS", transform=ax.transAxes, ha="center", va="bottom", fontsize=9)
    ax.legend(handles=[
        Patch(facecolor="#f1f1f1", edgecolor="#999999", label="Absent"),
        Patch(facecolor="#1f6f8b", label="Present")], frameon=False,
        loc="lower right", bbox_to_anchor=(1.0, -0.24), ncol=2, fontsize=8)


def prepare_volcano(stats):
    df = stats.copy()
    df["odds_ratio_num"] = pd.to_numeric(df["odds_ratio"], errors="coerce")
    finite = df.loc[np.isfinite(df["odds_ratio_num"]) & (df["odds_ratio_num"] > 0), "odds_ratio_num"]
    finite_logs = np.log2(finite)
    max_abs = max(4, math.ceil(np.nanmax(np.abs(finite_logs))) + 1)
    min_log, max_log = -max_abs, max_abs

    def log_or(x):
        if x == 0: return min_log
        if np.isinf(x): return max_log
        return np.log2(x) if x > 0 else np.nan
    df["log2_or"] = df["odds_ratio_num"].map(log_or)
    q = pd.to_numeric(df["q_value"], errors="coerce")
    min_pos = q[q > 0].min()
    df["neglog10_q"] = -np.log10(q.mask(q <= 0, min_pos / 10 if pd.notna(min_pos) else 1e-6))
    df["direction"] = "Not significant"
    df.loc[(q <= 0.05) & (df["log2_or"] > 0), "direction"] = "ISS-enriched"
    df.loc[(q <= 0.05) & (df["log2_or"] < 0), "direction"] = "Earth-enriched"
    return df, min_log, max_log


def plot_volcano(ax, df, min_log, max_log):
    colors = {"ISS-enriched": "#1f6f8b", "Earth-enriched": "#c85a3a", "Not significant": "#b8b8b8"}
    for direction, sub in df.groupby("direction"):
        ax.scatter(sub["log2_or"], sub["neglog10_q"], s=28, alpha=0.85, c=colors[direction], edgecolors="none", label=direction)
    ax.axhline(-np.log10(0.05), color="#555555", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_xlim(min_log - 0.5, max_log + 0.5)
    ax.set_xlabel("log2 odds ratio (ISS / Earth)")
    ax.set_ylabel("-log10(q value)")
    ax.text(0.0, 1.02, "B", transform=ax.transAxes, ha="left", va="bottom", fontsize=11, weight="bold")
    ax.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, borderaxespad=0, handletextpad=0.4, columnspacing=0.9)
    offsets = {"rfbA": (-28, -18), "rfbB": (6, 6), "rfbC": (8, -18), "rfbD": (-6, -32), MERGED_ID: (8, -32), "wzm": (6, 4), "wzt": (6, 4), "rfaZ": (6, -2)}
    for label, offset in offsets.items():
        hit = df[df["gene_id"] == label]
        if hit.empty: continue
        row = hit.iloc[0]
        ax.annotate("waaL/rfaL (KO:K02847)" if label == MERGED_ID else label,
                    xy=(row["log2_or"], row["neglog10_q"]), xytext=offset,
                    textcoords="offset points", fontsize=7,
                    ha="left" if offset[0] >= 0 else "right", va="bottom",
                    arrowprops={"arrowstyle": "-", "lw": 0.4, "color": "#555555"})
    ax.text(max_log, 0.08, "inf", fontsize=7, ha="center", va="bottom")
    ax.text(min_log, 0.08, "0", fontsize=7, ha="center", va="bottom")


def main():
    global LPS, FIGDIR, OUT_TSV
    ap = argparse.ArgumentParser(description="Reproduce the merged-ligase manuscript Figure 3")
    ap.add_argument("--analysis-dir", type=Path, required=True,
                    help="Directory containing lps_strict_presence*.tsv")
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()
    LPS = args.analysis_dir.resolve()
    FIGDIR = args.outdir.resolve()
    OUT_TSV = FIGDIR / "figure3_waaL_merged_fisher.tsv"
    FIGDIR.mkdir(parents=True, exist_ok=True)
    mat, n_earth, stats = load_tables()
    volcano, min_log, max_log = prepare_volcano(stats)
    fig = plt.figure(figsize=(12.8, 7.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.32)
    plot_heatmap(fig.add_subplot(gs[0, 0]), mat, n_earth)
    plot_volcano(fig.add_subplot(gs[0, 1]), volcano, min_log, max_log)
    fig.savefig(FIGDIR / "Figure_3.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGDIR / "Figure_3.pdf", bbox_inches="tight")
    print(FIGDIR / "Figure_3.png")
    print(FIGDIR / "Figure_3.pdf")
    print(OUT_TSV)


if __name__ == "__main__":
    main()
