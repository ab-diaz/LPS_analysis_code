from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent
OUTDIR = BASE / "figures"
SENS = BASE / "relax_sensitivity_20260530" / "RELAX_sensitivity_key_comparison.tsv"
MUT = BASE / "mutation_density_iss.tsv"
ORDER = BASE / "mutation_density_isolate_order.tsv"


def short_isolate(row):
    year = ""
    if pd.notna(row["year"]):
        year = str(int(float(row["year"])))
    label = row["display_isolate"].replace("Pantoea_piersonii_", "")
    return f"{year} {label}".strip()


def main():
    OUTDIR.mkdir(exist_ok=True)

    sens = pd.read_csv(SENS, sep="\t")
    sens = sens[sens["gene"].isin(["rfbD", "waaL"])].copy()
    sens["gene"] = pd.Categorical(sens["gene"], ["rfbD", "waaL"], ordered=True)
    sens = sens.sort_values("gene")

    mut = pd.read_csv(MUT, sep="\t")
    order = pd.read_csv(ORDER, sep="\t")
    genes = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL"]
    isolates = order["isolate"].tolist()
    labels = dict(zip(order["isolate"], order.apply(short_isolate, axis=1)))

    heat = (
        mut.pivot_table(
            index="isolate",
            columns="gene",
            values="mismatches_per_kb",
            aggfunc="first",
        )
        .reindex(index=isolates, columns=genes)
        .fillna(0)
    )

    fig = plt.figure(figsize=(9.2, 12.0))
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1.05, 2.55], hspace=0.34)

    ax = fig.add_subplot(gs[0])
    x = np.arange(len(sens))
    w = 0.24
    bars = [
        ("original_k", "original_p", "Original\nISS terminals", "#b23a48"),
        ("minus_divergent_k", "minus_divergent_p", "Minus divergent\nISS branch", "#8b9797"),
        ("single_divergent_k", "single_divergent_p", "Divergent branch\nonly", "#3778b8"),
    ]
    for i, (kcol, pcol, label, color) in enumerate(bars):
        vals = sens[kcol].astype(float).to_numpy()
        ps = sens[pcol].astype(float).to_numpy()
        xpos = x + (i - 1) * w
        ax.bar(xpos, vals, width=w, color=color, edgecolor="#222222", linewidth=0.9, label=label)
        for xi, val, p in zip(xpos, vals, ps):
            txt = "p<0.001" if p < 0.001 else f"p={p:.3g}"
            ax.text(xi, val * 1.09, txt, ha="center", va="bottom", fontsize=8)

    ax.axhline(1, color="#333333", linestyle="--", linewidth=1.0)
    ax.text(1.43, 1.06, "k = 1", fontsize=8, ha="left", va="bottom")
    ax.set_yscale("log")
    ax.set_ylim(0.7, 45)
    ax.set_xticks(x)
    ax.set_xticklabels(sens["gene"], fontsize=11)
    ax.set_ylabel("RELAX selection intensity k (log scale)", fontsize=10)
    ax.text(0.0, 1.04, "A", transform=ax.transAxes, ha="left", va="bottom", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax2 = fig.add_subplot(gs[1])
    cmap = LinearSegmentedColormap.from_list(
        "founder_mismatches",
        ["#f2f2f2", "#b9d7dc", "#2f6f9f", "#8b1e4b"],
    )
    vmax = max(160, float(np.nanmax(heat.to_numpy())))
    im = ax2.imshow(heat.to_numpy(), aspect="auto", cmap=cmap, vmin=0, vmax=vmax)
    ax2.set_xticks(np.arange(len(genes)))
    ax2.set_yticks(np.arange(len(heat.index)))
    ax2.set_xticks(np.arange(-0.5, len(genes), 1), minor=True)
    ax2.set_yticks(np.arange(-0.5, len(heat.index), 1), minor=True)
    ax2.grid(which="minor", color="#ffffff", linewidth=0.35)
    ax2.tick_params(which="minor", bottom=False, left=False)
    cbar = fig.colorbar(im, ax=ax2, shrink=0.68, pad=0.02)
    cbar.set_label("Mismatches/kb vs founder")
    ax2.text(0.0, 1.04, "B", transform=ax2.transAxes, ha="left", va="bottom", fontsize=13, fontweight="bold")
    ax2.set_xlabel("")
    ax2.set_ylabel("")
    ax2.set_xticklabels(genes, rotation=0, fontsize=10)
    ax2.xaxis.tick_top()
    ax2.tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False)
    ax2.set_yticklabels([labels.get(i, i) for i in heat.index], rotation=0, fontsize=7)

    # Match the year breaks used in the earlier figure.
    for y in [3, 10, 38]:
        ax2.axhline(y, color="#333333", linewidth=0.8)

    fig.subplots_adjust(left=0.19, right=0.92, top=0.96, bottom=0.04)
    for ext in ["png", "pdf"]:
        out = OUTDIR / f"figure7_relax_gard_sensitivity.{ext}"
        fig.savefig(out, dpi=300 if ext == "png" else None)
        print(out)


if __name__ == "__main__":
    main()
