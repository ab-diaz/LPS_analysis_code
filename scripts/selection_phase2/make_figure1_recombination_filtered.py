#!/usr/bin/env python3
"""Render the manuscript's 57-genome Gubbins-filtered Figure 1.

The script only renders an already inferred tree; it does not alter topology
or infer ancestral states. All inputs are explicit so it runs outside the
original workstation.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL", "wzm", "wzt", "rfaZ"]
SOURCE_COLORS = {"ISS": "#277da8", "Earth": "#a4a8ad"}
YEAR_COLORS = {"2016": "#7b3294", "2018": "#008837", "2021": "#e66101", "2022": "#5e3c99", "": "#ffffff"}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def support(name: str | None) -> tuple[float, float] | None:
    try:
        values = [float(v) for v in str(name).split("/")]
        return tuple(values[-2:]) if len(values) >= 2 else None
    except (TypeError, ValueError):
        return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=False)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    meta = {r["genome"]: dict(r) for r in rows(args.metadata)}
    if args.validation:
        for r in rows(args.validation):
            m = meta.setdefault(r["isolate"], {})
            m.update({"source": r["source"], "year": r.get("year", "")})
            if r.get("gene") in GENES:
                m[r["gene"]] = r["present"]

    tree = Phylo.read(args.tree, "newick")
    if any(t.name == "Reference" for t in tree.get_terminals()):
        raise ValueError("The matched-tree input must already exclude the technical Reference")
    tree.ladderize()
    names = [t.name for t in tree.get_terminals()]
    if len(names) != 57 or len(set(names)) != 57:
        raise ValueError(f"Expected 57 unique P. piersonii tips, found {len(names)}")
    missing = sorted(set(names) - set(meta))
    if missing:
        raise ValueError(f"Tree tips missing from metadata: {missing[:5]}")

    depths = tree.depths()
    ypos = {tip: i for i, tip in enumerate(tree.get_terminals())}
    for clade in tree.get_nonterminals(order="postorder"):
        ypos[clade] = float(np.mean([ypos[c] for c in clade.clades]))

    fig = plt.figure(figsize=(11.5, 12.5))
    gs = fig.add_gridspec(1, 4, width_ratios=[3.8, 3.0, .7, 2.0], wspace=.06,
                          left=.05, right=.98, top=.94, bottom=.12)
    ax, labels, strip, heat = [fig.add_subplot(gs[0, i]) for i in range(4)]
    for clade in tree.get_nonterminals():
        ax.plot([depths[clade]] * 2,
                [min(ypos[c] for c in clade.clades), max(ypos[c] for c in clade.clades)],
                color="#424b55", lw=.65)
        for child in clade.clades:
            ax.plot([depths[clade], depths[child]], [ypos[child]] * 2,
                    color="#424b55", lw=.65)
        s = support(clade.name)
        if s and s[0] >= 80 and s[1] >= 95:
            ax.scatter(depths[clade], ypos[clade], s=7, color="#a65025", zorder=3)
    for i, name in enumerate(names):
        labels.text(0, i, name.removeprefix("Pantoea_piersonii_"), fontsize=7, va="center")
    right = max(depths[t] for t in tree.get_terminals()) * 1.02
    for tip in tree.get_terminals():
        ax.plot([depths[tip], right], [ypos[tip]] * 2, color="#c6cbd0", lw=.4, ls=":", zorder=0)

    source = [[matplotlib.colors.to_rgb(SOURCE_COLORS.get(meta[n].get("source", ""), "#a4a8ad")),
               matplotlib.colors.to_rgb(YEAR_COLORS.get(str(meta[n].get("year", "")).removesuffix(".0"), "#ffffff"))]
              for n in names]
    strip.imshow(source, aspect="auto", origin="upper")
    strip.set_xticks([0, 1], ["Source", "Year"], rotation=90, fontsize=7)
    strip.xaxis.tick_top()
    matrix = np.array([[int(meta[n].get(g, 0)) for g in GENES] for n in names])
    heat.imshow(matrix, aspect="auto", origin="upper", cmap=ListedColormap(["#f0f1f2", "#277da8"]), vmin=0, vmax=1)
    heat.set_xticks(np.arange(-.5, len(GENES), 1), minor=True)
    heat.set_yticks(np.arange(-.5, len(names), 1), minor=True)
    heat.grid(which="minor", color="white", lw=.25)
    heat.tick_params(which="minor", length=0)
    heat.set_xticks(range(len(GENES)), GENES, rotation=55, ha="left", fontsize=8, fontstyle="italic")
    heat.xaxis.tick_top()
    for panel in [ax, labels, strip, heat]:
        panel.set_ylim(len(names) - .4, -.6)
    labels.axis("off")
    for panel in [strip, heat]:
        panel.tick_params(length=0)
        for spine in panel.spines.values(): spine.set_visible(False)
    for spine in ["top", "right", "left"]: ax.spines[spine].set_visible(False)
    ax.tick_params(axis="x", labelsize=7)
    ax.set_xlabel("Substitutions per site", fontsize=8)
    fig.suptitle("Recombination-filtered phylogeny of 57 Pantoea piersonii genomes", fontsize=12, y=.985)
    handles = [Patch(facecolor=c, label=k) for k, c in SOURCE_COLORS.items()]
    handles += [Patch(facecolor=c, label=k) for k, c in YEAR_COLORS.items() if k]
    handles += [Patch(facecolor="white", edgecolor="#aaa", label="Year unavailable / not applicable")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .035), ncol=4, frameon=False, fontsize=8)
    fig.text(.5, .012, "Gene matrix: blue = detected; light grey = not detected. Orange nodes: SH-aLRT ≥80 and UFBoot ≥95.", ha="center", fontsize=8)
    args.outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.outdir / "Figure_1.png", dpi=300, bbox_inches="tight")
    fig.savefig(args.outdir / "Figure_1.pdf", bbox_inches="tight")


if __name__ == "__main__":
    main()
