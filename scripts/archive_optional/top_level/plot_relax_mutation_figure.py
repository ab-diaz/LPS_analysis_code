#!/usr/bin/env python3
"""Plot RELAX k values and ISS mutation density for rfb/waaL genes."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

from Bio import SeqIO


GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL"]
VALID_BASES = set("ACGT")
METADATA_ALIASES = {
    "Pantoea_piersonii_F8_6S_D2_EB4_spades": "Pantoea_piersonii_F9_6S_D2_EB4_spades",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase2-dir", default="lps_analysis/selection_phase2")
    parser.add_argument("--metadata", default="supplTbl4.csv")
    parser.add_argument("--out-prefix", default="lps_analysis/selection_phase2/relax_k_mutation_heatmap")
    return parser.parse_args()


def clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).replace("\xa0", " ").strip() for c in df.columns]
    return df


def metadata_key(isolate: str) -> str:
    return METADATA_ALIASES.get(isolate, isolate)


def display_isolate(isolate: str) -> str:
    return metadata_key(isolate)


def extract_year(value: object) -> float:
    if pd.isna(value):
        return np.nan
    match = re.search(r"(19|20)\d{2}", str(value))
    return float(match.group(0)) if match else np.nan


def shorten_isolate(isolate: str) -> str:
    isolate = display_isolate(isolate)
    prefix = "Pantoea_piersonii_"
    return isolate[len(prefix) :] if isolate.startswith(prefix) else isolate


def p_label(p_value: float) -> str:
    if not np.isfinite(p_value) or p_value >= 0.05:
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    return "*"


def format_p_value(p_value: float) -> str:
    if not np.isfinite(p_value):
        return "p = NA"
    if p_value == 0 or p_value < 0.001:
        return "p < 0.001"
    return f"p = {p_value:.3f}"


def count_mismatches(seq: str, founder: str) -> tuple[int, int, float]:
    mismatches = 0
    comparable = 0
    for base, ref_base in zip(seq.upper(), founder.upper()):
        if base in VALID_BASES and ref_base in VALID_BASES:
            comparable += 1
            mismatches += base != ref_base
    density = (mismatches / comparable * 1000.0) if comparable else np.nan
    return int(mismatches), int(comparable), float(density)


def load_metadata(path: Path) -> pd.DataFrame:
    metadata = clean_columns(pd.read_csv(path, sep="\t"))
    metadata["year"] = metadata["Collection date"].map(extract_year)
    return metadata[["Sample ID", "Collection date", "year", "Category"]]


def build_mutation_table(phase2_dir: Path, metadata_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    id_map = pd.read_csv(phase2_dir / "sequence_id_map.tsv", sep="\t")
    metadata = load_metadata(metadata_path)
    meta_by_isolate = metadata.set_index("Sample ID", drop=False)

    rows: list[dict[str, object]] = []
    founder_rows: list[dict[str, object]] = []

    for gene in GENES:
        gene_map = id_map[(id_map["gene"] == gene) & (id_map["group"] == "ISS")].copy()
        gene_map["metadata_key"] = gene_map["isolate"].map(metadata_key)
        gene_map["display_isolate"] = gene_map["isolate"].map(display_isolate)
        gene_map["year"] = gene_map["metadata_key"].map(meta_by_isolate["year"])
        gene_map["collection_date"] = gene_map["metadata_key"].map(meta_by_isolate["Collection date"])

        records = {rec.id: str(rec.seq) for rec in SeqIO.parse(phase2_dir / f"{gene}.codon.aln.fasta", "fasta")}
        available = gene_map[gene_map["safe_id"].isin(records)].copy()
        founder_candidates = available.dropna(subset=["year"]).sort_values(["year", "display_isolate", "safe_id"])
        if founder_candidates.empty:
            founder_candidates = available.sort_values(["display_isolate", "safe_id"])
        if founder_candidates.empty:
            continue

        founder = founder_candidates.iloc[0]
        founder_seq = records[founder["safe_id"]]
        founder_rows.append(
            {
                "gene": gene,
                "founder_safe_id": founder["safe_id"],
                "founder_isolate": founder["display_isolate"],
                "source_isolate_id": founder["isolate"],
                "founder_year": founder["year"],
                "founder_collection_date": founder["collection_date"],
                "alignment_bases": len(founder_seq),
            }
        )

        for _, row in available.iterrows():
            mismatches, comparable, density = count_mismatches(records[row["safe_id"]], founder_seq)
            rows.append(
                {
                    "gene": gene,
                    "safe_id": row["safe_id"],
                    "isolate": row["display_isolate"],
                    "source_isolate_id": row["isolate"],
                    "collection_date": row["collection_date"],
                    "year": row["year"],
                    "founder_isolate": founder["display_isolate"],
                    "mismatches": mismatches,
                    "comparable_bases": comparable,
                    "mismatches_per_kb": density,
                }
            )

    mutation = pd.DataFrame(rows)
    founders = pd.DataFrame(founder_rows)

    unique_isolates = (
        id_map[id_map["group"] == "ISS"][["isolate"]]
        .drop_duplicates()
        .assign(
            metadata_key=lambda d: d["isolate"].map(metadata_key),
            display_isolate=lambda d: d["isolate"].map(display_isolate),
        )
        .assign(
            year=lambda d: d["metadata_key"].map(meta_by_isolate["year"]),
            collection_date=lambda d: d["metadata_key"].map(meta_by_isolate["Collection date"]),
            metadata_matched=lambda d: d["metadata_key"].isin(meta_by_isolate.index),
        )
        .sort_values(["year", "display_isolate"], na_position="last")
    )
    return mutation, founders, unique_isolates


def build_position_table(phase2_dir: Path, mutation: pd.DataFrame, founders: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    founder_by_gene = founders.set_index("gene")

    for gene in GENES:
        records = {rec.id: str(rec.seq) for rec in SeqIO.parse(phase2_dir / f"{gene}.codon.aln.fasta", "fasta")}
        if gene not in founder_by_gene.index:
            continue
        founder_id = founder_by_gene.loc[gene, "founder_safe_id"]
        if founder_id not in records:
            continue
        founder_seq = records[founder_id].upper()

        for _, row in mutation[mutation["gene"] == gene].iterrows():
            safe_id = row["safe_id"]
            if safe_id not in records:
                continue
            seq = records[safe_id].upper()
            for idx, (base, ref_base) in enumerate(zip(seq, founder_seq), start=1):
                if base in VALID_BASES and ref_base in VALID_BASES and base != ref_base:
                    rows.append(
                        {
                            "gene": gene,
                            "aligned_position": idx,
                            "isolate": row["isolate"],
                            "source_isolate_id": row["source_isolate_id"],
                            "safe_id": safe_id,
                            "year": row["year"],
                            "collection_date": row["collection_date"],
                            "founder_isolate": row["founder_isolate"],
                            "founder_base": ref_base,
                            "isolate_base": base,
                            "substitution": f"{ref_base}>{base}",
                        }
                    )

    return pd.DataFrame(rows)


def plot_figure(relax: pd.DataFrame, mutation: pd.DataFrame, isolates: pd.DataFrame, out_prefix: Path) -> None:
    relax = relax.set_index("gene").loc[GENES].reset_index()
    pivot = mutation.pivot_table(index="isolate", columns="gene", values="mismatches_per_kb", aggfunc="first")
    pivot = pivot.reindex(index=isolates["display_isolate"], columns=GENES)
    heat = pivot.to_numpy(dtype=float)
    masked_heat = np.ma.masked_invalid(heat)

    plt.rcParams.update(
        {
            "axes.titlesize": 13,
            "axes.labelsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 8,
            "font.size": 9,
        }
    )
    fig = plt.figure(figsize=(9.6, 12.0), constrained_layout=True)
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1.0, 3.0])

    ax_bar = fig.add_subplot(gs[0])
    colors = ["#4f5b66" if p >= 0.05 else "#b23a48" for p in relax["p_value"]]
    log_floor = 0.65
    ax_bar.bar(
        relax["gene"],
        relax["k"] - log_floor,
        bottom=log_floor,
        color=colors,
        edgecolor="#1f252b",
        linewidth=0.8,
    )
    ax_bar.set_yscale("log")
    ax_bar.axhline(1.0, color="#222222", linewidth=1.0, linestyle="--")
    ax_bar.text(4.48, 1.03, "k = 1", ha="right", va="bottom", fontsize=8, color="#222222")
    ax_bar.set_ylim(log_floor, 45)
    ax_bar.set_yticks([0.7, 1, 2, 5, 10, 20, 40])
    ax_bar.set_yticklabels(["0.7", "1", "2", "5", "10", "20", "40"])
    ax_bar.set_ylabel("RELAX selection intensity k (log scale)")
    ax_bar.set_title("A. RELAX detects intensified selection for rfbD and waaL", loc="left", fontweight="bold")
    ax_bar.spines[["top", "right"]].set_visible(False)
    for i, row in relax.iterrows():
        label = p_label(float(row["p_value"]))
        text = f"k = {row['k']:.2f}\n{format_p_value(float(row['p_value']))}"
        if label:
            text = f"{label}\n{text}"
            y_text = row["k"] / 1.10
            va = "top"
            weight = "bold"
            color = "white"
        else:
            y_text = row["k"] * 1.12
            va = "bottom"
            weight = "normal"
            color = "black"
        ax_bar.text(
            i,
            y_text,
            text,
            ha="center",
            va=va,
            fontsize=8,
            fontweight=weight,
            color=color,
            linespacing=1.05,
        )
    legend_handles = [
        Patch(facecolor="#b23a48", edgecolor="#1f252b", label="Significant RELAX shift"),
        Patch(facecolor="#4f5b66", edgecolor="#1f252b", label="Not significant"),
    ]
    ax_bar.legend(handles=legend_handles, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.98), fontsize=8)

    ax_heat = fig.add_subplot(gs[1])
    cmap = LinearSegmentedColormap.from_list("mutation_density", ["#f7f7f7", "#93b7be", "#2f6690", "#7d1d3f"])
    cmap.set_bad("#d8d8d8")
    image = ax_heat.imshow(masked_heat, aspect="auto", interpolation="nearest", cmap=cmap)
    ax_heat.set_xticks(np.arange(len(GENES)))
    ax_heat.set_xticklabels(GENES)
    labels = []
    for _, row in isolates.iterrows():
        year = "" if pd.isna(row["year"]) else f"{int(row['year'])} "
        labels.append(f"{year}{shorten_isolate(row['display_isolate'])}")
    ax_heat.set_yticks(np.arange(len(labels)))
    ax_heat.set_yticklabels(labels, fontsize=6)
    ax_heat.set_title("B. Founder-relative nucleotide divergence in ISS isolates", loc="left", fontweight="bold")
    ax_heat.tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False)
    ax_heat.set_xticks(np.arange(-0.5, len(GENES), 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, len(labels), 1), minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=0.6)
    ax_heat.tick_params(which="minor", bottom=False, left=False)
    ax_heat.spines[:].set_visible(False)

    years = isolates["year"].to_numpy()
    for idx in range(1, len(years)):
        if not pd.isna(years[idx]) and not pd.isna(years[idx - 1]) and years[idx] != years[idx - 1]:
            ax_heat.axhline(idx - 0.5, color="#333333", linewidth=0.8)

    cbar = fig.colorbar(image, ax=ax_heat, fraction=0.025, pad=0.02)
    cbar.set_label("Mismatches/kb vs founder")

    fig.savefig(out_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(out_prefix.with_suffix(".pdf"))
    plt.close(fig)


def plot_longitudinal_map(
    phase2_dir: Path,
    mutation: pd.DataFrame,
    position_events: pd.DataFrame,
    isolates: pd.DataFrame,
    out_prefix: Path,
) -> None:
    gene_lengths = {}
    for gene in GENES:
        first_record = next(SeqIO.parse(phase2_dir / f"{gene}.codon.aln.fasta", "fasta"))
        gene_lengths[gene] = len(first_record.seq)

    isolate_order = list(isolates["display_isolate"])
    y_lookup = {isolate: idx for idx, isolate in enumerate(isolate_order)}

    fig = plt.figure(figsize=(13.5, 12.0), constrained_layout=True)
    gs = fig.add_gridspec(1, len(GENES), width_ratios=[gene_lengths[g] for g in GENES], wspace=0.04)
    axes = [fig.add_subplot(gs[0, i]) for i in range(len(GENES))]

    for ax, gene in zip(axes, GENES):
        ax.set_xlim(0.5, gene_lengths[gene] + 0.5)
        ax.set_ylim(len(isolate_order) - 0.5, -0.5)
        ax.set_title(gene, fontweight="bold")
        ax.set_xlabel("Aligned nt")
        ax.spines[["top", "right"]].set_visible(False)

        present = set(mutation.loc[mutation["gene"] == gene, "isolate"])
        for isolate in isolate_order:
            y = y_lookup[isolate]
            if isolate not in present:
                ax.axhspan(y - 0.5, y + 0.5, color="#e3e3e3", zorder=0)

        gene_events = position_events[position_events["gene"] == gene].copy()
        if not gene_events.empty:
            gene_events["y"] = gene_events["isolate"].map(y_lookup)
            ax.scatter(
                gene_events["aligned_position"],
                gene_events["y"],
                marker="|",
                s=34,
                linewidths=1.1,
                color="#b23a48",
                alpha=0.9,
                zorder=3,
            )

        years = isolates["year"].to_numpy()
        for idx in range(1, len(years)):
            if not pd.isna(years[idx]) and not pd.isna(years[idx - 1]) and years[idx] != years[idx - 1]:
                ax.axhline(idx - 0.5, color="#333333", linewidth=0.8, zorder=2)

        ax.grid(axis="x", color="#ededed", linewidth=0.7)
        ax.set_xticks([1, gene_lengths[gene]])
        ax.set_xticklabels(["1", str(gene_lengths[gene])], fontsize=7)

    labels = []
    for _, row in isolates.iterrows():
        year = "" if pd.isna(row["year"]) else f"{int(row['year'])} "
        labels.append(f"{year}{shorten_isolate(row['display_isolate'])}")
    axes[0].set_yticks(np.arange(len(labels)))
    axes[0].set_yticklabels(labels, fontsize=6)
    axes[0].set_ylabel("ISS isolate ordered by collection date")
    for ax in axes[1:]:
        ax.set_yticks(np.arange(len(labels)))
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)

    legend_handles = [
        plt.Line2D([0], [0], marker="|", color="#b23a48", linestyle="None", markersize=11, markeredgewidth=1.3, label="Nucleotide mismatch"),
        Patch(facecolor="#e3e3e3", edgecolor="none", label="Gene absent from curated set"),
    ]
    fig.legend(handles=legend_handles, loc="upper right", bbox_to_anchor=(0.995, 0.995), frameon=False)
    fig.suptitle(
        "Aligned gene-coordinate mutation map relative to earliest available ISS sequence",
        x=0.01,
        ha="left",
        fontweight="bold",
    )
    fig.savefig(out_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(out_prefix.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    phase2_dir = Path(args.phase2_dir)
    metadata_path = Path(args.metadata)
    out_prefix = Path(args.out_prefix)

    mutation, founders, isolates = build_mutation_table(phase2_dir, metadata_path)
    relax = pd.read_csv(phase2_dir / "relax_results_summary.tsv", sep="\t")

    mutation.to_csv(phase2_dir / "mutation_density_iss.tsv", sep="\t", index=False)
    founders.to_csv(phase2_dir / "mutation_density_founders.tsv", sep="\t", index=False)
    isolates.to_csv(phase2_dir / "mutation_density_isolate_order.tsv", sep="\t", index=False)
    position_events = build_position_table(phase2_dir, mutation, founders)
    position_events.to_csv(phase2_dir / "longitudinal_mutation_events.tsv", sep="\t", index=False)
    plot_figure(relax, mutation, isolates, out_prefix)
    plot_longitudinal_map(phase2_dir, mutation, position_events, isolates, phase2_dir / "longitudinal_mutation_map")

    unmatched = isolates[~isolates["metadata_matched"]]
    print(f"Wrote {out_prefix.with_suffix('.png')}")
    print(f"Wrote {out_prefix.with_suffix('.pdf')}")
    print(f"Wrote {phase2_dir / 'longitudinal_mutation_map.png'}")
    print(f"Wrote {phase2_dir / 'longitudinal_mutation_map.pdf'}")
    print(f"Wrote {phase2_dir / 'mutation_density_iss.tsv'}")
    print(f"Wrote {phase2_dir / 'longitudinal_mutation_events.tsv'}")
    if not unmatched.empty:
        print("ISS isolates without SupplTbl4 metadata:")
        for isolate in unmatched["isolate"]:
            print(f"  {isolate}")


if __name__ == "__main__":
    main()
