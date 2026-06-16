#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow
import numpy as np
import pandas as pd


BASE = Path("results/lps_analysis")
PHASE = BASE / "selection_phase2"
OUT = PHASE / "pantoea_route1_strengthening"
FIG = PHASE / "figures"
RESULTS = Path("results")

STRICT_LONG = BASE / "lps_strict_presence_long.tsv"
STRICT_HITS = BASE / "lps_strict_hits.tsv"
ISS_ORDER = PHASE / "mutation_density_isolate_order.tsv"

TARGET_GENES = ["rfbA", "rfbB", "rfbC", "rfbD", "waaL", "wzm", "wzt", "rfaZ"]
TARGET_LOWER = [g.lower() for g in TARGET_GENES]
PRIMARY_OPERON = ["rfba", "rfbb", "rfbc", "rfbd"]


def fisher_two_sided(a, b, c, d):
    n = a + b + c + d
    row1 = a + b
    col1 = a + c

    def hypergeom(x):
        if x < max(0, row1 - (n - col1)) or x > min(row1, col1):
            return 0.0
        return math.comb(col1, x) * math.comb(n - col1, row1 - x) / math.comb(n, row1)

    obs = hypergeom(a)
    lo = max(0, row1 - (n - col1))
    hi = min(row1, col1)
    return min(1.0, sum(hypergeom(x) for x in range(lo, hi + 1) if hypergeom(x) <= obs + 1e-15))


def bh_qvalues(pvals):
    pvals = np.asarray(pvals, dtype=float)
    order = np.argsort(pvals)
    q = np.empty_like(pvals)
    prev = 1.0
    m = len(pvals)
    for rank, idx in enumerate(order[::-1], start=1):
        true_rank = m - rank + 1
        val = min(prev, pvals[idx] * m / true_rank)
        q[idx] = val
        prev = val
    return q


def odds_ratio(a, b, c, d):
    if b == 0 or c == 0:
        if a > 0 and d > 0:
            return math.inf
        return 0.0
    return (a * d) / (b * c)


def normalize_gene(g):
    if pd.isna(g):
        return ""
    return str(g).strip().lower()


def short_species(isolate):
    parts = isolate.split("_")
    return "_".join(parts[:2]) if len(parts) >= 2 else isolate


def infer_iss_year(isolate, year_map):
    if isolate in year_map:
        return int(year_map[isolate])
    if "IIIF" in isolate:
        return 2016
    if "_F3" in isolate or "_F8" in isolate:
        return 2018
    if "_F9" in isolate or "_F10" in isolate:
        return 2021
    if "_F11" in isolate:
        return 2022
    return None


def load_presence():
    long = pd.read_csv(STRICT_LONG, sep="\t")
    long["gene_id"] = long["gene_id"].str.lower()
    mat = (
        long.pivot_table(index=["isolate", "group"], columns="gene_id", values="present", fill_value=0)
        .reset_index()
    )
    return long, mat


def exact_profile_clusters(mat, genes):
    cols = [g for g in genes if g in mat.columns]
    rows = []
    for _, row in mat.iterrows():
        sig = "".join(str(int(row[g])) for g in cols)
        profile = hashlib.sha1(sig.encode()).hexdigest()[:10]
        rows.append(
            {
                "isolate": row["isolate"],
                "group": row["group"],
                "profile_cluster": profile,
                "profile_signature": sig,
                "profile_genes": ",".join(cols),
            }
        )
    return pd.DataFrame(rows)


def cluster_level_tests(mat, clusters, genes):
    merged = mat.merge(clusters[["isolate", "profile_cluster"]], on="isolate")
    rows = []
    for gene in genes:
        g = gene.lower()
        if g not in merged.columns:
            continue
        ctab = (
            merged.groupby(["group", "profile_cluster"])[g]
            .max()
            .reset_index()
            .groupby("group")[g]
            .agg(["sum", "count"])
        )
        iss_present = int(ctab.loc["ISS", "sum"]) if "ISS" in ctab.index else 0
        iss_total = int(ctab.loc["ISS", "count"]) if "ISS" in ctab.index else 0
        earth_present = int(ctab.loc["Earth", "sum"]) if "Earth" in ctab.index else 0
        earth_total = int(ctab.loc["Earth", "count"]) if "Earth" in ctab.index else 0
        a = iss_present
        b = iss_total - iss_present
        c = earth_present
        d = earth_total - earth_present
        rows.append(
            {
                "gene_id": gene,
                "iss_profile_clusters_present": a,
                "iss_profile_clusters_absent": b,
                "earth_profile_clusters_present": c,
                "earth_profile_clusters_absent": d,
                "odds_ratio": odds_ratio(a, b, c, d),
                "p_value": fisher_two_sided(a, b, c, d),
            }
        )
    out = pd.DataFrame(rows)
    out["q_value"] = bh_qvalues(out["p_value"]) if len(out) else []
    return out


def year_persistence(mat, genes):
    order = pd.read_csv(ISS_ORDER, sep="\t")
    year_map = dict(zip(order["isolate"], order["year"]))
    rows = []
    iss = mat[mat["group"] == "ISS"].copy()
    iss["year"] = iss["isolate"].apply(lambda x: infer_iss_year(x, year_map))
    for gene in genes:
        g = gene.lower()
        if g not in iss.columns:
            continue
        for year, sub in iss.groupby("year", dropna=False):
            rows.append(
                {
                    "gene_id": gene,
                    "year": int(year) if pd.notna(year) else "",
                    "iss_isolates": len(sub),
                    "present": int(sub[g].sum()),
                    "fraction_present": float(sub[g].mean()) if len(sub) else 0,
                }
            )
    return pd.DataFrame(rows)


def leave_year_out_tests(mat, genes):
    order = pd.read_csv(ISS_ORDER, sep="\t")
    year_map = dict(zip(order["isolate"], order["year"]))
    m = mat.copy()
    m["year"] = m.apply(
        lambda r: infer_iss_year(r["isolate"], year_map) if r["group"] == "ISS" else None,
        axis=1,
    )
    rows = []
    for excluded_year in sorted(y for y in m["year"].dropna().unique()):
        sub = m[(m["group"] != "ISS") | (m["year"] != excluded_year)]
        for gene in genes:
            g = gene.lower()
            if g not in sub.columns:
                continue
            iss = sub[sub["group"] == "ISS"]
            earth = sub[sub["group"] == "Earth"]
            a = int(iss[g].sum())
            b = len(iss) - a
            c = int(earth[g].sum())
            d = len(earth) - c
            rows.append(
                {
                    "excluded_iss_year": int(excluded_year),
                    "gene_id": gene,
                    "iss_present": a,
                    "iss_absent": b,
                    "earth_present": c,
                    "earth_absent": d,
                    "odds_ratio": odds_ratio(a, b, c, d),
                    "p_value": fisher_two_sided(a, b, c, d),
                }
            )
    out = pd.DataFrame(rows)
    if len(out):
        out["q_value_within_all_leave_year_tests"] = bh_qvalues(out["p_value"])
    return out


def read_bakta_tsv(isolate):
    path = RESULTS / isolate / "bakta_output" / f"{isolate}.tsv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", comment="#", names=[
        "contig",
        "type",
        "start",
        "end",
        "strand",
        "locus_tag",
        "gene",
        "product",
        "dbxrefs",
    ])


def build_locus_tables():
    hits = pd.read_csv(STRICT_HITS, sep="\t")
    hits["gene_norm"] = hits["preferred_name"].map(normalize_gene)
    hits = hits[hits["gene_norm"].isin(TARGET_LOWER)].copy()

    records = []
    neighborhoods = []
    for isolate, sub in hits.groupby("isolate"):
        bakta = read_bakta_tsv(isolate)
        if bakta.empty:
            continue
        cds = bakta[bakta["type"].eq("cds")].copy().reset_index(drop=True)
        cds["gene_norm"] = cds["gene"].fillna("").map(normalize_gene)
        by_locus = {str(r.locus_tag): i for i, r in cds.iterrows()}
        for _, h in sub.iterrows():
            idx = by_locus.get(str(h["query"]))
            if idx is None:
                # Some eggNOG query IDs are not exact Bakta locus tags.
                candidates = cds[cds["gene_norm"].eq(h["gene_norm"])]
                if candidates.empty:
                    continue
                idx = int(candidates.index[0])
            r = cds.loc[idx]
            records.append(
                {
                    "isolate": isolate,
                    "group": "ISS" if isolate in set(pd.read_csv(STRICT_LONG, sep="\t").query("group == 'ISS'")["isolate"]) else "Earth",
                    "gene": h["gene_norm"],
                    "locus_tag": r["locus_tag"],
                    "contig": r["contig"],
                    "start": int(r["start"]),
                    "end": int(r["end"]),
                    "strand": r["strand"],
                    "product": r["product"],
                    "eggnog_description": h["description"],
                }
            )
            lo, hi = max(0, idx - 6), min(len(cds), idx + 7)
            for j in range(lo, hi):
                nr = cds.loc[j]
                neighborhoods.append(
                    {
                        "anchor_gene": h["gene_norm"],
                        "anchor_locus": r["locus_tag"],
                        "isolate": isolate,
                        "neighbor_index": j - idx,
                        "contig": nr["contig"],
                        "start": int(nr["start"]),
                        "end": int(nr["end"]),
                        "strand": nr["strand"],
                        "locus_tag": nr["locus_tag"],
                        "gene": normalize_gene(nr["gene"]),
                        "product": nr["product"],
                    }
                )

    loc = pd.DataFrame(records).drop_duplicates()
    neigh = pd.DataFrame(neighborhoods).drop_duplicates()
    return loc, neigh


def summarize_loci(loc):
    rows = []
    for isolate, sub in loc.groupby("isolate"):
        group = sub["group"].iloc[0]
        genes = set(sub["gene"])
        rfb = sub[sub["gene"].isin(PRIMARY_OPERON)]
        complete_rfb = all(g in genes for g in PRIMARY_OPERON)
        rfb_colocalized = False
        rfb_span = ""
        rfb_contig = ""
        rfb_order = ""
        if complete_rfb:
            for contig, cs in rfb.groupby("contig"):
                if set(cs["gene"]) >= set(PRIMARY_OPERON):
                    starts = cs.drop_duplicates("gene").sort_values("start")
                    span = int(starts["end"].max() - starts["start"].min() + 1)
                    if span <= 15000:
                        rfb_colocalized = True
                        rfb_span = span
                        rfb_contig = contig
                        rfb_order = "-".join(starts["gene"])
                        break
        wzm_wzt_pair = False
        if {"wzm", "wzt"} <= genes:
            wz = sub[sub["gene"].isin(["wzm", "wzt"])]
            for _, cs in wz.groupby("contig"):
                if set(cs["gene"]) == {"wzm", "wzt"}:
                    if int(cs["end"].max() - cs["start"].min() + 1) <= 10000:
                        wzm_wzt_pair = True
                        break
        rows.append(
            {
                "isolate": isolate,
                "group": group,
                "genes_detected": ",".join(sorted(genes)),
                "complete_rfbABCD": int(complete_rfb),
                "rfbABCD_colocalized_15kb": int(rfb_colocalized),
                "rfb_contig": rfb_contig,
                "rfb_order_by_coordinate": rfb_order,
                "rfb_span_bp": rfb_span,
                "has_waaL": int("waal" in genes),
                "has_wzm": int("wzm" in genes),
                "has_wzt": int("wzt" in genes),
                "wzm_wzt_colocalized_10kb": int(wzm_wzt_pair),
                "has_rfaZ": int("rfaz" in genes),
            }
        )
    return pd.DataFrame(rows)


def representative_isolates(summary):
    candidates = []
    for iso in [
        "Pantoea_piersonii_F10_5S-D1_P5",
        "Pantoea_piersonii_F10_5S_D1_EP3",
        "Pantoea_piersonii_F11_7S_D1_EB1",
    ]:
        if iso in set(summary["isolate"]):
            candidates.append(iso)
            break
    earth_with = summary[
        (summary["group"] == "Earth")
        & (summary["rfbABCD_colocalized_15kb"] == 1)
    ]["isolate"].tolist()
    earth_without = summary[
        (summary["group"] == "Earth")
        & (summary["complete_rfbABCD"] == 0)
    ]["isolate"].tolist()
    if earth_with:
        candidates.append(earth_with[0])
    if earth_without:
        candidates.append(earth_without[0])
    return candidates[:3]


def plot_locus_architecture(loc, summary):
    reps = representative_isolates(summary)
    if not reps:
        return
    colors = {
        "rfba": "#4c78a8",
        "rfbb": "#f58518",
        "rfbc": "#54a24b",
        "rfbd": "#b279a2",
        "waal": "#e45756",
        "wzm": "#72b7b2",
        "wzt": "#ff9da6",
        "rfaz": "#9d755d",
    }
    fig, axes = plt.subplots(len(reps), 1, figsize=(9, 1.45 * len(reps)), sharex=False)
    if len(reps) == 1:
        axes = [axes]
    for ax, iso in zip(axes, reps):
        sub = loc[loc["isolate"].eq(iso)].copy()
        rfb_row = summary[summary["isolate"].eq(iso)]
        if not rfb_row.empty and rfb_row.iloc[0]["rfb_contig"]:
            contig = rfb_row.iloc[0]["rfb_contig"]
            sub = sub[sub["contig"].eq(contig)]
        elif not sub.empty:
            contig = sub["contig"].mode().iloc[0]
            sub = sub[sub["contig"].eq(contig)]
        if sub.empty:
            continue
        start0 = int(sub["start"].min()) - 1500
        xmax = int(sub["end"].max()) - start0 + 1500
        for _, r in sub.sort_values("start").iterrows():
            x = int(r["start"]) - start0
            width = max(250, int(r["end"]) - int(r["start"]))
            dx = width if r["strand"] == "+" else -width
            x0 = x if r["strand"] == "+" else x + width
            ax.add_patch(
                FancyArrow(
                    x0,
                    0,
                    dx,
                    0,
                    width=0.22,
                    length_includes_head=True,
                    head_width=0.45,
                    head_length=min(250, width * 0.35),
                    color=colors.get(r["gene"], "#777777"),
                    edgecolor="#222222",
                )
            )
            ax.text(x + width / 2, 0.34, r["gene"], ha="center", va="bottom", fontsize=8)
        ax.set_ylim(-0.6, 0.8)
        ax.set_xlim(0, xmax)
        ax.set_yticks([])
        ax.set_title(iso.replace("Pantoea_piersonii_", ""), loc="left", fontsize=9)
        ax.spines[["left", "right", "top"]].set_visible(False)
        ax.set_xlabel("relative position within detected O-antigen/LPS locus (bp)")
    fig.suptitle("Representative O-antigen/LPS target-gene architecture", x=0.01, ha="left", fontsize=12, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"figure_oantigen_locus_architecture_candidate.{ext}", dpi=300 if ext == "png" else None)


def plot_robustness(cluster_tests, year_df, locus_summary):
    genes = [g for g in TARGET_GENES if g in set(cluster_tests["gene_id"])]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))

    ct = cluster_tests.set_index("gene_id").reindex(genes).reset_index()
    x = np.arange(len(ct))
    iss_frac = ct["iss_profile_clusters_present"] / (
        ct["iss_profile_clusters_present"] + ct["iss_profile_clusters_absent"]
    )
    earth_frac = ct["earth_profile_clusters_present"] / (
        ct["earth_profile_clusters_present"] + ct["earth_profile_clusters_absent"]
    )
    axes[0].bar(x - 0.18, iss_frac, width=0.36, label="ISS profile clusters", color="#b23a48")
    axes[0].bar(x + 0.18, earth_frac, width=0.36, label="Earth profile clusters", color="#3d6f9f")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(genes, rotation=45, ha="right")
    axes[0].set_ylabel("Fraction of exact LPS-profile clusters present")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_title("A. Cluster-collapsed gene presence", loc="left", weight="bold")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].spines[["top", "right"]].set_visible(False)

    subset = year_df[year_df["gene_id"].isin(["rfbA", "rfbB", "rfbC", "rfbD", "waaL", "wzm"])].copy()
    pivot = subset.pivot_table(index="gene_id", columns="year", values="fraction_present").reindex(["rfbA", "rfbB", "rfbC", "rfbD", "waaL", "wzm"])
    im = axes[1].imshow(pivot.fillna(0).to_numpy(), aspect="auto", vmin=0, vmax=1, cmap="Reds")
    axes[1].set_yticks(np.arange(len(pivot.index)))
    axes[1].set_yticklabels(pivot.index)
    axes[1].set_xticks(np.arange(len(pivot.columns)))
    axes[1].set_xticklabels([str(int(c)) for c in pivot.columns], rotation=0)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            v = pivot.iloc[i, j]
            axes[1].text(j, i, "" if pd.isna(v) else f"{v:.2f}", ha="center", va="center", fontsize=8)
    axes[1].set_title("B. ISS year-level persistence", loc="left", weight="bold")
    cbar = fig.colorbar(im, ax=axes[1], shrink=0.82)
    cbar.set_label("fraction present")

    fig.tight_layout()
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"figure6_founder_year_persistence.{ext}", dpi=300 if ext == "png" else None)


def main():
    OUT.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    _, mat = load_presence()

    strict_genes = sorted([c for c in mat.columns if c not in {"isolate", "group"}])
    clusters_all = exact_profile_clusters(mat, strict_genes)
    clusters_target = exact_profile_clusters(mat, TARGET_LOWER)
    clusters_all.to_csv(OUT / "exact_strict_lps_profile_clusters.tsv", sep="\t", index=False)
    clusters_target.to_csv(OUT / "exact_target_gene_profile_clusters.tsv", sep="\t", index=False)

    cluster_tests_all = cluster_level_tests(mat, clusters_all, TARGET_GENES)
    cluster_tests_target = cluster_level_tests(mat, clusters_target, TARGET_GENES)
    cluster_tests_all.to_csv(OUT / "cluster_collapsed_strict_lps_fisher.tsv", sep="\t", index=False)
    cluster_tests_target.to_csv(OUT / "cluster_collapsed_target_profile_fisher.tsv", sep="\t", index=False)

    year_df = year_persistence(mat, TARGET_GENES)
    year_df.to_csv(OUT / "iss_year_persistence_target_genes.tsv", sep="\t", index=False)

    leave_year = leave_year_out_tests(mat, TARGET_GENES)
    leave_year.to_csv(OUT / "leave_one_iss_year_out_fisher.tsv", sep="\t", index=False)

    loc, neigh = build_locus_tables()
    loc.to_csv(OUT / "target_gene_bakta_locations.tsv", sep="\t", index=False)
    neigh.to_csv(OUT / "target_gene_bakta_neighborhoods_pm6.tsv", sep="\t", index=False)
    locus_summary = summarize_loci(loc)
    locus_summary.to_csv(OUT / "oantigen_locus_architecture_summary.tsv", sep="\t", index=False)

    plot_locus_architecture(loc, locus_summary)
    plot_robustness(cluster_tests_all, year_df, locus_summary)

    with open(OUT / "route1_strengthening_summary.txt", "w", encoding="utf-8") as handle:
        handle.write("Route 1 Pantoea-focused strengthening outputs\\n")
        handle.write(f"Total isolates in strict matrix: {len(mat)}\\n")
        handle.write(f"ISS isolates: {(mat['group'] == 'ISS').sum()}\\n")
        handle.write(f"Earth isolates: {(mat['group'] == 'Earth').sum()}\\n")
        handle.write(f"Exact strict LPS profile clusters: {clusters_all['profile_cluster'].nunique()}\\n")
        handle.write(
            f"ISS exact strict LPS profile clusters: {clusters_all[clusters_all['group']=='ISS']['profile_cluster'].nunique()}\\n"
        )
        handle.write(
            f"Earth exact strict LPS profile clusters: {clusters_all[clusters_all['group']=='Earth']['profile_cluster'].nunique()}\\n"
        )
        handle.write("\\nCluster-collapsed strict-profile Fisher tests:\\n")
        handle.write(cluster_tests_all.to_string(index=False))
        handle.write("\\n\\nISS year persistence:\\n")
        handle.write(year_df.to_string(index=False))
        handle.write("\\n\\nLocus architecture group summary:\\n")
        if len(locus_summary):
            handle.write(
                locus_summary.groupby("group")[
                    [
                        "complete_rfbABCD",
                        "rfbABCD_colocalized_15kb",
                        "has_waaL",
                        "has_wzm",
                        "has_wzt",
                        "wzm_wzt_colocalized_10kb",
                        "has_rfaZ",
                    ]
                ]
                .agg(["sum", "count"])
                .to_string()
            )

    print(OUT)
    print(FIG / "figure6_founder_year_persistence.png")
    print(FIG / "figure_oantigen_locus_architecture_candidate.png")


if __name__ == "__main__":
    main()
