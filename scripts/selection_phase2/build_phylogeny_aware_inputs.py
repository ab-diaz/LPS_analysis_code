#!/usr/bin/env python3
"""Build clean inputs for phylogeny-aware O-antigen/LPS association tests."""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "phylogeny_aware_inputs"
PROJECT_RESULTS = ROOT.parents[1]

TARGET_PROFILE = ROOT / "pantoea_route1_strengthening" / "exact_target_gene_profile_clusters.tsv"
STRICT_PROFILE = ROOT / "pantoea_route1_strengthening" / "exact_strict_lps_profile_clusters.tsv"
LOCUS_SUMMARY = ROOT / "pantoea_route1_strengthening" / "oantigen_locus_architecture_summary.tsv"
ISS_YEAR_METADATA = ROOT / "mutation_density_isolate_order.tsv"
DIAA_FFN = ROOT.parent / "lps_shared_genes" / "diaA.ffn"
DIAA_FAA = ROOT.parent / "lps_shared_genes" / "diaA.faa"
SEQUENCE_ID_MAP = ROOT / "sequence_id_map.tsv"

TARGET_GENES = ["rfba", "rfbb", "rfbc", "rfbd", "waal", "wzm", "wzt", "rfaz"]
OUT_GENE_NAMES = {
    "rfba": "rfbA",
    "rfbb": "rfbB",
    "rfbc": "rfbC",
    "rfbd": "rfbD",
    "waal": "waaL",
    "wzm": "wzm",
    "wzt": "wzt",
    "rfaz": "rfaZ",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def species_from_isolate(isolate: str) -> str:
    parts = isolate.split("_")
    if len(parts) >= 2 and parts[1] != "sp":
        return "_".join(parts[:2])
    if len(parts) >= 3 and parts[1] == "sp":
        return "_".join(parts[:3])
    return isolate


def infer_iss_year(isolate: str) -> tuple[str, str]:
    """Infer missing ISS year from the established sample-name/year mapping."""
    if "IIIF" in isolate:
        return "2016", "inferred_from_IIIF_name"
    if isolate.startswith("Pantoea_piersonii_F3") or isolate.startswith("Pantoea_piersonii_F8_5S") or isolate.startswith("Pantoea_piersonii_F8_8S"):
        return "2018", "inferred_from_F3_or_F8_2018_name"
    if isolate.startswith("Pantoea_piersonii_F9") or isolate.startswith("Pantoea_piersonii_F10") or isolate.startswith("Pantoea_piersonii_F8_6S"):
        return "2021", "inferred_from_F9_F10_or_F8_6S_2021_name"
    if isolate.startswith("Pantoea_piersonii_F11"):
        return "2022", "inferred_from_F11_name"
    return "", "not_applicable"


def load_years() -> dict[str, tuple[str, str, str]]:
    years: dict[str, tuple[str, str, str]] = {}
    if not ISS_YEAR_METADATA.exists():
        return years
    for row in read_tsv(ISS_YEAR_METADATA):
        year = row.get("year", "")
        if year.endswith(".0"):
            year = year[:-2]
        years[row["isolate"]] = (year, row.get("collection_date", ""), "metadata_table")
    return years


def fasta_names(path: Path) -> list[str]:
    names = []
    if not path.exists():
        return names
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                names.append(line[1:].strip().split()[0].split("|")[0])
    return names


def newick_tip_tokens(path: Path) -> list[str]:
    if not path.exists():
        return []
    text = path.read_text(errors="ignore")
    tokens = re.findall(r"(?<=[(,])\s*([^():,;]+?)\s*(?=[:),;])", text)
    cleaned = []
    for token in tokens:
        token = token.strip("'\"")
        token = token.split("{", 1)[0]
        token = token.split("[", 1)[0]
        if token:
            cleaned.append(token)
    return cleaned


def build_matrix() -> list[dict[str, object]]:
    target_rows = read_tsv(TARGET_PROFILE)
    strict_by_isolate = {r["isolate"]: r for r in read_tsv(STRICT_PROFILE)}
    locus_by_isolate = {r["isolate"]: r for r in read_tsv(LOCUS_SUMMARY)}
    years = load_years()

    rows: list[dict[str, object]] = []
    for row in target_rows:
        isolate = row["isolate"]
        signature = row["profile_signature"]
        genes = [g.strip().lower() for g in row["profile_genes"].split(",")]
        if genes != TARGET_GENES:
            raise ValueError(f"Unexpected target gene order for {isolate}: {genes}")
        if len(signature) != len(TARGET_GENES):
            raise ValueError(f"Unexpected signature length for {isolate}: {signature}")

        year, collection_date, year_source = "", "", "not_applicable"
        if row["group"] == "ISS":
            if isolate in years:
                year, collection_date, year_source = years[isolate]
            else:
                year, year_source = infer_iss_year(isolate)

        out: dict[str, object] = {
            "genome": isolate,
            "source": row["group"],
            "year": year,
            "collection_date": collection_date,
            "year_source": year_source,
            "species": species_from_isolate(isolate),
            "target_profile_cluster": row["profile_cluster"],
            "target_profile_signature": signature,
            "strict_lps_profile_cluster": strict_by_isolate.get(isolate, {}).get("profile_cluster", ""),
            "strict_lps_profile_signature": strict_by_isolate.get(isolate, {}).get("profile_signature", ""),
        }
        for gene, value in zip(TARGET_GENES, signature):
            out[OUT_GENE_NAMES[gene]] = int(value)

        locus = locus_by_isolate.get(isolate, {})
        out.update(
            {
                "complete_rfbABCD": locus.get("complete_rfbABCD", ""),
                "rfbABCD_colocalized_15kb": locus.get("rfbABCD_colocalized_15kb", ""),
                "rfb_order_by_coordinate": locus.get("rfb_order_by_coordinate", ""),
                "rfb_span_bp": locus.get("rfb_span_bp", ""),
                "wzm_wzt_colocalized_10kb": locus.get("wzm_wzt_colocalized_10kb", ""),
                "has_locus_coordinate_summary": int(isolate in locus_by_isolate),
            }
        )
        rows.append(out)
    return rows


def build_tree_inventory(matrix_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    genomes = {str(r["genome"]) for r in matrix_rows}
    seq_map = {}
    if SEQUENCE_ID_MAP.exists():
        for row in read_tsv(SEQUENCE_ID_MAP):
            seq_map[row["safe_id"]] = row["isolate"]

    candidates: list[tuple[str, Path, str]] = [
        ("diaA_nucleotide_sequences", DIAA_FFN, "fasta_sequence_set"),
        ("diaA_protein_sequences", DIAA_FAA, "fasta_sequence_set"),
    ]
    for gene in ["rfbA", "rfbB", "rfbC", "rfbD", "waaL"]:
        candidates.append((f"{gene}_gene_tree", ROOT / f"{gene}.codon.aln.fasta.treefile", "gene_tree"))
        candidates.append((f"{gene}_relax_labeled_tree", ROOT / f"{gene}.relax_labeled.tree", "gene_tree_labeled"))

    inventory = []
    for name, path, kind in candidates:
        if kind == "fasta_sequence_set":
            raw_tips = fasta_names(path)
            isolates = raw_tips
        else:
            raw_tips = newick_tip_tokens(path)
            isolates = [seq_map.get(t, t) for t in raw_tips]
        isolate_set = set(isolates)
        matched = sorted(isolate_set & genomes)
        missing = sorted(genomes - isolate_set)
        extra = sorted(isolate_set - genomes)
        group_counts = Counter(str(r["source"]) for r in matrix_rows if str(r["genome"]) in isolate_set)
        inventory.append(
            {
                "candidate": name,
                "kind": kind,
                "path": str(path),
                "exists": int(path.exists()),
                "raw_tip_or_sequence_count": len(raw_tips),
                "unique_isolate_count": len(isolate_set),
                "matched_to_92_genome_matrix": len(matched),
                "matched_ISS": group_counts.get("ISS", 0),
                "matched_Earth": group_counts.get("Earth", 0),
                "missing_from_candidate_count": len(missing),
                "extra_not_in_matrix_count": len(extra),
                "recommended_for_phylogeny_aware_test": "yes_if_tree_is_built" if name.startswith("diaA_") and len(matched) >= 70 else "no",
                "notes": (
                    "Sequence set, not a tree; covers only genomes with diaA sequence and would need alignment/tree construction."
                    if name.startswith("diaA_")
                    else "Gene-specific tree; unsuitable for whole-dataset ISS/Earth permutation because it includes only genomes with this gene."
                ),
            }
        )
    return inventory


def first_existing(paths: list[Path]) -> str:
    for path in paths:
        if path.exists():
            return str(path)
    return ""


def build_genome_file_index(matrix_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in matrix_rows:
        genome = str(row["genome"])
        result_dir = PROJECT_RESULTS / genome
        bakta = result_dir / "bakta_output"
        eggnog = result_dir / "eggnog_out"
        rows.append(
            {
                "genome": genome,
                "source": row["source"],
                "year": row["year"],
                "species": row["species"],
                "result_dir": str(result_dir) if result_dir.exists() else "",
                "fna": first_existing([bakta / f"{genome}.fna"]),
                "faa": first_existing([bakta / f"{genome}.faa"]),
                "ffn": first_existing([bakta / f"{genome}.ffn"]),
                "gff3": first_existing([bakta / f"{genome}.gff3"]),
                "gbff": first_existing([bakta / f"{genome}.gbff"]),
                "bakta_tsv": first_existing([bakta / f"{genome}.tsv"]),
                "eggnog_annotations": first_existing(
                    [
                        eggnog / f"{genome}.emapper.annotations",
                        eggnog / "emapper.annotations",
                    ]
                ),
                "has_fna": int((bakta / f"{genome}.fna").exists()),
                "has_bakta_tsv": int((bakta / f"{genome}.tsv").exists()),
                "has_eggnog_annotations": int(
                    (eggnog / f"{genome}.emapper.annotations").exists()
                    or (eggnog / "emapper.annotations").exists()
                ),
            }
        )
    return rows


def write_summary(
    matrix_rows: list[dict[str, object]],
    inventory: list[dict[str, object]],
    genome_index: list[dict[str, object]],
) -> None:
    source_counts = Counter(str(r["source"]) for r in matrix_rows)
    species_counts = Counter(str(r["species"]) for r in matrix_rows)
    year_counts = Counter(str(r["year"]) for r in matrix_rows if r["source"] == "ISS")
    target_clusters = {(r["source"], r["target_profile_cluster"]) for r in matrix_rows}
    strict_clusters = {(r["source"], r["strict_lps_profile_cluster"]) for r in matrix_rows}

    lines = [
        "# Phylogeny-aware input inventory",
        "",
        "## Matrix",
        f"- Genomes in focal-gene matrix: {len(matrix_rows)}",
        f"- Source counts: {dict(source_counts)}",
        f"- ISS year counts: {dict(sorted(year_counts.items()))}",
        f"- Species represented: {len(species_counts)}",
        f"- Target genes: {', '.join(OUT_GENE_NAMES[g] for g in TARGET_GENES)}",
        f"- Exact target-profile clusters by source: ISS={len({c for s, c in target_clusters if s == 'ISS'})}, Earth={len({c for s, c in target_clusters if s == 'Earth'})}",
        f"- Exact strict-LPS-profile clusters by source: ISS={len({c for s, c in strict_clusters if s == 'ISS'})}, Earth={len({c for s, c in strict_clusters if s == 'Earth'})}",
        f"- Genome FASTA files indexed: {sum(int(r['has_fna']) for r in genome_index)}/{len(genome_index)}",
        f"- Bakta TSV files indexed: {sum(int(r['has_bakta_tsv']) for r in genome_index)}/{len(genome_index)}",
        f"- eggNOG annotation files indexed: {sum(int(r['has_eggnog_annotations']) for r in genome_index)}/{len(genome_index)}",
        "",
        "## Tree/Input Candidates",
    ]
    for row in inventory:
        lines.append(
            f"- {row['candidate']}: matched {row['matched_to_92_genome_matrix']}/92 genomes "
            f"({row['matched_ISS']} ISS, {row['matched_Earth']} Earth). {row['notes']}"
        )
    lines.extend(
        [
            "",
            "## Decision For Next Step",
            "- No whole-genome/core-genome tree covering all 92 genomes was found in `selection_phase2`.",
            "- The existing rfb/waaL trees are gene-specific and are not appropriate for a whole-dataset source-label permutation test.",
            "- The diaA sequence set covers a subset of the matrix and could be used only for a pruned/single-gene exploratory tree.",
            "- For the strongest phylogeny-aware association test, build or locate a core-genome/ANI tree covering the 92-genome matrix.",
        ]
    )
    (OUT / "input_inventory_summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    matrix_rows = build_matrix()
    matrix_fields = [
        "genome",
        "source",
        "year",
        "collection_date",
        "year_source",
        "species",
        "target_profile_cluster",
        "target_profile_signature",
        "strict_lps_profile_cluster",
        "strict_lps_profile_signature",
        "rfbA",
        "rfbB",
        "rfbC",
        "rfbD",
        "waaL",
        "wzm",
        "wzt",
        "rfaZ",
        "complete_rfbABCD",
        "rfbABCD_colocalized_15kb",
        "rfb_order_by_coordinate",
        "rfb_span_bp",
        "wzm_wzt_colocalized_10kb",
        "has_locus_coordinate_summary",
    ]
    write_tsv(OUT / "focal_gene_matrix.tsv", matrix_rows, matrix_fields)

    inventory = build_tree_inventory(matrix_rows)
    genome_index = build_genome_file_index(matrix_rows)
    inventory_fields = [
        "candidate",
        "kind",
        "path",
        "exists",
        "raw_tip_or_sequence_count",
        "unique_isolate_count",
        "matched_to_92_genome_matrix",
        "matched_ISS",
        "matched_Earth",
        "missing_from_candidate_count",
        "extra_not_in_matrix_count",
        "recommended_for_phylogeny_aware_test",
        "notes",
    ]
    write_tsv(OUT / "tree_candidate_inventory.tsv", inventory, inventory_fields)
    genome_index_fields = [
        "genome",
        "source",
        "year",
        "species",
        "result_dir",
        "fna",
        "faa",
        "ffn",
        "gff3",
        "gbff",
        "bakta_tsv",
        "eggnog_annotations",
        "has_fna",
        "has_bakta_tsv",
        "has_eggnog_annotations",
    ]
    write_tsv(OUT / "genome_file_index.tsv", genome_index, genome_index_fields)
    write_summary(matrix_rows, inventory, genome_index)
    print(OUT / "focal_gene_matrix.tsv")
    print(OUT / "tree_candidate_inventory.tsv")
    print(OUT / "genome_file_index.tsv")
    print(OUT / "input_inventory_summary.md")


if __name__ == "__main__":
    main()
