#!/usr/bin/env python3
"""Download selected Pseudomonas fulva Earth/reference assemblies."""

from __future__ import annotations

import csv
import gzip
import shutil
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SELECTED = ROOT / "pfulva_earth_references_selected.tsv"
GENOMES = ROOT / "genomes"


def safe_folder(accession: str, isolate_id: str) -> str:
    return f"{accession}_{isolate_id}"


def main() -> None:
    if not SELECTED.exists():
        raise SystemExit(f"Missing {SELECTED}. Run select_pfulva_earth_references.py first.")

    with SELECTED.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    for row in rows:
        accession = row["assembly_accession"]
        isolate_id = row["isolate_id"]
        url = row["genomic_fna_url"]
        if not url:
            print(f"Skipping {accession}: no genomic_fna_url")
            continue
        outdir = GENOMES / safe_folder(accession, isolate_id)
        outdir.mkdir(parents=True, exist_ok=True)
        fasta = outdir / f"{isolate_id}.fasta"
        gz_path = outdir / f"{isolate_id}.fasta.gz"
        tmp = outdir / f"{isolate_id}.fasta.gz.tmp"
        if fasta.exists() and fasta.stat().st_size > 0:
            print(f"Already exists: {fasta}")
            continue
        print(f"Downloading {accession} -> {gz_path}")
        urllib.request.urlretrieve(url, tmp)
        tmp.replace(gz_path)
        with gzip.open(gz_path, "rb") as src, fasta.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        if fasta.stat().st_size == 0:
            raise RuntimeError(f"Downloaded empty FASTA for {accession}")

    print("Downloaded selected Pseudomonas fulva Earth/reference assemblies.")


if __name__ == "__main__":
    main()
