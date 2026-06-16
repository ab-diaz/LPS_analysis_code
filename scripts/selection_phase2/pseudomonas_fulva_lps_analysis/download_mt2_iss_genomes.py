#!/usr/bin/env python3
"""Download MT-2 ISS Pseudomonas fulva genome FASTAs from NCBI FTP."""

from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
import time
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parent
META = ROOT / "genome_metadata_pfulva_expanded_with_mt2_iss.tsv"
GENOMES = ROOT / "genomes"
NCBI_EMAIL = os.environ.get("NCBI_EMAIL")


if not NCBI_EMAIL:
    raise SystemExit("Set NCBI_EMAIL before querying NCBI, e.g. export NCBI_EMAIL=name@example.org")


def eutils_url(endpoint: str, params: dict[str, str]) -> str:
    base = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}.fcgi"
    params = dict(params)
    params.setdefault("tool", "kalamiella_pfulva_mt2_download")
    params.setdefault("email", NCBI_EMAIL)
    return base + "?" + urllib.parse.urlencode(params)


def fetch_assembly_summary(accessions: list[str]) -> dict[str, dict]:
    out = {}
    for i in range(0, len(accessions), 100):
        chunk = accessions[i : i + 100]
        search = eutils_url(
            "esearch",
            {
                "db": "assembly",
                "term": " OR ".join(f"{acc}[Assembly Accession]" for acc in chunk),
                "retmode": "json",
                "retmax": "200",
            },
        )
        with urllib.request.urlopen(search, timeout=90) as handle:
            ids = json.load(handle)["esearchresult"].get("idlist", [])
        if not ids:
            continue
        summary = eutils_url(
            "esummary",
            {
                "db": "assembly",
                "id": ",".join(ids),
                "retmode": "json",
            },
        )
        with urllib.request.urlopen(summary, timeout=90) as handle:
            data = json.load(handle)["result"]
        for uid in data.get("uids", []):
            row = data[uid]
            syn = row.get("synonym", {}) or {}
            for acc in [syn.get("refseq", ""), syn.get("genbank", ""), row.get("assemblyaccession", "")]:
                if acc:
                    out[acc] = row
        time.sleep(0.34)
    return out


def fna_url(summary: dict) -> str:
    ftp = summary.get("ftppath_refseq") or summary.get("ftppath_genbank") or ""
    if not ftp:
        return ""
    https = ftp.replace("ftp://", "https://").rstrip("/")
    base = https.split("/")[-1]
    return f"{https}/{base}_genomic.fna.gz"


def main() -> None:
    with META.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    targets = [row for row in rows if row["group"] == "ISS"]
    accessions = [row["accession"] for row in targets if row["accession"].startswith(("GCF_", "GCA_"))]
    summaries = fetch_assembly_summary(accessions)
    missing_summary = []

    for row in targets:
        accession = row["accession"]
        isolate_id = row["isolate_id"]
        outdir = GENOMES / f"{accession}_{isolate_id}"
        fasta = outdir / f"{isolate_id}.fasta"
        if fasta.exists() and fasta.stat().st_size > 0:
            print(f"Already exists: {fasta}")
            continue
        summary = summaries.get(accession)
        if not summary:
            missing_summary.append(accession)
            continue
        url = fna_url(summary)
        if not url:
            missing_summary.append(accession)
            continue
        outdir.mkdir(parents=True, exist_ok=True)
        gz_path = outdir / f"{isolate_id}.fasta.gz"
        tmp = outdir / f"{isolate_id}.fasta.gz.tmp"
        print(f"Downloading {accession} -> {fasta}")
        urllib.request.urlretrieve(url, tmp)
        tmp.replace(gz_path)
        with gzip.open(gz_path, "rb") as src, fasta.open("wb") as dst:
            shutil.copyfileobj(src, dst)

    if missing_summary:
        print("WARNING: Missing assembly summaries or FASTA URLs for:")
        for accession in missing_summary:
            print(accession)
    print("MT-2 ISS download step complete.")


if __name__ == "__main__":
    main()
