#!/usr/bin/env python3
"""
Extract MT-2 Pseudomonas fulva SRA/BioSample metadata from NCBI.

The script uses public NCBI endpoints:
1. SRA RunInfo for PRJNA690512 / SRP303310.
2. BioSample XML for each linked BioSample accession.
3. Assembly E-utilities for BioSample-linked WGS/assembly accessions when available.

It writes a table with sample/isolate names, SRA runs, BioSample accessions,
collection dates, ISS locations, and assembly/WGS accessions.
"""
from __future__ import annotations

import csv
import io
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


OUTDIR = Path(__file__).resolve().parent
OUT = OUTDIR / "mt2_pseudomonas_fulva_biosample_metadata.tsv"
SUMMARY = OUTDIR / "mt2_pseudomonas_fulva_biosample_metadata_summary.md"

RUNINFO_URL = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIOPROJECT = "PRJNA690512"
STUDY = "SRP303310"
ORGANISM = "Pseudomonas fulva"
SLEEP_SECONDS = 0.34


def fetch_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "pfulva-metadata-extractor/1.0"})
    with urllib.request.urlopen(req, timeout=60) as handle:
        return handle.read().decode("utf-8", errors="replace")


def get_runinfo() -> list[dict[str, str]]:
    query = urllib.parse.urlencode({"acc": STUDY})
    text = fetch_text(f"{RUNINFO_URL}?{query}")
    rows = list(csv.DictReader(io.StringIO(text)))
    return [
        row
        for row in rows
        if row.get("ScientificName", "").strip().lower() == ORGANISM.lower()
        or row.get("Organism", "").strip().lower() == ORGANISM.lower()
    ]


def get_biosample_attributes(accession: str) -> dict[str, str]:
    if not accession:
        return {}
    query = urllib.parse.urlencode({"db": "biosample", "id": accession, "retmode": "xml"})
    xml_text = fetch_text(f"{EUTILS}/efetch.fcgi?{query}")
    root = ET.fromstring(xml_text)
    attrs: dict[str, str] = {}
    sample = root.find(".//BioSample")
    if sample is not None:
        attrs["biosample_accession"] = sample.attrib.get("accession", accession)
        attrs["biosample_id"] = sample.attrib.get("id", "")
        attrs["biosample_publication_date"] = sample.attrib.get("publication_date", "")
        attrs["biosample_last_update"] = sample.attrib.get("last_update", "")
    for attr in root.findall(".//Attribute"):
        key = attr.attrib.get("attribute_name") or attr.attrib.get("harmonized_name") or ""
        key = key.strip().lower().replace(" ", "_")
        if key:
            attrs[key] = (attr.text or "").strip()
    return attrs


def assembly_accessions_for_biosample(biosample: str) -> str:
    if not biosample:
        return ""
    term = f"{biosample}[BioSample]"
    query = urllib.parse.urlencode({"db": "assembly", "term": term, "retmode": "json", "retmax": 20})
    try:
        search_text = fetch_text(f"{EUTILS}/esearch.fcgi?{query}")
        import json

        ids = json.loads(search_text).get("esearchresult", {}).get("idlist", [])
        if not ids:
            return ""
        time.sleep(SLEEP_SECONDS)
        summary_query = urllib.parse.urlencode({"db": "assembly", "id": ",".join(ids), "retmode": "json"})
        summary_text = fetch_text(f"{EUTILS}/esummary.fcgi?{summary_query}")
        data = json.loads(summary_text).get("result", {})
        accs = []
        for uid in data.get("uids", []):
            rec = data.get(uid, {})
            acc = rec.get("assemblyaccession", "")
            if acc:
                accs.append(acc)
        return ";".join(sorted(set(accs)))
    except Exception:
        return ""


def first_present(row: dict[str, str], keys: list[str]) -> str:
    for key in keys:
        if row.get(key):
            return row[key]
    return ""


def main() -> None:
    rows = []
    runinfo_rows = get_runinfo()
    for i, run in enumerate(runinfo_rows, start=1):
        biosample = run.get("BioSample", "")
        attrs = get_biosample_attributes(biosample)
        if i < len(runinfo_rows):
            time.sleep(SLEEP_SECONDS)
        assembly_acc = assembly_accessions_for_biosample(biosample)
        time.sleep(SLEEP_SECONDS)

        sample_name = first_present(
            attrs,
            ["isolate", "strain", "sample_name", "title"],
        ) or run.get("SampleName", "")
        collection_date = first_present(attrs, ["collection_date", "collection_date_sam", "collection-date"])
        isolation_source = first_present(attrs, ["isolation_source", "source_material_id"])
        geo_loc = first_present(attrs, ["geo_loc_name", "geographic_location"])

        rows.append(
            {
                "sample_name": sample_name,
                "run": run.get("Run", ""),
                "experiment": run.get("Experiment", ""),
                "sra_sample": run.get("Sample", ""),
                "biosample": biosample,
                "bioproject": run.get("BioProject", BIOPROJECT),
                "assembly_accessions": assembly_acc,
                "collection_date": collection_date,
                "flight": "F8" if sample_name.startswith("F8_") else "",
                "iss_location_from_runinfo": run.get("LibraryName", ""),
                "isolation_source": isolation_source,
                "geo_loc_name": geo_loc,
                "organism": run.get("ScientificName", ORGANISM),
                "bases": run.get("bases", ""),
                "spots": run.get("spots", ""),
                "published": run.get("ReleaseDate", ""),
            }
        )

    rows = sorted(rows, key=lambda r: (r["sample_name"], r["run"]))
    fieldnames = [
        "sample_name",
        "run",
        "experiment",
        "sra_sample",
        "biosample",
        "bioproject",
        "assembly_accessions",
        "collection_date",
        "flight",
        "iss_location_from_runinfo",
        "isolation_source",
        "geo_loc_name",
        "organism",
        "bases",
        "spots",
        "published",
    ]
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    dates = sorted({r["collection_date"] for r in rows if r["collection_date"]})
    SUMMARY.write_text(
        "\n".join(
            [
                "# MT-2 Pseudomonas fulva BioSample metadata",
                "",
                f"- BioProject queried: `{BIOPROJECT}`",
                f"- SRA study queried: `{STUDY}`",
                f"- Organism filter: `{ORGANISM}`",
                f"- Rows: {len(rows)}",
                f"- Unique BioSamples: {len({r['biosample'] for r in rows if r['biosample']})}",
                f"- Collection dates found: {', '.join(dates) if dates else 'none'}",
                f"- Output: `{OUT.name}`",
                "",
            ]
        )
    )
    print(OUT)
    print(SUMMARY)


if __name__ == "__main__":
    main()
