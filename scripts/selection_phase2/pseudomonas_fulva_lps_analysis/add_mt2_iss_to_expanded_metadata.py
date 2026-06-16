#!/usr/bin/env python3
"""Add MT-2 ISS Pseudomonas fulva assemblies to the expanded metadata table."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parent
BASE_META = ROOT / "genome_metadata_pfulva_expanded.tsv"
MT2 = ROOT / "mt2_pseudomonas_fulva_biosample_metadata.tsv"
OUT = ROOT / "genome_metadata_pfulva_expanded_with_mt2_iss.tsv"
MISSING = ROOT / "mt2_pfulva_without_assembly.tsv"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    base = read_tsv(BASE_META)
    mt2 = read_tsv(MT2)
    fields = ["accession", "strain", "isolate_id", "group", "source_note"]

    seen_accessions = {row["accession"] for row in base}
    seen_isolates = {row["isolate_id"] for row in base}

    added = []
    missing = []
    for row in mt2:
        sample = row["sample_name"]
        accession_field = row["assembly_accessions"].strip()
        if not accession_field:
            missing.append(row)
            continue
        accessions = [x.strip() for x in accession_field.replace(",", ";").split(";") if x.strip()]
        accession = accessions[0]
        isolate_id = f"Pseudomonas_fulva_{sample}"
        if accession in seen_accessions or isolate_id in seen_isolates:
            continue
        source_bits = [
            "International Space Station MT-2 isolate",
            f"sample={sample}",
            f"BioSample={row['biosample']}",
            f"collection_date={row['collection_date']}",
            f"flight={row['flight']}",
        ]
        if row.get("iss_location_from_runinfo"):
            source_bits.append(f"runinfo_location_code={row['iss_location_from_runinfo']}")
        added.append(
            {
                "accession": accession,
                "strain": f"Pseudomonas fulva strain {sample}",
                "isolate_id": isolate_id,
                "group": "ISS",
                "source_note": "; ".join(source_bits),
            }
        )
        seen_accessions.add(accession)
        seen_isolates.add(isolate_id)

    write_tsv(OUT, base + added, fields)
    if missing:
        write_tsv(MISSING, missing, list(missing[0].keys()))

    print(OUT)
    print(f"Base rows: {len(base)}")
    print(f"MT-2 ISS rows added: {len(added)}")
    print(f"MT-2 rows without assembly accession: {len(missing)}")
    if missing:
        print(MISSING)


if __name__ == "__main__":
    main()
