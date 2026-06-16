#!/usr/bin/env python3
"""Select Earth/reference Pseudomonas fulva assemblies for a side analysis.

The script queries NCBI Assembly and BioSample metadata, excludes ISS/spacecraft
records from the Earth-reference panel, ranks assemblies by quality/metadata,
and writes both the full candidate table and a conservative selected set.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT_ALL = ROOT / "pfulva_earth_reference_candidates.tsv"
OUT_SELECTED = ROOT / "pfulva_earth_references_selected.tsv"
OUT_METADATA = ROOT / "genome_metadata_pfulva_expanded.tsv"
MAX_SELECTED = 15

NCBI_TOOL = "kalamiella_pfulva_reference_selection"
NCBI_EMAIL = "your.email@example.com"
ASSEMBLY_TERM = 'txid47880[Organism:exp] AND "Pseudomonas fulva"[Organism]'

SPACE_TERMS = re.compile(
    r"iss|international space station|spacecraft|mars odyssey|odyssey|cleanroom|clean room|jpl|jet propulsion",
    re.I,
)
BAD_TERMS = re.compile(r"metagenome|metagenomic|uncultured|MAG|TPA_asm", re.I)

MANUAL_PRIORITY = {
    "NBRC": "type strain / taxonomic anchor",
    "DSM 17717": "type strain / taxonomic anchor",
    "12-X": "published complete/reference-quality environmental strain",
    "MTT5": "published complete maize phyllosphere strain",
    "FDAARGOS_167": "complete clinically associated Earth isolate; useful but check ANI carefully",
}

# A balanced Earth/reference panel is preferable to the top 15 by assembly
# completeness, because many complete P. fulva genomes are human-associated.
# These accessions prioritize source diversity while retaining several complete
# reference-quality genomes.
CURATED_EARTH_ACCESSIONS = [
    "GCF_050517415.1",  # MTT5, maize phyllosphere, complete
    "GCF_000213805.1",  # 12-X, complete/reference-quality historical comparator
    "GCF_002688705.1",  # SB1, fungal/plant-associated environment, complete
    "GCF_023517795.1",  # ZJU1, silkworm/mulberry-associated, complete
    "GCF_040212335.1",  # OsEnb_ALM_C17, environmental/plant-associated metadata
    "GCF_050155465.1",  # T23, Amaranthus-associated, complete
    "GCF_030077585.1",  # OS-1, India, chromosome-level
    "GCF_003205395.1",  # LB-090714, Chicago environmental/freshwater set
    "GCF_003936165.1",  # PF_122, Pakistan
    "GCF_000834555.1",  # MEJ086, Arabidopsis-associated
    "GCF_011777485.1",  # PS9.1, maize, South Africa
    "GCF_047468835.1",  # Pf_FISH/3a, Antarctica
    "GCF_051216675.1",  # SMV3, citrus-associated
    "GCF_002951475.1",  # FDAARGOS_167, complete clinical Earth comparator
    "GCF_015679205.1",  # ZDHY316, complete clinical Earth comparator
]


def ncbi_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=90) as handle:
        return json.load(handle)


def ncbi_text(url: str) -> str:
    with urllib.request.urlopen(url, timeout=90) as handle:
        return handle.read().decode("utf-8", errors="replace")


def eutils_url(endpoint: str, params: dict[str, str]) -> str:
    base = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}.fcgi"
    params = dict(params)
    params.setdefault("tool", NCBI_TOOL)
    params.setdefault("email", NCBI_EMAIL)
    return base + "?" + urllib.parse.urlencode(params)


def esearch_assembly_ids() -> list[str]:
    url = eutils_url(
        "esearch",
        {
            "db": "assembly",
            "term": ASSEMBLY_TERM,
            "retmode": "json",
            "retmax": "500",
        },
    )
    data = ncbi_json(url)
    return data["esearchresult"].get("idlist", [])


def esummary_assemblies(ids: list[str]) -> list[dict]:
    rows = []
    for i in range(0, len(ids), 100):
        chunk = ids[i : i + 100]
        url = eutils_url(
            "esummary",
            {
                "db": "assembly",
                "id": ",".join(chunk),
                "retmode": "json",
            },
        )
        data = ncbi_json(url)
        result = data.get("result", {})
        for uid in result.get("uids", []):
            rows.append(result[uid])
        time.sleep(0.34)
    return rows


def biosample_metadata(accessions: list[str]) -> dict[str, dict[str, str]]:
    meta = {}
    wanted = [a for a in accessions if a]
    for i in range(0, len(wanted), 50):
        chunk = wanted[i : i + 50]
        url = eutils_url(
            "efetch",
            {
                "db": "biosample",
                "id": ",".join(chunk),
                "rettype": "full",
                "retmode": "xml",
            },
        )
        try:
            xml = ncbi_text(url)
        except Exception:
            time.sleep(0.34)
            continue
        try:
            root = ET.fromstring(xml)
        except ET.ParseError:
            time.sleep(0.34)
            continue
        for sample in root.findall(".//BioSample"):
            accession = sample.attrib.get("accession", "")
            attrs = {}
            for attr in sample.findall(".//Attribute"):
                name = attr.attrib.get("attribute_name") or attr.attrib.get("harmonized_name") or ""
                if name:
                    attrs[name.lower()] = (attr.text or "").strip()
            desc = " ".join((el.text or "").strip() for el in sample.findall(".//Description/*") if el.text)
            meta[accession] = {
                "strain": attrs.get("strain", ""),
                "isolate": attrs.get("isolate", ""),
                "sample_name": attrs.get("sample name", ""),
                "isolation_source": attrs.get("isolation source", ""),
                "geo_loc_name": attrs.get("geo_loc_name", ""),
                "collection_date": attrs.get("collection date", ""),
                "host": attrs.get("host", ""),
                "description": desc,
            }
        time.sleep(0.34)
    return meta


def clean_id(text: str) -> str:
    text = re.sub(r"^Pseudomonas fulva\s*", "", text, flags=re.I)
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_")
    text = re.sub(r"_+", "_", text)
    return text or "unknown"


def strain_from_organism(organism: str) -> str:
    text = organism
    text = re.sub(r"\s*\([^)]*\)\s*$", "", text)
    text = re.sub(r"^Pseudomonas fulva\s*", "", text, flags=re.I).strip()
    if not text or text.lower() in {"strain", "isolate"}:
        return ""
    return text


def is_assembly_placeholder(value: str) -> bool:
    return bool(re.fullmatch(r"ASM\d+v\d+|ASM\d+|GCA_\d+\.\d+|GCF_\d+\.\d+", value or ""))


def normalize_strain_key(value: str) -> str:
    value = value.lower()
    value = re.sub(r"^pseudomonas fulva\s*", "", value)
    value = re.sub(r"[^a-z0-9]+", "", value)
    return value


def add_type_strain_aliases(keys: set[str], value: str) -> None:
    norm = normalize_strain_key(value)
    if not norm:
        return
    if "dsm17717" in norm or "nbrc16637" in norm:
        keys.update(
            {
                "dsm17717",
                "nbrc16637",
                "nbrc16637dsm17717",
                "pseudomonasfulvanbrc16637dsm17717",
            }
        )


def ftp_to_genomic_fna_url(ftp_path: str) -> str:
    if not ftp_path:
        return ""
    https = ftp_path.replace("ftp://", "https://")
    base = https.rstrip("/").split("/")[-1]
    return f"{https}/{base}_genomic.fna.gz"


def quality_score(row: dict) -> tuple:
    level_rank = {
        "Complete Genome": 0,
        "Chromosome": 1,
        "Scaffold": 2,
        "Contig": 3,
    }.get(row["assembly_level"], 4)
    source_rank = 0 if row["refseq_accession"] else 1
    manual_rank = 0 if row["manual_priority"] else 1
    metadata_rank = 0 if (row["isolation_source"] or row["geo_loc_name"]) else 1
    length = int(row["total_length"] or 0)
    contigs = int(row["contigs"] or 999999)
    return (manual_rank, level_rank, source_rank, metadata_rank, contigs, -length, row["assembly_accession"])


def classify_group(row: dict) -> str:
    text = " ".join(
        [
            row["organism"],
            row["strain"],
            row["assembly_name"],
            row["isolation_source"],
            row["geo_loc_name"],
            row["description"],
            row["submitter"],
        ]
    )
    if SPACE_TERMS.search(text):
        if re.search(r"iss|international space station", text, re.I):
            return "ISS"
        return "Earth_spacecraft_associated"
    return "Earth_reference"


def manual_priority(row: dict) -> str:
    text = " ".join([row["organism"], row["strain"], row["assembly_name"]])
    for key, reason in MANUAL_PRIORITY.items():
        if key.lower() in text.lower():
            return reason
    return ""


def row_from_summary(summary: dict, bs_meta: dict[str, dict[str, str]]) -> dict:
    biosample = summary.get("biosampleaccn", "")
    bio = bs_meta.get(biosample, {})
    synonym = summary.get("synonym", {}) or {}
    stats = summary.get("assemblystats", {}) or {}
    refseq = synonym.get("refseq", "") or ""
    genbank = synonym.get("genbank", "") or summary.get("assemblyaccession", "")
    accession = refseq or genbank
    ftp_refseq = summary.get("ftppath_refseq", "") or ""
    ftp_genbank = summary.get("ftppath_genbank", "") or ""
    organism = summary.get("organism", "")
    assembly_name = summary.get("assemblyname", "")
    infraspecies = summary.get("infraspecies", "")
    infraspecies = re.sub(r"^strain=", "", infraspecies)
    strain = (
        infraspecies
        or bio.get("strain", "")
        or bio.get("isolate", "")
        or strain_from_organism(organism)
        or bio.get("sample_name", "")
        or assembly_name
        or accession
    )
    if is_assembly_placeholder(strain) and strain_from_organism(organism):
        strain = strain_from_organism(organism)
    row = {
        "assembly_accession": accession,
        "refseq_accession": refseq,
        "genbank_accession": genbank,
        "organism": organism,
        "strain": strain,
        "assembly_name": assembly_name,
        "assembly_level": summary.get("assemblystatus", ""),
        "refseq_category": summary.get("refseq_category", ""),
        "biosample": biosample,
        "bioproject": summary.get("bioprojectaccn", ""),
        "submitter": summary.get("submitterorganization", ""),
        "total_length": str(stats.get("totalsequencelength", "")),
        "contigs": str(stats.get("numberofcontigs", "")),
        "gc_percent": str(stats.get("gcpercent", "")),
        "ftp_path": ftp_refseq or ftp_genbank,
        "genomic_fna_url": ftp_to_genomic_fna_url(ftp_refseq or ftp_genbank),
        "isolation_source": bio.get("isolation_source", ""),
        "geo_loc_name": bio.get("geo_loc_name", ""),
        "collection_date": bio.get("collection_date", ""),
        "host": bio.get("host", ""),
        "biosample_strain": bio.get("strain", ""),
        "biosample_isolate": bio.get("isolate", ""),
        "biosample_sample_name": bio.get("sample_name", ""),
        "description": bio.get("description", ""),
        "exclude_reason": "",
        "manual_priority": "",
        "group": "",
        "isolate_id": "",
        "source_note": "",
    }
    row["manual_priority"] = manual_priority(row)
    row["group"] = classify_group(row)
    id_source = row["strain"]
    if is_assembly_placeholder(id_source):
        id_source = row["assembly_accession"]
    row["isolate_id"] = "Pseudomonas_fulva_" + clean_id(id_source)
    note_bits = [x for x in [row["geo_loc_name"], row["isolation_source"], row["host"], row["manual_priority"]] if x]
    row["source_note"] = "; ".join(note_bits)
    text = " ".join(row.values())
    if BAD_TERMS.search(text):
        row["exclude_reason"] = "metagenome_or_uncultured_record"
    elif "Pseudomonas fulva" not in organism:
        row["exclude_reason"] = "organism_name_not_pseudomonas_fulva"
    elif not row["genomic_fna_url"]:
        row["exclude_reason"] = "missing_genomic_fna_url"
    return row


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    fieldnames = [
        "assembly_accession",
        "refseq_accession",
        "genbank_accession",
        "organism",
        "strain",
        "assembly_name",
        "assembly_level",
        "refseq_category",
        "biosample",
        "bioproject",
        "submitter",
        "total_length",
        "contigs",
        "gc_percent",
        "group",
        "isolate_id",
        "source_note",
        "isolation_source",
        "geo_loc_name",
        "collection_date",
        "host",
        "biosample_strain",
        "biosample_isolate",
        "biosample_sample_name",
        "manual_priority",
        "exclude_reason",
        "genomic_fna_url",
        "ftp_path",
        "description",
    ]

    try:
        ids = esearch_assembly_ids()
        summaries = esummary_assemblies(ids)
        biosamples = sorted({s.get("biosampleaccn", "") for s in summaries if s.get("biosampleaccn")})
        bs_meta = biosample_metadata(biosamples)
        rows = [row_from_summary(s, bs_meta) for s in summaries]
        rows.sort(key=quality_score)
        write_tsv(OUT_ALL, rows, fieldnames)
    except Exception as exc:
        if not OUT_ALL.exists():
            raise
        print(f"WARNING: NCBI query failed ({exc}); using cached {OUT_ALL}")
        rows = load_tsv(OUT_ALL)

    original_meta = ROOT / "genome_metadata.tsv"
    existing_accessions = set()
    existing_strain_keys = set()
    with original_meta.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            existing_accessions.add(row["accession"])
            existing_strain_keys.add(normalize_strain_key(row["strain"]))
            existing_strain_keys.add(normalize_strain_key(row["isolate_id"]))
            add_type_strain_aliases(existing_strain_keys, row["strain"])
            add_type_strain_aliases(existing_strain_keys, row["isolate_id"])

    selected = []
    seen_ids = set()
    seen_strain_keys = set(existing_strain_keys)

    rows_by_accession = {row["assembly_accession"]: row for row in rows}
    ordered_rows = []
    for accession in CURATED_EARTH_ACCESSIONS:
        row = rows_by_accession.get(accession)
        if row:
            ordered_rows.append(row)
    ordered_rows.extend(row for row in rows if row["assembly_accession"] not in CURATED_EARTH_ACCESSIONS)

    for row in ordered_rows:
        if row["exclude_reason"] or row["group"] != "Earth_reference":
            continue
        if row["assembly_accession"] in existing_accessions:
            continue
        strain_key = normalize_strain_key(row["strain"])
        if strain_key and strain_key in seen_strain_keys:
            continue
        if row["isolate_id"] in seen_ids:
            continue
        selected.append(row)
        seen_ids.add(row["isolate_id"])
        if strain_key:
            seen_strain_keys.add(strain_key)
        if len(selected) >= MAX_SELECTED:
            break
    write_tsv(OUT_SELECTED, selected, fieldnames)

    expanded_rows = []
    with original_meta.open(newline="", encoding="utf-8") as handle:
        expanded_rows.extend(csv.DictReader(handle, delimiter="\t"))
    known_accessions = {r["accession"] for r in expanded_rows}
    for row in selected:
        if row["assembly_accession"] in known_accessions:
            continue
        expanded_rows.append(
            {
                "accession": row["assembly_accession"],
                "strain": row["organism"],
                "isolate_id": row["isolate_id"],
                "group": "Earth_reference",
                "source_note": row["source_note"] or "Earth/reference Pseudomonas fulva assembly",
            }
        )
    write_tsv(
        OUT_METADATA,
        expanded_rows,
        ["accession", "strain", "isolate_id", "group", "source_note"],
    )

    print(OUT_ALL)
    print(OUT_SELECTED)
    print(OUT_METADATA)
    print(f"Selected Earth/reference assemblies: {len(selected)}")


if __name__ == "__main__":
    main()
