#!/usr/bin/env python3
import csv
import gzip
import json
import os
import shutil
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DOWNLOADS = ROOT / "genbank"
META = ROOT / "expanded_genbank_metadata.tsv"

# NCBI query based on the collaborator's suggested dataset:
# taxon 47880 plus Jet Propulsion Laboratory / JPL text. The query is kept
# editable because NCBI metadata wording can vary between records.
ASSEMBLY_QUERY = (
    'txid47880[Organism:exp] AND '
    '("Jet Propulsion Laboratory"[All Fields] OR "Jet Propulsion lab"[All Fields] OR JPL[All Fields])'
)


def fetch_json(url):
    with urllib.request.urlopen(url, timeout=60) as handle:
        return json.loads(handle.read().decode("utf-8"))


def fetch_text(url):
    with urllib.request.urlopen(url, timeout=60) as handle:
        return handle.read().decode("utf-8")


def ncbi_get(path, params):
    params = dict(params)
    params.setdefault("retmode", "json")
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" + path + "?" + urllib.parse.urlencode(params)
    return fetch_json(url)


def esearch_assemblies(query):
    data = ncbi_get(
        "esearch.fcgi",
        {
            "db": "assembly",
            "term": query,
            "retmax": 500,
        },
    )
    return data.get("esearchresult", {}).get("idlist", [])


def esummary_assemblies(ids):
    if not ids:
        return []
    data = ncbi_get(
        "esummary.fcgi",
        {
            "db": "assembly",
            "id": ",".join(ids),
            "retmax": len(ids),
        },
    )
    result = data.get("result", {})
    return [result[i] for i in result.get("uids", [])]


def safe_name(text):
    keep = []
    for char in text:
        if char.isalnum() or char in "._-":
            keep.append(char)
        else:
            keep.append("_")
    return "_".join("".join(keep).split("_"))


def ftp_to_https(path):
    if path.startswith("ftp://"):
        return "https://" + path[len("ftp://") :]
    return path


def download_file(url, out):
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    if out.exists() and out.stat().st_size > 0:
        return
    with urllib.request.urlopen(url, timeout=120) as response, open(tmp, "wb") as handle:
        shutil.copyfileobj(response, handle)
    tmp.replace(out)


def decompress_gzip(src, dest):
    if dest.exists() and dest.stat().st_size > 0:
        return
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    with gzip.open(src, "rb") as inp, open(tmp, "wb") as out:
        shutil.copyfileobj(inp, out)
    tmp.replace(dest)


def assembly_files(summary):
    ftp = summary.get("ftppath_genbank") or summary.get("ftppath_refseq") or ""
    if not ftp:
        return None, None
    base = ftp.rstrip("/").split("/")[-1]
    https = ftp_to_https(ftp.rstrip("/"))
    gbff = f"{https}/{base}_genomic.gbff.gz"
    fna = f"{https}/{base}_genomic.fna.gz"
    return gbff, fna


def main():
    query = os.environ.get("NCBI_ASSEMBLY_QUERY", ASSEMBLY_QUERY)
    print(f"NCBI Assembly query: {query}")
    ids = esearch_assemblies(query)
    print(f"Found {len(ids)} assembly records")
    summaries = esummary_assemblies(ids)

    rows = []
    DOWNLOADS.mkdir(exist_ok=True)
    for item in summaries:
        accession = item.get("assemblyaccession", "")
        organism = item.get("organism", "")
        biosample = item.get("biosampleaccn", "")
        submitter = item.get("submitterorganization", "")
        species = item.get("speciesname", "")
        strain = item.get("infraspecifickey", "") + " " + item.get("infraspecificname", "")
        strain = strain.strip() or item.get("assemblyname", "")
        isolate_id = safe_name(f"{accession}_{organism}_{strain}")[:160]
        gbff_url, fna_url = assembly_files(item)
        if not gbff_url:
            print(f"Skipping {accession}: no GenBank/RefSeq FTP path", file=sys.stderr)
            continue

        outdir = DOWNLOADS / isolate_id
        gbff_gz = outdir / f"{isolate_id}.gbff.gz"
        gbff = outdir / f"{isolate_id}.gbff"
        fna_gz = outdir / f"{isolate_id}.fna.gz"

        print(f"Downloading {accession}: {organism}")
        download_file(gbff_url, gbff_gz)
        decompress_gzip(gbff_gz, gbff)
        try:
            download_file(fna_url, fna_gz)
        except Exception as exc:
            print(f"Warning: FASTA download failed for {accession}: {exc}", file=sys.stderr)

        rows.append(
            {
                "assembly_accession": accession,
                "biosample": biosample,
                "organism": organism,
                "species": species,
                "strain": strain,
                "isolate_id": isolate_id,
                "group": "ISS_or_JPL_Pseudomonas",
                "source_note": submitter,
                "gbff": str(gbff),
                "fna_gz": str(fna_gz) if fna_gz.exists() else "",
                "ftp_path": item.get("ftppath_genbank") or item.get("ftppath_refseq") or "",
            }
        )
        time.sleep(0.34)

    with open(META, "w", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "assembly_accession",
            "biosample",
            "organism",
            "species",
            "strain",
            "isolate_id",
            "group",
            "source_note",
            "gbff",
            "fna_gz",
            "ftp_path",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {META}")


if __name__ == "__main__":
    main()
