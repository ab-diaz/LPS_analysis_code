#!/usr/bin/env python3
import argparse
import csv
import json
import re
from pathlib import Path


ISS_RE = re.compile(r"\bISS\b|International Space Station|space station", re.IGNORECASE)
MISSING_VALUES = {"", "missing", "not applicable", "unknown", "na", "n/a", "-"}


def clean(value):
    value = "" if value is None else str(value).strip()
    return "" if value.lower() in MISSING_VALUES else value


def load_targets(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def attr_map(biosample):
    return {
        (item.get("name") or "").strip(): clean(item.get("value"))
        for item in biosample.get("attributes", [])
        if item.get("name")
    }


def biosample_text(obj, biosample, attrs):
    description = biosample.get("description") or {}
    bioprojects = (
        obj.get("assemblyInfo", {})
        .get("bioprojectLineage", [{}])[0]
        .get("bioprojects", [])
    )
    project_titles = " ".join(clean(p.get("title")) for p in bioprojects)
    parts = [
        biosample.get("geoLocName"),
        biosample.get("isolationSource"),
        biosample.get("host"),
        description.get("title"),
        description.get("comment"),
        project_titles,
        " ".join(f"{k}={v}" for k, v in attrs.items() if v),
    ]
    return " | ".join(clean(p) for p in parts if clean(p))


def load_assembly_report(path):
    by_biosample = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            obj = json.loads(line)
            assembly_info = obj.get("assemblyInfo", {})
            biosample = assembly_info.get("biosample") or {}
            accession = clean(biosample.get("accession"))
            if not accession or accession in by_biosample:
                continue
            attrs = attr_map(biosample)
            by_biosample[accession] = {
                "object": obj,
                "biosample": biosample,
                "attrs": attrs,
                "text": biosample_text(obj, biosample, attrs),
            }
    return by_biosample


def classify_group(text):
    matches = sorted(set(m.group(0) for m in ISS_RE.finditer(text)))
    if matches:
        return "ISS", "high", ";".join(matches)
    return "Earth", "medium", ""


def choose_source(biosample, attrs):
    return (
        clean(biosample.get("isolationSource"))
        or clean(attrs.get("isolation_source"))
        or clean(attrs.get("env_material"))
        or clean(attrs.get("source_material_id"))
        or clean(attrs.get("host"))
        or clean(biosample.get("host"))
    )


def choose_category(group, source, host):
    text = f"{source} {host}".lower()
    if group == "ISS":
        return "ISS"
    if any(term in text for term in ["homo sapiens", "human", "patient", "biopsy", "clinical"]):
        return "Clinical"
    if any(term in text for term in ["soil", "root", "leaf", "leaves", "plant", "seed", "fruit", "mushroom"]):
        return "Environmental"
    return "Environmental"


def load_existing_metadata(path):
    if not path or not Path(path).exists():
        return [], [], {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    return fieldnames, rows, {row.get("Sample ID", ""): row for row in rows}


def write_group_metadata(rows, path):
    fieldnames = [
        "sample",
        "group",
        "group_confidence",
        "matched_terms",
        "biosample",
        "assembly_accession",
        "organism",
        "strain",
        "location",
        "source",
        "host",
        "collection_date",
        "category",
        "notes",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_merged_plot_metadata(existing_fieldnames, existing_rows, existing_by_sample, rows, path):
    if not path:
        return
    fieldnames = existing_fieldnames or [
        "Sample ID",
        "Location",
        "Source",
        "Collection date",
        "Patient",
        "Category",
    ]
    merged = list(existing_rows)
    for row in rows:
        sample = row["sample"]
        if sample in existing_by_sample:
            continue
        out = {field: "" for field in fieldnames}
        out["Sample ID"] = sample
        if "Location" in out:
            out["Location"] = row["location"]
        if "Source" in out:
            out["Source"] = row["source"]
        if "Collection date" in out:
            out["Collection date"] = row["collection_date"]
        if "Patient" in out:
            out["Patient"] = row["host"]
        if "Category" in out:
            out["Category"] = row["category"]
        merged.append(out)

    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(merged)


def main():
    parser = argparse.ArgumentParser(
        description="Build reviewable ISS/Earth metadata from local NCBI assembly BioSample reports."
    )
    parser.add_argument("--targets", default="ncbi_annotation_targets_all_unique.tsv")
    parser.add_argument(
        "--assembly-report",
        default="resources/ncbi_dataset/data/assembly_data_report.jsonl",
    )
    parser.add_argument("--existing-metadata", default="lps_analysis/metadata.csv")
    parser.add_argument("--out", default="lps_analysis/ncbi_sample_group_metadata.tsv")
    parser.add_argument("--merged-metadata", default="lps_analysis/metadata_with_ncbi.tsv")
    args = parser.parse_args()

    targets = load_targets(args.targets)
    report = load_assembly_report(args.assembly_report)
    existing_fieldnames, existing_rows, existing_by_sample = load_existing_metadata(args.existing_metadata)

    rows = []
    for target in targets:
        sample = clean(target.get("sample"))
        biosample_acc = clean(target.get("BioSample"))
        record = report.get(biosample_acc)
        if not record:
            rows.append(
                {
                    "sample": sample,
                    "group": "Earth",
                    "group_confidence": "low",
                    "matched_terms": "",
                    "biosample": biosample_acc,
                    "assembly_accession": clean(target.get("Assembly Accession")),
                    "organism": clean(target.get("Organism Scientific Name")),
                    "strain": clean(target.get("strain")),
                    "location": "",
                    "source": "",
                    "host": "",
                    "collection_date": "",
                    "category": "Environmental",
                    "notes": "BioSample not found in local assembly report; defaulted to Earth.",
                }
            )
            continue

        obj = record["object"]
        biosample = record["biosample"]
        attrs = record["attrs"]
        group, confidence, matched_terms = classify_group(record["text"])
        source = choose_source(biosample, attrs)
        host = clean(biosample.get("host")) or clean(attrs.get("host"))
        category = choose_category(group, source, host)
        rows.append(
            {
                "sample": sample,
                "group": group,
                "group_confidence": confidence,
                "matched_terms": matched_terms,
                "biosample": biosample_acc,
                "assembly_accession": clean(target.get("Assembly Accession")) or clean(obj.get("accession")),
                "organism": clean(target.get("Organism Scientific Name")),
                "strain": clean(target.get("strain")) or clean(biosample.get("strain")),
                "location": clean(biosample.get("geoLocName")) or clean(attrs.get("geo_loc_name")),
                "source": source,
                "host": host,
                "collection_date": clean(biosample.get("collectionDate")) or clean(attrs.get("collection_date")),
                "category": category,
                "notes": "",
            }
        )

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    write_group_metadata(rows, args.out)
    write_merged_plot_metadata(
        existing_fieldnames,
        existing_rows,
        existing_by_sample,
        rows,
        args.merged_metadata,
    )

    counts = {}
    for row in rows:
        counts[row["group"]] = counts.get(row["group"], 0) + 1
    print(f"Wrote {args.out}")
    if args.merged_metadata:
        print(f"Wrote {args.merged_metadata}")
    print("Group counts:", ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))


if __name__ == "__main__":
    main()
