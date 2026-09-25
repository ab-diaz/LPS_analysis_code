#!/usr/bin/env python3
"""Validate a protein FASTA alignment and write interleaved NEXUS."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


VALID = set("ABCDEFGHIKLMNPQRSTVWXYZ-?")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    with path.open() as handle:
        for line_no, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts).upper()))
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"Empty FASTA identifier at line {line_no}")
                parts = []
            elif name is None:
                raise ValueError(f"Sequence before first FASTA header at line {line_no}")
            else:
                parts.append(line)
    if name is not None:
        records.append((name, "".join(parts).upper()))
    if not records:
        raise ValueError("No FASTA records found")
    return records


def nexus_quote(text: str) -> str:
    return "'" + text.replace("'", "''") + "'"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--width", type=int, default=80)
    args = parser.parse_args()

    records = read_fasta(args.input)
    names = [name for name, _ in records]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate taxon identifiers in FASTA")
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError(f"Alignment is ragged: observed lengths {sorted(lengths)}")
    nchar = lengths.pop()
    invalid = sorted({char for _, seq in records for char in seq if char not in VALID})
    if invalid:
        raise ValueError(f"Invalid protein characters: {invalid}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as out:
        out.write("#NEXUS\n\nBEGIN DATA;\n")
        out.write(f"  DIMENSIONS NTAX={len(records)} NCHAR={nchar};\n")
        out.write("  FORMAT DATATYPE=PROTEIN MISSING=? GAP=- INTERLEAVE=YES;\n")
        out.write("  MATRIX\n")
        for start in range(0, nchar, args.width):
            stop = min(nchar, start + args.width)
            for name, seq in records:
                out.write(f"    {nexus_quote(name)} {seq[start:stop]}\n")
            out.write("\n")
        out.write("  ;\nEND;\n")

    metadata = {
        "source": str(args.input),
        "source_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "output": str(args.output),
        "taxa": len(records),
        "sites": nchar,
        "datatype": "protein",
        "gap_character": "-",
        "missing_character": "?",
        "taxon_names": names,
    }
    args.metadata.parent.mkdir(parents=True, exist_ok=True)
    args.metadata.write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
