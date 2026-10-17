#!/usr/bin/env python3
"""Replace S2's full-92 Fisher sheet with the exact merged Figure 3 tests.

The other S2 worksheets are copied byte-for-byte. The input table is produced
by make_figure3_lps_gene_enrichment.py. This script never overwrites inputs.
"""

import argparse
import csv
from pathlib import Path
from zipfile import ZipFile

from lxml import etree


NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
Q = lambda name: f"{{{NS}}}{name}"
SHEET = "xl/worksheets/sheet2.xml"
FIELDS = (
    "gene_id", "kegg_ko", "description", "iss_present", "iss_absent",
    "earth_present", "earth_absent", "p_value", "q_value",
)


def set_cell(row, column, number, value, style=None):
    cell = etree.SubElement(row, Q("c"), r=f"{column}{number}")
    if style is not None:
        cell.set("s", str(style))
    if column in "ABC":
        cell.set("t", "inlineStr")
        etree.SubElement(etree.SubElement(cell, Q("is")), Q("t")).text = value or ""
    else:
        etree.SubElement(cell, Q("v")).text = str(value)


def replace_sheet(xml, records):
    root = etree.fromstring(xml)
    data = root.find(Q("sheetData"))
    if data is None:
        raise ValueError("Missing worksheet data")
    rows = data.findall(Q("row"))
    if len(rows) != 52 or rows[3].get("r") != "5":
        raise ValueError("Unexpected S2 template: expected 48 tests below row 5")
    for row in rows[4:]:
        data.remove(row)

    def text_in_row(number, value):
        row = next(x for x in rows if x.get("r") == str(number))
        for cell in list(row):
            row.remove(cell)
        set_cell(row, "A", number, value, style=1)

    text_in_row(2, "40 ISS versus 52 Earth genomes; two-sided Fisher tests and BH correction across 46 nonredundant features after merging waaL, rfaL and KO:K02847 annotations.")
    text_in_row(3, "This sheet contains the exact tests plotted in Figure 3B; q values are unrounded numeric values.")
    for index, record in enumerate(records, start=6):
        row = etree.SubElement(data, Q("row"), r=str(index), spans="1:9")
        for column, field in zip("ABCDEFGHI", FIELDS):
            set_cell(row, column, index, record[field], style=1 if column == "C" else 4 if column in "HI" else None)
    end = 5 + len(records)
    root.find(Q("dimension")).set("ref", f"A1:I{end}")
    root.find(Q("autoFilter")).set("ref", f"A5:I{end}")
    return etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-workbook", type=Path, required=True)
    parser.add_argument("--figure3-table", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.resolve() == args.source_workbook.resolve():
        parser.error("Output already exists or would overwrite the source workbook")
    with args.figure3_table.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "gene_id", "kegg_ko", "description", "iss_present", "iss_absent",
            "earth_present", "earth_absent", "odds_ratio", "p_value", "q_value",
        ):
            parser.error("Unexpected Figure 3 table columns")
        records = list(reader)
    names = [row["gene_id"] for row in records]
    if len(records) != 46 or len(set(names)) != 46:
        parser.error("Expected 46 distinct Figure 3 tests")
    if names.count("waaL/rfaL (KO:K02847)") != 1 or any(alias in names for alias in ("waaL", "rfaL", "ko:K02847")):
        parser.error("Ligase labels were not merged as expected")
    rfb_c = next(row for row in records if row["gene_id"] == "rfbC")
    if tuple(rfb_c[field] for field in ("iss_present", "iss_absent", "earth_present", "earth_absent")) != ("39", "1", "9", "43"):
        parser.error("Unexpected rfbC contingency table")
    if not (4.34e-16 < float(rfb_c["p_value"]) < 4.36e-16):
        parser.error("rfbC p value is not the independently audited value")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(args.source_workbook) as source, ZipFile(args.output, "w") as target:
        for member in source.infolist():
            contents = source.read(member.filename)
            if member.filename == SHEET:
                contents = replace_sheet(contents, records)
            target.writestr(member, contents)
    with ZipFile(args.output) as check:
        if check.testzip() is not None:
            raise ValueError("Output workbook failed ZIP integrity test")
        root = etree.fromstring(check.read(SHEET))
        assert len(root.find(Q("sheetData")).findall(Q("row"))) == 50
    print(args.output)


if __name__ == "__main__":
    main()
