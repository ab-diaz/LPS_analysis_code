#!/usr/bin/env python3
"""Correct the verified F8/F9 typo in a rendered Supplementary Figure S1 PDF."""

import argparse
from pathlib import Path
import fitz

WRONG = "Pantoea_piersonii_F8_6S_D2_EB4_spades"
RIGHT = "Pantoea_piersonii_F9_6S_D2_EB4_spades"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input_pdf", type=Path)
    ap.add_argument("output_pdf", type=Path)
    args = ap.parse_args()
    if args.output_pdf.exists():
        raise SystemExit(f"Refusing to overwrite {args.output_pdf}")
    doc = fitz.open(args.input_pdf)
    if len(doc) != 1:
        raise SystemExit("Expected a one-page network figure")
    page = doc[0]
    spans = [span for block in page.get_text("dict")["blocks"] if "lines" in block
             for line in block["lines"] for span in line["spans"] if WRONG in span["text"]]
    if len(spans) != 1 or spans[0]["text"] != WRONG:
        raise SystemExit("Expected the old label in exactly one isolated text span")
    span = spans[0]
    page.add_redact_annot(fitz.Rect(span["bbox"]), fill=(1, 1, 1))
    page.apply_redactions(images=0, graphics=0)
    page.insert_text(span["origin"], RIGHT, fontsize=span["size"], fontname="helv",
                     color=(0, 0, 0), overlay=True)
    if WRONG in page.get_text() or RIGHT not in page.get_text():
        raise SystemExit("Label correction failed validation")
    doc.save(args.output_pdf, garbage=4, deflate=True)


if __name__ == "__main__":
    main()
