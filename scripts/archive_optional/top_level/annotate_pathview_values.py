#!/usr/bin/env python3
import argparse
import csv
import os
import re
import xml.etree.ElementTree as ET

from PIL import Image, ImageDraw, ImageFont


KO_RE = re.compile(r"K\d{5}")


def load_ko_values(path, shared_only=False):
    values = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ko = row["ko"].strip()
            earth_count = int(float(row["Earth_count"]))
            iss_count = int(float(row["ISS_count"]))
            if shared_only and (earth_count == 0 or iss_count == 0):
                continue
            values[ko] = float(row["ISS_frac"]) - float(row["Earth_frac"])
    return values


def load_font(size, bold=False):
    candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
    ]
    for path in candidates:
        if os.path.exists(path):
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


def extract_entries(kgml_path, values):
    tree = ET.parse(kgml_path)
    root = tree.getroot()
    entries = []
    for entry in root.findall("entry"):
        if entry.get("type") != "ortholog":
            continue
        graphics = entry.find("graphics")
        if graphics is None or graphics.get("type") != "rectangle":
            continue

        kos = set(KO_RE.findall(entry.get("name", "")))
        kos.update(KO_RE.findall(graphics.get("name", "")))
        matched = sorted(ko for ko in kos if ko in values)
        if not matched:
            continue

        # Pathview's default node.sum is "sum", so multiple KOs on the same
        # KEGG node are collapsed before coloring the rectangle.
        node_value = sum(values[ko] for ko in matched)
        label = f"{node_value:+.2f}"

        entries.append(
            {
                "x": int(float(graphics.get("x"))),
                "y": int(float(graphics.get("y"))),
                "w": int(float(graphics.get("width"))),
                "h": int(float(graphics.get("height"))),
                "label": label,
                "value": node_value,
            }
        )
    return entries


def label_color(value):
    if value > 0.001:
        return (160, 0, 0)
    if value < -0.001:
        return (0, 115, 0)
    return (45, 45, 45)


def draw_centered_label(draw, xy, label, font, fill):
    x, y = xy
    bbox = draw.textbbox((0, 0), label, font=font)
    width = bbox[2] - bbox[0]
    height = bbox[3] - bbox[1]
    pad_x = 2
    pad_y = 1
    rect = [
        x - width / 2 - pad_x,
        y - height / 2 - pad_y,
        x + width / 2 + pad_x,
        y + height / 2 + pad_y,
    ]
    draw.rectangle(rect, fill=(255, 255, 255, 230))
    draw.text((x - width / 2, y - height / 2 - 1), label, font=font, fill=fill)


def annotate(image_path, kgml_path, table_path, output_path, shared_only=False, font_size=12):
    values = load_ko_values(table_path, shared_only=shared_only)
    entries = extract_entries(kgml_path, values)

    image = Image.open(image_path).convert("RGBA")
    overlay = Image.new("RGBA", image.size, (255, 255, 255, 0))
    draw = ImageDraw.Draw(overlay)
    font = load_font(font_size, bold=True)

    for entry in entries:
        x = entry["x"]
        y = max(8, entry["y"] - entry["h"] / 2 - font_size / 2 - 4)
        draw_centered_label(draw, (x, y), entry["label"], font, label_color(entry["value"]))

    image = Image.alpha_composite(image, overlay).convert("RGB")
    image.save(output_path)
    return len(entries)


def main():
    parser = argparse.ArgumentParser(
        description="Add numeric ISS-minus-Earth node values above KO boxes in Pathview PNGs."
    )
    parser.add_argument("pathways", nargs="+", help="Pathway IDs such as ko00520")
    parser.add_argument("--image-template", default="{pathway}.ISS_minus_Earth.png")
    parser.add_argument("--kgml-template", default="{pathway}.kgml")
    parser.add_argument("--table-template", default="lps_analysis/pathview_tables/{pathway}_ko.tsv")
    parser.add_argument("--out-template", default="{pathway}.ISS_minus_Earth.values.png")
    parser.add_argument("--font-size", type=int, default=12)
    parser.add_argument(
        "--shared-only",
        action="store_true",
        help="Only include KOs present in both Earth and ISS before summing node values",
    )
    args = parser.parse_args()

    for pathway in args.pathways:
        image_path = args.image_template.format(pathway=pathway)
        kgml_path = args.kgml_template.format(pathway=pathway)
        table_path = args.table_template.format(pathway=pathway)
        output_path = args.out_template.format(pathway=pathway)
        count = annotate(
            image_path,
            kgml_path,
            table_path,
            output_path,
            shared_only=args.shared_only,
            font_size=args.font_size,
        )
        print(f"{output_path}\t{count} labels")


if __name__ == "__main__":
    main()
