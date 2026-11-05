#!/usr/bin/env python3
"""Collect every needs_review row from the *-hgnc-map.csv files into one review list.

Rows are grouped per source by (item, gene used, issue) so subject matter experts
see each distinct issue once, with how many rows it affects and an example.
CTann rows are left out: their flags only repeat ASCT+B issues (fix those at the
source). Each item gets a coarse `category` to help triage.

Usage: python3 build_review_items.py [-o output/review-items.csv] [--maps-dir output]
"""
import argparse
import csv
import re
import sys
from collections import OrderedDict
from pathlib import Path

# source name -> (map file, column naming the item, column with an example record)
SOURCES = {
    "hpa": ("hpa-hgnc-map.csv", "gene_symbol", "ensembl"),
    "deepcell-types": ("deepcell-types-hgnc-map.csv", "channel", "marker"),
    "avr": ("avr-hgnc-map.csv", "target_symbol", "antibody_id"),
    "asctb": ("asctb-hgnc-map.csv", "biomarker_label", "biomarker_record"),
    "omap": ("omap-hgnc-map.csv", "antibody_label", "used_antibody"),
}

# Checked in order; first match wins.
CATEGORIES = [
    ("wrong-id-likely", r"approved symbol of|truncated|disagrees with|suspected wrong target|"
                        r"antibody name mentions"),
    ("withdrawn-id", r"withdrawn|not a known HGNC ID"),
    ("ensembl-mismatch", r"Ensembl ID not in HGNC|Ensembl ID maps to several"),
    ("uniprot-mismatch", r"UniProt"),
    ("ambiguous-alias", r"ambiguous alias|alias symbol"),
    ("label-check", r"not a recognized name|empty label|previous symbol of|prev symbol of"),
    ("missing-target", r"no target gene"),
    ("non-gene-type", r"biomarker type is"),
    ("curated-mapping", r"."),  # DeepCell/HPA rows flagged by design
]

OUTPUT_COLUMNS = ["source", "category", "item", "hgnc_symbol", "suggested_hgnc_symbol",
                  "row_count", "example_record", "note"]


def issue_text(note):
    """Drop per-row free text (ASCT+B curator notes) so identical issues group together."""
    return "; ".join(p for p in note.split("; ") if not p.startswith("ASCT+B notes:"))


def categorize(note):
    for name, pattern in CATEGORIES:
        if re.search(pattern, note, re.IGNORECASE):
            return name
    return "other"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-o", "--output", type=Path, default=Path("output/review-items.csv"))
    ap.add_argument("--maps-dir", type=Path, default=Path("output"))
    args = ap.parse_args()

    items = OrderedDict()
    for source, (fname, item_col, example_col) in SOURCES.items():
        with open(args.maps_dir / fname, encoding="utf-8", newline="") as f:
            for r in csv.DictReader(f):
                if r.get("needs_review") != "true":
                    continue
                note = issue_text(r.get("note", ""))
                key = (source, r[item_col], r.get("hgnc_symbol", ""), note)
                if key not in items:
                    items[key] = {
                        "source": source,
                        "category": categorize(note),
                        "item": r[item_col],
                        "hgnc_symbol": r.get("hgnc_symbol", ""),
                        "suggested_hgnc_symbol": r.get("suggested_hgnc_symbol", ""),
                        "row_count": 0,
                        "example_record": r[example_col],
                        "note": note,
                    }
                items[key]["row_count"] += 1

    order = {name: i for i, (name, _) in enumerate(CATEGORIES)}
    rows = sorted(items.values(), key=lambda r: (list(SOURCES).index(r["source"]),
                                                 order.get(r["category"], 99), -r["row_count"]))
    with open(args.output, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, OUTPUT_COLUMNS)
        w.writeheader()
        w.writerows(rows)

    summary = {}
    for r in rows:
        summary.setdefault(r["source"], {}).setdefault(r["category"], 0)
        summary[r["source"]][r["category"]] += 1
    print(f"Wrote {len(rows)} review items to {args.output}: {summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
