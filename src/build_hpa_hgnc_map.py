#!/usr/bin/env python3
"""Join the Human Protein Atlas gene table with HGNC IDs and write a CSV.

Downloads proteinatlas.tsv.zip from the HPA and the HGNC complete set, then
matches each HPA gene to an Approved HGNC record by Ensembl gene ID. When the
Ensembl ID is not in HGNC, falls back to the HPA gene name as an approved HGNC
symbol (hgnc_match = "symbol", flagged for review because the Ensembl IDs differ).

Usage: python3 build_hpa_hgnc_map.py [-o OUTPUT.csv] [--download-dir DIR] [--refresh]
"""
import argparse
import csv
import io
import sys
import urllib.request
import zipfile
from collections import defaultdict
from pathlib import Path

HPA_URL = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

OUTPUT_COLUMNS = [
    "ensembl", "gene_symbol", "hgnc_id", "hgnc_symbol", "hgnc_match",
    "uniprot", "gene_description", "protein_class", "evidence", "needs_review", "note",
]


def download(url, dest, refresh):
    if dest.exists() and not refresh:
        print(f"Using cached {dest}", file=sys.stderr)
        return dest
    print(f"Downloading {url}", file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as resp, open(tmp, "wb") as f:
        while chunk := resp.read(1 << 20):
            f.write(chunk)
    tmp.replace(dest)
    return dest


def read_hpa(zip_path):
    with zipfile.ZipFile(zip_path) as zf:
        name = next(n for n in zf.namelist() if n.endswith(".tsv"))
        with zf.open(name) as raw:
            yield from csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8"), delimiter="\t")


def read_hgnc(path):
    by_ensembl, by_symbol = defaultdict(list), {}
    with open(path, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["status"] != "Approved":
                continue
            by_symbol[r["symbol"]] = r
            if r["ensembl_gene_id"]:
                by_ensembl[r["ensembl_gene_id"]].append(r)
    return by_ensembl, by_symbol


def build_rows(hpa_rows, hgnc_by_ensembl, hgnc_by_symbol):
    for r in hpa_rows:
        hits = hgnc_by_ensembl.get(r["Ensembl"], [])
        note, review = "", False
        if len(hits) > 1:
            match = "ensembl_multiple"
            note, review = "Ensembl ID maps to several HGNC genes", True
        elif hits:
            match = "ensembl"
            if hits[0]["symbol"] != r["Gene"]:
                note = f"HPA gene name {r['Gene']} differs from current HGNC symbol {hits[0]['symbol']}"
        elif r["Gene"] in hgnc_by_symbol:
            hits, match = [hgnc_by_symbol[r["Gene"]]], "symbol"
            other = hits[0]["ensembl_gene_id"] or "none"
            note = (f"HPA Ensembl ID not in HGNC; matched by symbol {r['Gene']}, "
                    f"whose HGNC Ensembl ID is {other}")
            review = True
        else:
            match = ""
            note = "no HGNC record for this Ensembl ID or gene name"
        yield {
            "ensembl": r["Ensembl"],
            "gene_symbol": r["Gene"],
            "hgnc_id": "|".join(h["hgnc_id"] for h in hits),
            "hgnc_symbol": "|".join(h["symbol"] for h in hits),
            "hgnc_match": match,
            "uniprot": r["Uniprot"],
            "gene_description": r["Gene description"],
            "protein_class": r["Protein class"],
            "evidence": r["Evidence"],
            "needs_review": str(review).lower(),
            "note": note,
        }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-o", "--output", type=Path, default=Path("output/hpa-hgnc-map.csv"))
    ap.add_argument("--download-dir", type=Path, default=Path("downloads"))
    ap.add_argument("--refresh", action="store_true", help="re-download even if cached files exist")
    args = ap.parse_args()

    args.download_dir.mkdir(parents=True, exist_ok=True)
    hpa_zip = download(HPA_URL, args.download_dir / "proteinatlas.tsv.zip", args.refresh)
    hgnc_txt = download(HGNC_URL, args.download_dir / "hgnc_complete_set.txt", args.refresh)

    hgnc_by_ensembl, hgnc_by_symbol = read_hgnc(hgnc_txt)
    total = matched = renamed = 0
    with open(args.output, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, OUTPUT_COLUMNS)
        w.writeheader()
        for row in build_rows(read_hpa(hpa_zip), hgnc_by_ensembl, hgnc_by_symbol):
            w.writerow(row)
            total += 1
            if row["hgnc_id"]:
                matched += 1
                if row["hgnc_symbol"] != row["gene_symbol"]:
                    renamed += 1

    print(f"Wrote {total} genes to {args.output}: {matched} with HGNC IDs, "
          f"{total - matched} without, {renamed} with a different HGNC symbol", file=sys.stderr)


if __name__ == "__main__":
    main()
