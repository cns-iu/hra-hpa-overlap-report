#!/usr/bin/env python3
"""Map Antibody Validation Report (AVR) antibodies to HGNC and Ensembl gene IDs.

Reads the AVR graph (avr.ttl) and, for every antibody, resolves its HGNC ID in
the HGNC complete set to get the approved symbol, Ensembl gene ID and gene
name. The AVR target symbol and UniProt accession are cross-checked against
the HGNC record and any disagreement is reported in the `note` column.

Some AVR records have a target_symbol that disagrees with their hgnc_id. When
that happens the AVR PDF filename (e.g. "CD3E_OMAP9.pdf") is used as the tie
breaker: if it names the target symbol, the target symbol's HGNC record is used
(hgnc_match = "target_symbol"); otherwise the HGNC ID is trusted
(hgnc_match = "hgnc_id").

AVR records lack antibody names, so wrong targets that the AVR PDF filename
repeats (e.g. a CD123 antibody filed as CDC123_OMAP17.pdf) cannot be caught
here. When --omap points at omap-hgnc-map.csv, AVR rows for an (RRID, gene)
that the OMAP map flags inherit that flag and its suggested_hgnc_symbol.

Usage: python3 build_avr_hgnc_map.py [-o OUTPUT.csv]
       [--avr downloads/avr.ttl] [--hgnc downloads/hgnc_complete_set.txt]
       [--omap output/omap-hgnc-map.csv]
"""
import argparse
import csv
import re
import sys
from pathlib import Path

from rdflib import RDF, Graph, URIRef

from hgnc_lookup import Hgnc, split_list

AVR = "https://purl.humanatlas.io/schema/avr#"
FIELDS = [
    "antibody_uuid", "target_symbol", "hgnc_id", "uniprot_accession_number", "rrid",
    "clone_id", "catalog_number", "vendor_name", "host", "clonality", "conjugate",
    "method", "organ", "omap", "avr_pdf_filename", "recommended", "previous_version_id",
    "next_version_id",
]

OUTPUT_COLUMNS = [
    "antibody_id", "rrid", "target_symbol", "ensembl", "hgnc_id", "hgnc_symbol",
    "hgnc_match", "avr_hgnc_id", "uniprot", "gene_description", "clone_id",
    "catalog_number", "vendor_name", "host", "clonality", "conjugate", "method",
    "organ", "omap", "recommended", "previous_version_id", "next_version_id",
    "avr_pdf_filename", "needs_review", "suggested_hgnc_symbol", "note",
]


def read_omap_flags(path):
    """(RRID without prefix, gene symbol) -> (suggested symbol, OMAP note) for flagged OMAP rows."""
    flags = {}
    if path is None or not Path(path).exists():
        return flags
    with open(path, encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            if r.get("needs_review") == "true" and r["hgnc_symbol"]:
                flags.setdefault((r["rrid"], r["hgnc_symbol"]), (r["suggested_hgnc_symbol"], r["note"]))
    return flags


def local_id(uri, sep="/"):
    return str(uri).rsplit(sep, 1)[-1] if uri is not None else ""


def read_antibodies(path):
    g = Graph()
    g.parse(path, format="turtle")
    for ab in g.subjects(RDF.type, URIRef(AVR + "Antibody")):
        rec = {f: g.value(ab, URIRef(AVR + f)) for f in FIELDS}
        rec["antibody_id"] = local_id(ab, "#")
        yield rec


def filename_tokens(name):
    return {t.upper() for t in re.split(r"[-_. ]+", name) if t}


def build_row(rec, hgnc, omap_flags=None):
    target = str(rec["target_symbol"] or "")
    avr_hgnc_id = "HGNC:" + local_id(rec["hgnc_id"]) if rec["hgnc_id"] is not None else ""
    uniprot = local_id(rec["uniprot_accession_number"])
    pdf = str(rec["avr_pdf_filename"] or "")
    target_rec = hgnc.by_symbol.get(target.upper())
    h = hgnc.by_id.get(avr_hgnc_id)
    match = "hgnc_id"
    notes, review = [], False
    if h is None:
        notes.append("AVR " + hgnc.describe_missing(avr_hgnc_id) if avr_hgnc_id
                     else "AVR record has no HGNC ID")
        review = True
        if target_rec:
            h, match = target_rec, "target_symbol"
            notes.append(f"using target symbol {target}")
    elif target.upper() != h["symbol"].upper():
        review = True
        if hgnc.is_known_name(target, h):
            notes.append(f"AVR target symbol {target} is an outdated/alias symbol of {h['symbol']}")
            review = False
        elif target_rec and target.upper() in filename_tokens(pdf):
            notes.append(f"AVR HGNC ID {avr_hgnc_id} ({h['symbol']}) disagrees with target "
                         f"symbol {target}; the AVR PDF filename names {target}, so using "
                         f"{target}. AVR record should be corrected")
            h, match = target_rec, "target_symbol"
        else:
            notes.append(f"AVR target symbol {target} disagrees with HGNC ID {avr_hgnc_id}; "
                         f"the AVR PDF filename ({pdf}) does not name {target}, so using the HGNC "
                         f"ID ({h['symbol']}). AVR record should be corrected")
    if h is None:
        match = ""
    else:
        if not h["ensembl_gene_id"]:
            notes.append("HGNC record has no Ensembl gene ID (not a single gene)")
        base_uniprot = uniprot.split("-")[0]
        if base_uniprot and base_uniprot not in split_list(h["uniprot_ids"]):
            notes.append(f"AVR UniProt {uniprot} is not {h['symbol']}'s UniProt "
                         f"({h['uniprot_ids'] or 'none'})")
            review = True
    suggested = ""
    rrid = local_id(rec["rrid"]).removeprefix("RRID:")
    omap_flag = (omap_flags or {}).get((rrid, h["symbol"])) if h else None
    if omap_flag:
        suggested = omap_flag[0]
        notes.append(f"the OMAP map flags this antibody: {omap_flag[1]}")
        review = True
    recommended = "" if rec["recommended"] is None else str(rec["recommended"]).lower()
    if recommended == "false":
        notes.append("AVR marks this antibody as not recommended")
    h = h or {}
    return {
        "antibody_id": rec["antibody_id"],
        "rrid": local_id(rec["rrid"]),
        "target_symbol": target,
        "ensembl": h.get("ensembl_gene_id", ""),
        "hgnc_id": h.get("hgnc_id", ""),
        "hgnc_symbol": h.get("symbol", ""),
        "hgnc_match": match,
        "avr_hgnc_id": avr_hgnc_id,
        "uniprot": uniprot,
        "gene_description": h.get("name", ""),
        **{f: "" if rec[f] is None else str(rec[f]) for f in (
            "clone_id", "catalog_number", "vendor_name", "host", "clonality",
            "conjugate", "method", "organ", "omap", "previous_version_id",
            "next_version_id", "avr_pdf_filename")},
        "recommended": recommended,
        "needs_review": str(review).lower(),
        "suggested_hgnc_symbol": suggested,
        "note": "; ".join(notes),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-o", "--output", type=Path, default=Path("output/avr-hgnc-map.csv"))
    ap.add_argument("--avr", type=Path, default=Path("downloads/avr.ttl"))
    ap.add_argument("--hgnc", type=Path, default=Path("downloads/hgnc_complete_set.txt"))
    ap.add_argument("--omap", type=Path, default=Path("output/omap-hgnc-map.csv"),
                    help="OMAP map whose antibody-level flags AVR rows inherit (skipped if missing)")
    args = ap.parse_args()

    hgnc = Hgnc(args.hgnc)
    omap_flags = read_omap_flags(args.omap)
    rows = [build_row(r, hgnc, omap_flags) for r in read_antibodies(args.avr)]
    rows.sort(key=lambda r: (r["target_symbol"].upper(), r["omap"], r["antibody_id"]))
    with open(args.output, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, OUTPUT_COLUMNS)
        w.writeheader()
        w.writerows(rows)

    matched = sum(bool(r["hgnc_id"]) for r in rows)
    by_symbol = sum(r["hgnc_match"] == "target_symbol" for r in rows)
    flagged = sum(r["needs_review"] == "true" for r in rows)
    targets = len({r["hgnc_id"] for r in rows if r["hgnc_id"]})
    print(f"Wrote {len(rows)} antibodies ({targets} distinct HGNC genes) to {args.output}: "
          f"{matched} mapped to approved HGNC genes ({by_symbol} via target symbol "
          f"instead of AVR HGNC ID), {flagged} need review", file=sys.stderr)


if __name__ == "__main__":
    main()
