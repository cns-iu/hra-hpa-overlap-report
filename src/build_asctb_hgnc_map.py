#!/usr/bin/env python3
"""Map ASCT+B biomarker records in the HRA collection to HGNC and Ensembl gene IDs.

Reads every ccf:BiomarkerRecord (gene, protein, lipid, metabolite and proteoform
markers) from hra.ttl and keeps those whose source_concept references a gene:

* an HGNC IRI (http://identifiers.org/hgnc/N), resolved by HGNC ID, or
* an ASCTB-TEMP placeholder IRI whose label is a gene symbol or contains an
  Ensembl gene ID, resolved by approved symbol, then Ensembl ID, then previous
  symbol, then alias symbol.

Placeholder records that match no gene are left out (counted on stderr).
Writes one row per biomarker record.

HGNC IRIs are cross-checked against the record's label. When the label is the
exact approved symbol of a *different* gene (e.g. label SOX10 with HGNC:1119 =
BST2, a truncation of HGNC:11190), the label's gene is used instead
(hgnc_match = "label_symbol"); the original IRI stays in source_concept. All
such corrections, alias-based matches, unrecognized labels and non-gene
biomarker types that reference genes are flagged with needs_review = true.

Usage: python3 build_asctb_hgnc_map.py [-o OUTPUT.csv]
       [--hra downloads/hra.ttl] [--hgnc downloads/hgnc_complete_set.txt]
"""
import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from pyoxigraph import RdfFormat, Store

from hgnc_lookup import Hgnc, split_list

HGNC_IRI = "http://identifiers.org/hgnc/"
TEMP_IRI = "https://purl.org/ccf/ASCTB-TEMP_"

QUERY = """
PREFIX ccf: <http://purl.org/ccf/>
SELECT ?biomarker ?type ?label ?concept ?notes WHERE {
  ?biomarker a ccf:BiomarkerRecord ;
             ccf:ccf_biomarker_type ?type ;
             ccf:ccf_pref_label ?label ;
             ccf:source_concept ?concept .
  OPTIONAL { ?biomarker ccf:notes ?notes }
}
"""

RECORD_QUERY = """
PREFIX ccf: <http://purl.org/ccf/>
SELECT ?record ?biomarker WHERE {
  VALUES ?p { ccf:gene_marker ccf:protein_marker ccf:lipid_marker ccf:metabolite_marker ccf:proteoform_marker }
  ?record ?p ?biomarker .
}
"""

OUTPUT_COLUMNS = [
    "asctb_table", "asctb_version", "record_number", "column", "biomarker_type",
    "biomarker_label", "source_concept", "ensembl", "hgnc_id", "hgnc_symbol",
    "hgnc_match", "uniprot", "gene_description", "asctb_record", "biomarker_record",
    "needs_review", "suggested_hgnc_symbol", "note",
]

GENE_TYPES = {"gene", "protein", "proteoforms"}


def suggestion(hgnc, concept, label):
    """Likely intended gene(s) for a flagged HGNC IRI whose ID looks truncated."""
    if not concept.startswith(HGNC_IRI):
        return ""
    trunc = hgnc.truncation_candidates("HGNC:" + concept[len(HGNC_IRI):], label)
    return "|".join(x["symbol"] for x in trunc)


def check_label(hgnc, hgnc_id, h, label):
    """Compare an HGNC IRI's gene with the record label -> (hits, match, notes, review)."""
    if hgnc.is_known_name(label, h):
        return [h], "hgnc_id", [], False
    other = hgnc.by_symbol.get(label.strip().upper())
    if other:
        note = (f"source_concept {hgnc_id} is {h['symbol']}, but label {label} is the approved "
                f"symbol of {other['hgnc_id']}; using {other['symbol']}")
        if other["hgnc_id"].startswith(hgnc_id):
            note += f" ({hgnc_id} looks like a truncated {other['hgnc_id']})"
        return [other], "label_symbol", [note], True
    hits, how = hgnc.match_label(label)
    if hits:
        names = "|".join(f"{x['symbol']} ({x['hgnc_id']})" for x in hits)
        return [h], "hgnc_id", [f"source_concept {hgnc_id} is {h['symbol']}, but label {label} is "
                                f"an {how.replace('_', ' ')} of {names}; keeping {h['symbol']}"], True
    if not label.strip():
        return [h], "hgnc_id", ["biomarker has an empty label"], True
    note = f"label {label} is not a recognized name of {h['symbol']} ({hgnc_id}); keeping {h['symbol']}"
    trunc = hgnc.truncation_candidates(hgnc_id, label)
    if trunc:
        note += "; " + hgnc_id + " may be a truncated ID for " + \
            "|".join(f"{x['symbol']} ({x['hgnc_id']})" for x in trunc)
    return [h], "hgnc_id", [note], True


def parse_biomarker_iri(iri):
    # https://purl.humanatlas.io/asct-b/<table>/<version>#R<record>-<column><order>
    base, _, frag = iri.partition("#")
    table, version = base.rstrip("/").split("/")[-2:]
    m = re.fullmatch(r"R(\d+)-([A-Z]+)(\d+)", frag)
    record, column = (m.group(1), f"{m.group(2)}/{m.group(3)}") if m else ("", "")
    return table, version, record, column


def build_rows(store, hgnc):
    skipped = Counter()
    records = defaultdict(list)
    for sol in store.query(RECORD_QUERY):
        records[sol["biomarker"].value].append(sol["record"].value)
    for sol in store.query(QUERY):
        concept = sol["concept"].value
        label = sol["label"].value
        btype = sol["type"].value
        notes, review = [], False
        if concept.startswith(HGNC_IRI):
            hgnc_id = "HGNC:" + concept[len(HGNC_IRI):]
            h = hgnc.by_id.get(hgnc_id)
            if h is None:
                review = True
                notes.append(hgnc.describe_missing(hgnc_id))
                hits, match = hgnc.match_label(label)
                if hits:
                    notes.append(f"matched by label ({match})")
            else:
                hits, match, label_notes, review = check_label(hgnc, hgnc_id, h, label)
                notes += label_notes
        elif concept.startswith(TEMP_IRI):
            hits, match = hgnc.match_label(label)
            if not hits:
                skipped[btype] += 1
                continue
            notes.append(f"ASCT+B has no ID for this biomarker (placeholder IRI); matched label "
                         f"by {match.replace('_', ' ')}")
            if match.startswith("alias_symbol"):
                review = True
                if len(hits) > 1:
                    notes.append("ambiguous alias: an expert should pick the right gene")
        else:
            skipped["other source_concept"] += 1
            continue
        if hits and btype not in GENE_TYPES:
            review = True
            notes.append(f"biomarker type is {btype} but it resolves to a gene")

        biomarker = sol["biomarker"].value
        table, version, record, column = parse_biomarker_iri(biomarker)
        if any(not h["ensembl_gene_id"] for h in hits):
            notes.append("HGNC record has no Ensembl gene ID (not a single gene)")
        if sol["notes"] is not None and sol["notes"].value:
            notes.append("ASCT+B notes: " + sol["notes"].value)
        yield {
            "asctb_table": table,
            "asctb_version": version,
            "record_number": record,
            "column": column,
            "biomarker_type": btype,
            "biomarker_label": label,
            "source_concept": concept,
            "ensembl": "|".join(h["ensembl_gene_id"] for h in hits),
            "hgnc_id": "|".join(h["hgnc_id"] for h in hits),
            "hgnc_symbol": "|".join(h["symbol"] for h in hits),
            "hgnc_match": match,
            "uniprot": "|".join(h["uniprot_ids"] for h in hits),
            "gene_description": "|".join(h["name"] for h in hits),
            "asctb_record": "|".join(sorted(records.get(biomarker, []))),
            "biomarker_record": biomarker,
            "needs_review": str(review).lower(),
            "suggested_hgnc_symbol": suggestion(hgnc, concept, label) if match == "hgnc_id" else "",
            "note": "; ".join(notes),
        }
    build_rows.skipped = skipped


def sort_key(r):
    return (r["asctb_table"], int(r["record_number"] or 0), r["column"].split("/")[0],
            int(r["column"].split("/")[-1] or 0) if r["column"] else 0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-o", "--output", type=Path, default=Path("output/asctb-hgnc-map.csv"))
    ap.add_argument("--hra", type=Path, default=Path("downloads/hra.ttl"))
    ap.add_argument("--hgnc", type=Path, default=Path("downloads/hgnc_complete_set.txt"))
    args = ap.parse_args()

    store = Store()
    store.bulk_load(path=str(args.hra), format=RdfFormat.TURTLE)
    hgnc = Hgnc(args.hgnc)
    rows = sorted(build_rows(store, hgnc), key=sort_key)
    with open(args.output, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, OUTPUT_COLUMNS)
        w.writeheader()
        w.writerows(rows)

    matches = Counter(r["hgnc_match"] for r in rows)
    genes = len({g for r in rows for g in split_list(r["hgnc_id"])})
    skipped = build_rows.skipped
    review = sum(r["needs_review"] == "true" for r in rows)
    print(f"Wrote {len(rows)} biomarker records ({genes} distinct HGNC genes) to {args.output}; "
          f"matches: {dict(matches)}; {review} need review; "
          f"skipped (no gene match): {dict(skipped)}", file=sys.stderr)


if __name__ == "__main__":
    main()
