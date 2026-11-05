#!/usr/bin/env python3
"""Infer biomarker genes for CTann (cell type annotation) mappings via ASCT+B tables.

CTann mappings in hra.ttl only link an annotation tool's cell type label to a
Cell Ontology (CL) term. This script joins each mapping's CL term to the ASCT+B
records whose most specific cell type (the last CT column, which the record's
biomarkers characterize) is that same CL term, and collects those records'
gene-resolvable biomarkers (resolved as in build_asctb_hgnc_map.py).

Writes one row per (CTann mapping, gene). `same_organ` is true when at least
one supporting ASCT+B record lists the CTann mapping's organ among its
anatomical structures. Mappings with no inferred genes are kept with empty
gene columns.

Every gene here is inferred, not asserted by the CTann. needs_review is true
when *all* ASCT+B biomarker rows supporting a gene are themselves flagged for
review in asctb-hgnc-map.csv (see that file's notes for the reasons).

Usage: python3 build_ctann_hgnc_map.py [-o OUTPUT.csv]
       [--hra downloads/hra.ttl] [--hgnc downloads/hgnc_complete_set.txt]
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

from pyoxigraph import RdfFormat, Store

from build_asctb_hgnc_map import build_rows as build_asctb_rows
from hgnc_lookup import Hgnc

OBO = "http://purl.obolibrary.org/obo/"
PREFIXES = """
PREFIX ccf: <http://purl.org/ccf/>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX sssom: <https://w3id.org/sssom/>
"""

CTANN_QUERY = PREFIXES + """
SELECT ?mapping ?subject ?subject_label ?predicate ?object ?object_label ?organ ?organ_level WHERE {
  ?mapping a ccf:CtAnnMapping ;
           owl:annotatedSource ?subject ;
           owl:annotatedTarget ?object .
  OPTIONAL { ?mapping sssom:subject_label ?subject_label }
  OPTIONAL { ?mapping owl:annotatedProperty ?predicate }
  OPTIONAL { ?mapping sssom:object_label ?object_label }
  OPTIONAL { ?mapping ccf:ext_organ_id ?organ }
  OPTIONAL { ?mapping ccf:ext_organ_level ?organ_level }
}
"""

CELL_TYPE_QUERY = PREFIXES + """
SELECT ?record ?order ?concept WHERE {
  ?record a ccf:AsctbRecord ; ccf:cell_type ?ct .
  ?ct ccf:order_number ?order ; ccf:source_concept ?concept .
}
"""

ANATOMY_QUERY = PREFIXES + """
SELECT ?record ?concept WHERE {
  ?record a ccf:AsctbRecord ; ccf:anatomical_structure ?as .
  ?as ccf:source_concept ?concept .
}
"""

OUTPUT_COLUMNS = [
    "ctann", "ctann_version", "mapping_record", "subject_id", "subject_label", "predicate",
    "cl_id", "cl_label", "organ_level", "organ_uberon_id", "ensembl", "hgnc_id",
    "hgnc_symbol", "hgnc_match", "uniprot", "gene_description", "biomarker_types",
    "same_organ", "asctb_record_count", "asctb_tables", "needs_review", "note",
]


def val(term):
    return "" if term is None else term.value


def curie(iri):
    return iri[len(OBO):].replace("_", ":", 1) if iri.startswith(OBO) else iri


def parse_dataset_iri(iri):
    # https://purl.humanatlas.io/<type>/<name>/<version>#<fragment>
    base, _, frag = iri.partition("#")
    name, version = base.rstrip("/").split("/")[-2:]
    return name, version, frag


def primary_cell_types(store):
    """AsctbRecord IRI -> CURIE of its most specific (highest order) cell type."""
    best = {}
    for s in store.query(CELL_TYPE_QUERY):
        record, order = val(s["record"]), int(val(s["order"]))
        if record not in best or order > best[record][0]:
            best[record] = (order, curie(val(s["concept"])))
    return {r: c for r, (_, c) in best.items()}


def record_anatomy(store):
    anatomy = defaultdict(set)
    for s in store.query(ANATOMY_QUERY):
        anatomy[val(s["record"])].add(curie(val(s["concept"])))
    return anatomy


def genes_by_cell_type(store, hgnc):
    """CL CURIE -> {hgnc_id: {"gene": row, "types": set, "records": set}}."""
    primary = primary_cell_types(store)
    by_ct = defaultdict(dict)
    for row in build_asctb_rows(store, hgnc):
        if not row["hgnc_id"]:
            continue
        for record in row["asctb_record"].split("|"):
            ct = primary.get(record)
            if not ct:
                continue
            g = by_ct[ct].setdefault(row["hgnc_id"], {"gene": row, "types": set(),
                                                      "matches": set(), "records": set(),
                                                      "rows": 0, "flagged": 0})
            g["rows"] += 1
            g["flagged"] += row["needs_review"] == "true"
            g["types"].add(row["biomarker_type"])
            g["matches"].add(row["hgnc_match"])
            g["records"].add(record)
    return by_ct


def build_rows(store, hgnc):
    by_ct = genes_by_cell_type(store, hgnc)
    anatomy = record_anatomy(store)
    stats = {"mappings": 0, "with_genes": 0}
    for s in store.query(CTANN_QUERY):
        mapping = val(s["mapping"])
        name, version, frag = parse_dataset_iri(mapping)
        cl_id = val(s["object"])
        organ = curie(val(s["organ"]))
        base = {
            "ctann": name,
            "ctann_version": version,
            "mapping_record": frag,
            "subject_id": val(s["subject"]),
            "subject_label": val(s["subject_label"]),
            "predicate": val(s["predicate"]),
            "cl_id": cl_id,
            "cl_label": val(s["object_label"]),
            "organ_level": val(s["organ_level"]),
            "organ_uberon_id": organ,
        }
        genes = by_ct.get(cl_id, {})
        stats["mappings"] += 1
        stats["with_genes"] += bool(genes)
        if not genes:
            yield {**base, **{c: "" for c in OUTPUT_COLUMNS if c not in base},
                   "needs_review": "false",
                   "note": "no ASCT+B record has this CL term as its most specific cell type"}
            continue
        for hgnc_id, g in genes.items():
            gene = g["gene"]
            yield {
                **base,
                "ensembl": gene["ensembl"],
                "hgnc_id": hgnc_id,
                "hgnc_symbol": gene["hgnc_symbol"],
                "hgnc_match": "|".join(sorted(g["matches"])),
                "uniprot": gene["uniprot"],
                "gene_description": gene["gene_description"],
                "biomarker_types": "|".join(sorted(g["types"])),
                "same_organ": str(any(organ in anatomy[r] for r in g["records"])).lower(),
                "asctb_record_count": len(g["records"]),
                "asctb_tables": "|".join(sorted({parse_dataset_iri(r)[0] for r in g["records"]})),
                "needs_review": str(g["flagged"] == g["rows"]).lower(),
                "note": (f"{g['flagged']} of {g['rows']} supporting ASCT+B biomarker rows are "
                         f"flagged for review" if g["flagged"] else ""),
            }
    build_rows.stats = stats


def sort_key(r):
    num = r["mapping_record"].lstrip("R")
    return (r["ctann"], int(num) if num.isdigit() else 0, r["mapping_record"],
            -int(r["asctb_record_count"] or 0), r["hgnc_symbol"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-o", "--output", type=Path, default=Path("output/ctann-hgnc-map.csv"))
    ap.add_argument("--hra", type=Path, default=Path("downloads/hra.ttl"))
    ap.add_argument("--hgnc", type=Path, default=Path("downloads/hgnc_complete_set.txt"))
    args = ap.parse_args()

    store = Store()
    store.bulk_load(path=str(args.hra), format=RdfFormat.TURTLE)
    rows = sorted(build_rows(store, Hgnc(args.hgnc)), key=sort_key)
    with open(args.output, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, OUTPUT_COLUMNS)
        w.writeheader()
        w.writerows(rows)

    st = build_rows.stats
    genes = len({r["hgnc_id"] for r in rows if r["hgnc_id"]})
    same = sum(r["same_organ"] == "true" for r in rows)
    print(f"Wrote {len(rows)} rows to {args.output}: {st['with_genes']} of {st['mappings']} CTann "
          f"mappings have inferred genes ({genes} distinct HGNC genes); "
          f"{same} gene rows supported by a same-organ ASCT+B record; "
          f"{sum(r['needs_review'] == 'true' for r in rows)} need review", file=sys.stderr)


if __name__ == "__main__":
    main()
