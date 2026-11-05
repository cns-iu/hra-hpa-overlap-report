#!/usr/bin/env python3
"""Map OMAP (Organ Mapping Antibody Panel) antibodies in the HRA collection to HGNC and Ensembl IDs.

Reads every ccf:ExperimentUsedAntibody from hra.ttl (one per antibody used in
an OMAP experiment) and resolves its ccf:detects targets, which are HGNC IRIs
(http://identifiers.org/hgnc/N), against the HGNC complete set. Antibodies that
detect several genes (e.g. pan-cytokeratin) list them pipe-separated; antibodies
with no detects target (mostly secondary antibodies) are kept with empty gene
columns. Antibody details come from the antibody's RRID class, which may be
shared by several OMAPs, so differing values are joined with "; ".

needs_review is set when (a) the antibody name mentions a gene (approved symbol
or CD name, outside brackets) that is not among its detects targets, (b) the
antibody is in SUSPECTED_TARGETS, a hand-curated list of targets that look wrong
from the antibody's protein name, or (c) a primary antibody has no target.

Usage: python3 build_omap_hgnc_map.py [-o OUTPUT.csv]
       [--hra downloads/hra.ttl] [--hgnc downloads/hgnc_complete_set.txt]
"""
import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from pyoxigraph import RdfFormat, Store

from hgnc_lookup import Hgnc

HGNC_IRI = "http://identifiers.org/hgnc/"
RRID_IRI = "https://identifiers.org/RRID:"
PREFIXES = """
PREFIX ccf: <http://purl.org/ccf/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX hgnc: <http://purl.bioontology.org/ontology/HGNC/>
"""

USED_QUERY = PREFIXES + """
SELECT ?used ?experiment ?cycle ?record ?core ?rrid WHERE {
  ?used a ccf:ExperimentUsedAntibody ;
        ccf:used_by_experiment ?experiment ;
        ccf:based_on ?instance .
  ?instance a ?rrid .
  ?rrid rdfs:subClassOf ccf:Antibody .
  OPTIONAL { ?used ccf:cycle_number ?cycle }
  OPTIONAL { ?used ccf:record_number ?record }
  OPTIONAL { ?used ccf:is_core_panel ?core }
}
"""

DETECTS_QUERY = PREFIXES + """
SELECT ?used ?target ?symbol WHERE {
  ?used a ccf:ExperimentUsedAntibody ; ccf:detects ?target .
  OPTIONAL { ?target hgnc:symbol ?symbol }
}
"""

EXPERIMENT_QUERY = PREFIXES + """
SELECT ?experiment ?method ?organ ?organ_label ?preservation WHERE {
  ?experiment a ccf:MultiplexedAntibodyBasedImagingExperiment .
  OPTIONAL { ?experiment ccf:study_method ?method }
  OPTIONAL { ?experiment ccf:sample_organ ?organ . OPTIONAL { ?organ rdfs:label ?organ_label } }
  OPTIONAL { ?experiment ccf:tissue_preservation ?preservation }
}
"""

ANTIBODY_FIELDS = {
    "antibody_label": "rdfs:label",
    "antibody_target": "rdfs:comment",
    "clone_id": "ccf:clone_id",
    "catalog_number": "ccf:catalog_number",
    "vendor_name": "ccf:antibody_manufacturer",
    "host": "ccf:host",
    "clonality": "ccf:clonality",
    "conjugate": "ccf:conjugate",
}
ANTIBODY_QUERY = PREFIXES + """
SELECT ?rrid ?field ?value WHERE {
  VALUES (?field ?p) { %s }
  ?rrid rdfs:subClassOf ccf:Antibody ; ?p ?value .
}
""" % " ".join(f'("{k}" {v})' for k, v in ANTIBODY_FIELDS.items())

# (RRID, detected symbol) -> (suggested symbol or None, why the target looks wrong). Curated by reviewing antibody
# names that the automatic check cannot catch (protein names, not gene symbols).
SUSPECTED_TARGETS = {
    ("AB_2562826", "CD74"): ("HLA-DRA", "HLA-DR antibody; HLA-DR is HLA-DRA/HLA-DRB*, CD74 is the invariant chain"),
    ("AB_893582", "CD74"): ("HLA-DRA", "HLA-DR antibody; HLA-DR is HLA-DRA/HLA-DRB*, CD74 is the invariant chain"),
    ("AB_10585132", "VIM"): ("CLU", "Apolipoprotein J is clusterin (CLU), not vimentin"),
    ("AB_2212386", "TNS1"): ("UMOD", "uromodulin is UMOD (OMAP-14 targets of TNS1/UMOD/SYNPO/PODXL look shuffled)"),
    ("AB_1130061", "TNC"): ("TRD", "TCR delta is the TRD locus, not tenascin C"),
    ("AB_2243952", "SNAP25"): ("CASP3", "active caspase-3 is CASP3, not SNAP25"),
    ("AB_2890241", "ATP1B3"): ("ATP1A1", "clone EP1845Y is a Na/K-ATPase alpha-1 antibody (other OMAPs map it to ATP1A1)"),
    ("AB_3477119", "LGALS1"): ("LGALS3", "Mac-2/galectin-3 is LGALS3, not LGALS1"),
    ("AB_3676530", "H3-4"): (None, "histone H3 is encoded by many genes; H3-4 alone may be too narrow"),
    ("AB_778073", "FUT4"): (None, "SSEA-4 is a glycolipid antigen; FUT4 makes SSEA-1/Lewis X, not SSEA-4"),
    ("AB_2857381", "CKAP4"): (None, "verify the antigen of this 'Plasma Cell' antibody"),
    ("AB_2631163", "CD3D"): ("CD3E", "generic CD3 antibody; most CD3 antibodies bind CD3E"),
}
SECONDARY_RE = re.compile(r"\bIg[GMYAE]\b|secondary|labeling kit", re.IGNORECASE)

OUTPUT_COLUMNS = [
    "omap", "omap_version", "organ", "organ_uberon_id", "method", "tissue_preservation",
    "cycle_number", "record_number", "is_core_panel", "rrid", "antibody_label",
    "detects_symbol", "ensembl", "hgnc_id", "hgnc_symbol", "hgnc_match", "uniprot",
    "gene_description", "clone_id", "catalog_number", "vendor_name", "host", "clonality",
    "conjugate", "used_antibody", "needs_review", "suggested_hgnc_symbol", "note",
]


def val(term):
    return "" if term is None else term.value


def parse_experiment_iri(iri):
    # https://purl.humanatlas.io/omap/<name>/<version>#OMAP-N
    base, _, frag = iri.partition("#")
    name, version = base.rstrip("/").split("/")[-2:]
    return frag or name, name, version


def curie(iri):
    for prefix, short in (("http://purl.obolibrary.org/obo/", ""),):
        if iri.startswith(prefix):
            return iri[len(prefix):].replace("_", ":", 1)
    return iri


def load(store):
    experiments = {val(s["experiment"]): s for s in store.query(EXPERIMENT_QUERY)}
    detects = defaultdict(list)
    for s in store.query(DETECTS_QUERY):
        detects[val(s["used"])].append((val(s["target"]), val(s["symbol"])))
    antibodies = defaultdict(lambda: defaultdict(set))
    for s in store.query(ANTIBODY_QUERY):
        antibodies[val(s["rrid"])][val(s["field"])].add(val(s["value"]).strip())
    return experiments, detects, antibodies


def build_rows(store, hgnc):
    experiments, detects, antibodies = load(store)
    for s in store.query(USED_QUERY):
        used, experiment, rrid = val(s["used"]), val(s["experiment"]), val(s["rrid"])
        notes = []
        exp = experiments.get(experiment)
        if exp is None and experiment.replace("#OMAP_", "#OMAP-") in experiments:
            notes.append(f"used_by_experiment {experiment} is undefined; "
                         f"using {experiment.replace('#OMAP_', '#OMAP-')}")
            experiment = experiment.replace("#OMAP_", "#OMAP-")
            exp = experiments[experiment]
        omap, name, version = parse_experiment_iri(experiment)
        ab = antibodies.get(rrid, {})
        joined = {k: "; ".join(sorted(v for v in ab.get(k, ()) if v)) for k in ANTIBODY_FIELDS}

        targets = sorted(detects.get(used, []), key=lambda t: int(t[0].rsplit("/", 1)[-1])
                         if t[0].startswith(HGNC_IRI) else 0)
        hits = []
        for iri, symbol in targets:
            if not iri.startswith(HGNC_IRI):
                notes.append(f"detects non-HGNC IRI {iri}")
                continue
            hgnc_id = "HGNC:" + iri[len(HGNC_IRI):]
            h = hgnc.by_id.get(hgnc_id)
            if h is None:
                notes.append(f"{hgnc_id} not an approved HGNC ID")
                continue
            hits.append(h)
            if symbol and symbol != h["symbol"]:
                notes.append(f"{hgnc_id} symbol in hra.ttl is {symbol}, HGNC now {h['symbol']}")
        review = any(n.startswith(("HGNC", "detects")) or "not an approved" in n for n in notes)
        if not targets:
            desc = joined["antibody_target"] or joined["antibody_label"]
            notes.append("no detects target" + (f" ({desc})" if desc else ""))
            if not SECONDARY_RE.search(joined["antibody_label"] + " " + joined["antibody_target"]):
                notes.append("primary antibody with no target gene in the OMAP")
                review = True
        suggested = []
        hit_ids = {h["hgnc_id"] for h in hits}
        named = hgnc.genes_named_in(joined["antibody_label"])
        if hits and named and not hit_ids & set(named.values()):
            others = ", ".join(f"{tok} = {hgnc.by_id[i]['symbol']}" for tok, i in named.items())
            notes.append(f"antibody name mentions {others}, but the OMAP says it detects "
                         f"{'|'.join(h['symbol'] for h in hits)}")
            suggested += [hgnc.by_id[i]["symbol"] for i in dict.fromkeys(named.values())]
            review = True
        for h in hits:
            suspect = SUSPECTED_TARGETS.get((rrid[len(RRID_IRI):], h["symbol"]))
            if suspect:
                notes.append(f"suspected wrong target {h['symbol']}: {suspect[1]}")
                suggested += [suspect[0]] if suspect[0] else []
                review = True
        if any(not h["ensembl_gene_id"] for h in hits):
            notes.append("HGNC record has no Ensembl gene ID")
        if hits:
            match = "hgnc_id_multiple" if len(hits) > 1 else "hgnc_id"
        else:
            match = ""

        yield {
            "omap": omap,
            "omap_version": version,
            "organ": val(exp["organ_label"]) if exp else "",
            "organ_uberon_id": curie(val(exp["organ"])) if exp else "",
            "method": val(exp["method"]) if exp else "",
            "tissue_preservation": val(exp["preservation"]) if exp else "",
            "cycle_number": val(s["cycle"]),
            "record_number": val(s["record"]),
            "is_core_panel": val(s["core"]),
            "rrid": rrid[len(RRID_IRI):] if rrid.startswith(RRID_IRI) else rrid,
            "antibody_label": joined["antibody_label"],
            "detects_symbol": "|".join(sym for _, sym in targets),
            "ensembl": "|".join(h["ensembl_gene_id"] for h in hits),
            "hgnc_id": "|".join(h["hgnc_id"] for h in hits),
            "hgnc_symbol": "|".join(h["symbol"] for h in hits),
            "hgnc_match": match,
            "uniprot": "|".join(h["uniprot_ids"] for h in hits),
            "gene_description": "|".join(h["name"] for h in hits),
            **{k: joined[k] for k in ("clone_id", "catalog_number", "vendor_name", "host",
                                      "clonality", "conjugate")},
            "used_antibody": used,
            "needs_review": str(review).lower(),
            "suggested_hgnc_symbol": "|".join(dict.fromkeys(suggested)),
            "note": "; ".join(notes),
        }


def sort_key(r):
    num = r["omap"].rsplit("-", 1)[-1]
    return (int(num) if num.isdigit() else 0, r["omap_version"],
            int(r["record_number"] or 0), r["used_antibody"])


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-o", "--output", type=Path, default=Path("output/omap-hgnc-map.csv"))
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

    matches = Counter(r["hgnc_match"] for r in rows)
    genes = len({g for r in rows for g in r["hgnc_id"].split("|") if g})
    omaps = len({r["omap"] for r in rows})
    review = sum(r["needs_review"] == "true" for r in rows)
    print(f"Wrote {len(rows)} OMAP antibodies ({omaps} OMAPs, {genes} distinct HGNC genes) to "
          f"{args.output}; matches: {dict(matches)}; {review} need review", file=sys.stderr)


if __name__ == "__main__":
    main()
