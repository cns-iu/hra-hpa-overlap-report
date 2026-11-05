#!/usr/bin/env python3
"""Map DeepCell Types channel names to HGNC gene symbols, HGNC IDs and Ensembl IDs.

Channel names are protein/marker names (e.g. CD20, PanCK), not gene symbols, so
each canonical marker is mapped to its gene(s) via a curated table below. IDs are
then looked up by approved symbol in the HGNC complete set. Markers encoded by
several genes list all of them pipe-separated, like hpa-hgnc-map.csv does.

The marker -> gene table was curated by hand (not taken from an external
source). Markers in NEEDS_REVIEW depend on which antibody/clone the dataset
used, so a subject matter expert should confirm them.

Usage: python3 build_deepcell_types_hgnc_map.py [-o OUTPUT.csv]
       [--channels downloads/deepcell-types-channels.yaml]
       [--hgnc downloads/hgnc_complete_set.txt]
"""
import argparse
import csv
import sys
from pathlib import Path

import yaml

# canonical marker -> (gene symbols, note)
MARKER_GENES = {
    "Beta-catenin": (["CTNNB1"], ""),
    "CD11b": (["ITGAM"], ""),
    "CD11c": (["ITGAX"], ""),
    "CD138": (["SDC1"], ""),
    "CD16": (["FCGR3A", "FCGR3B"], "CD16a and CD16b; most antibodies detect both"),
    "CD163": (["CD163"], ""),
    "CD20": (["MS4A1"], ""),
    "CD209": (["CD209"], ""),
    "CD3": (["CD3E", "CD3D", "CD3G"], "CD3 complex; most antibodies target CD3E"),
    "CD31": (["PECAM1"], ""),
    "CD4": (["CD4"], ""),
    "CD45": (["PTPRC"], ""),
    "CD45RO": (["PTPRC"], "CD45RO isoform of PTPRC"),
    "CD56": (["NCAM1"], ""),
    "CD63": (["CD63"], ""),
    "CD68": (["CD68"], ""),
    "CD8": (["CD8A", "CD8B"], "most CD8 antibodies target CD8A only"),
    "EGFR": (["EGFR"], ""),
    "FoxP3": (["FOXP3"], ""),
    "HLA-Class-2": (["HLA-DRA", "HLA-DRB1"],
                    "MHC class II; HLA-DR only (HLA-DP/DQ and other DRB genes not listed)"),
    "HLA-Class-1": (["HLA-A", "HLA-B", "HLA-C"],
                    "MHC class I heavy chains; some pan-class-I antibodies bind B2M instead"),
    "IDO": (["IDO1"], ""),
    "CK17": (["KRT17"], ""),
    "CK6": (["KRT6A", "KRT6B", "KRT6C"], "antibodies typically cross-react with CK6 isoforms"),
    "CK7": (["KRT7"], ""),
    "Ki67": (["MKI67"], ""),
    "LAG3": (["LAG3"], ""),
    "MPO": (["MPO"], ""),
    "OX40": (["TNFRSF4"], ""),
    "PDL1": (["CD274"], ""),
    "PD1": (["PDCD1"], ""),
    "PanCK": ([], "pan-cytokeratin antibody cocktail; gene set depends on clone (e.g. AE1/AE3 "
                  "detects KRT1-8,10,14-16,19 per OMAPs), so left unmapped"),
    "SMA": (["ACTA2"], ""),
    "Vimentin": (["VIM"], ""),
    "p53": (["TP53"], ""),
    "DEFA5": (["DEFA5"], ""),
    "Bcl-2": (["BCL2"], ""),
    "c-kit": (["KIT"], ""),
    "CD123": (["IL3RA"], ""),
    "CD127": (["IL7R"], ""),
    "CD15": (["FUT4"], "CD15 is the Lewis X carbohydrate epitope; FUT4 synthesizes it"),
    "CD161": (["KLRB1"], ""),
    "CD19": (["CD19"], ""),
    "CD206": (["MRC1"], ""),
    "CD21": (["CR2"], ""),
    "CD25": (["IL2RA"], ""),
    "CD34": (["CD34"], ""),
    "CD36": (["CD36"], ""),
    "CD38": (["CD38"], ""),
    "CD44": (["CD44"], ""),
    "CD57": (["B3GAT1"], "CD57 is the HNK-1 carbohydrate epitope; B3GAT1 synthesizes it"),
    "CD66": (["CEACAM1", "CEACAM3", "CEACAM5", "CEACAM6", "CEACAM8"],
             "CD66a-e (CEACAM family); antibody specificity unknown"),
    "CD69": (["CD69"], ""),
    "CD7": (["CD7"], ""),
    "CD90": (["THY1"], ""),
    "Col4": (["COL4A1", "COL4A2"], "collagen IV; main alpha1/alpha2 chains listed"),
    "FAP": (["FAP"], ""),
    "Intelectin-1": (["ITLN1"], ""),
    "MUC1": (["MUC1"], ""),
    "MUC2": (["MUC2"], ""),
    "NKG2D": (["KLRK1"], ""),
    "OLFM4": (["OLFM4"], ""),
    "PDPN": (["PDPN"], ""),
    "SOX9": (["SOX9"], ""),
    "Synaptophysin": (["SYP"], ""),
}

# Markers whose gene list depends on antibody specificity or biology choices.
NEEDS_REVIEW = {"CD16", "CD3", "CD8", "HLA-Class-2", "HLA-Class-1", "CK6", "PanCK", "CD15",
                "CD57", "CD66", "Col4"}

# Raw channel names whose mapping to a marker looks like a typo upstream.
CHANNEL_NOTES = {"NKG2G": "channel name NKG2G looks like a typo for NKG2D (as the YAML maps it); "
                          "NKG2G is not a known gene/protein"}

OUTPUT_COLUMNS = [
    "channel", "marker", "ensembl", "gene_symbol", "hgnc_id", "hgnc_symbol",
    "hgnc_match", "uniprot", "gene_description", "needs_review", "note",
]


def read_hgnc_by_symbol(path):
    with open(path, encoding="utf-8", newline="") as f:
        return {r["symbol"]: r for r in csv.DictReader(f, delimiter="\t")
                if r["status"] == "Approved"}


def build_rows(channels, hgnc):
    for channel, marker in channels.items():
        if marker not in MARKER_GENES:
            raise SystemExit(f"No curated gene mapping for marker {marker!r} (channel {channel!r})")
        symbols, note = MARKER_GENES[marker]
        missing = [s for s in symbols if s not in hgnc]
        if missing:
            raise SystemExit(f"Not an approved HGNC symbol: {missing} (marker {marker!r})")
        hits = [hgnc[s] for s in symbols]
        yield {
            "channel": channel,
            "marker": marker,
            "ensembl": "|".join(h["ensembl_gene_id"] for h in hits),
            "gene_symbol": "|".join(symbols),
            "hgnc_id": "|".join(h["hgnc_id"] for h in hits),
            "hgnc_symbol": "|".join(h["symbol"] for h in hits),
            "hgnc_match": "symbol_multiple" if len(hits) > 1 else "symbol" if hits else "",
            "uniprot": "|".join(h["uniprot_ids"] for h in hits),
            "gene_description": "|".join(h["name"] for h in hits),
            "needs_review": str(marker in NEEDS_REVIEW or channel in CHANNEL_NOTES).lower(),
            "note": "; ".join(n for n in (CHANNEL_NOTES.get(channel, ""), note) if n),
        }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-o", "--output", type=Path, default=Path("output/deepcell-types-hgnc-map.csv"))
    ap.add_argument("--channels", type=Path, default=Path("downloads/deepcell-types-channels.yaml"))
    ap.add_argument("--hgnc", type=Path, default=Path("downloads/hgnc_complete_set.txt"))
    args = ap.parse_args()

    channels = yaml.safe_load(args.channels.read_text(encoding="utf-8"))
    hgnc = read_hgnc_by_symbol(args.hgnc)
    rows = list(build_rows(channels, hgnc))
    with open(args.output, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, OUTPUT_COLUMNS)
        w.writeheader()
        w.writerows(rows)

    single = sum(r["hgnc_match"] == "symbol" for r in rows)
    multi = sum(r["hgnc_match"] == "symbol_multiple" for r in rows)
    print(f"Wrote {len(rows)} channels to {args.output}: {single} single-gene, "
          f"{multi} multi-gene, {len(rows) - single - multi} unmapped", file=sys.stderr)


if __name__ == "__main__":
    main()
