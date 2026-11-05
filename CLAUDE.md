# CLAUDE.md: notes for continuing this project

This file is for an AI assistant or developer picking up the work.
[README.md](README.md) is the public description; read it first for each
map's columns and the `hgnc_match` values. This file records what the README
doesn't: conventions, decisions, quirks in the source data, and how to extend
the code.

## What this project does

Repository: https://github.com/cns-iu/hra-hpa-overlap-report (default branch `main`).

It maps every gene, protein and antibody target in several Human Reference
Atlas (HRA) resources, the Human Protein Atlas (HPA) and DeepCell Types to
HGNC IDs, symbols, Ensembl IDs and UniProt accessions. It then reports how
the HPA overlaps the HRA. The owner is bherr (HuBMAP/HRA team, Indiana
University). Results are meant to be reviewed by subject matter experts
(SMEs) and browsed on GitHub.

## Pipeline

Always run from the repo root.

```
scripts/10-download-resources.sh   curl -fsSL everything into downloads/ (not committed)
scripts/20-generate-maps.sh        src/build_*_hgnc_map.py -> output/*-hgnc-map.csv, then review-items.csv
scripts/50-overlap-report.sh       src/build_overlap_report.py -> overlap-report/README.md, figures/, hpa-hra-overlap.csv
```

- **Step order inside `20-generate-maps.sh` matters.**
  - OMAP must run before AVR, because AVR inherits OMAP's flags through
    `--omap`.
  - CTann imports and calls `build_asctb_hgnc_map.build_rows`, so the ASCT+B
    matching logic lives in exactly one place.
  - `build_review_items.py` runs last.
- **Script names:** the 10/20/50 numbering is the owner's own choice. Keep new
  scripts in that style (e.g. `30-…`, `40-…`).
- **Environment:** Python 3.10. `pip install --user` is fine. Dependencies are
  pyoxigraph, rdflib (used only for the small `avr.ttl`), PyYAML and
  matplotlib.
- **Default paths:** scripts read from `downloads/` and write to `output/` by
  default.
- **Imports:** scripts run as `python3 src/x.py`, so
  `from hgnc_lookup import …` works without packaging.

## Code layout and conventions

- [`src/hgnc_lookup.py`](src/hgnc_lookup.py) holds the shared `Hgnc` class.
  Reuse it for any new source; don't write a new HGNC parser.
  - `by_id`, `by_symbol` (keys upper-cased), `by_ensembl`, and `by_prev` /
    `by_alias` (lists, because aliases aren't unique).
  - `match_label`: tries approved symbol, then Ensembl ID found in the
    label, then previous symbol, then alias. Returns the genes and the match
    type.
  - `is_known_name`: label comparison that ignores case, punctuation, a
    trailing +/- and `hi`/`lo` qualifiers. It also accepts the curated
    `PROTEIN_NAMES` table and genes found by `genes_named_in`.
  - `genes_named_in`: finds gene symbols or CD names in antibody names. It
    ignores text in brackets or parentheses (usually clone IDs) and the
    words in `NOT_GENE_WORDS` (KIT, APC, SP7, …).
  - `truncation_candidates`: given an ID and a label, looks for the ID plus
    one more digit whose symbol starts with the label (e.g. HGNC:1050 BID →
    HGNC:10500 S100B).
  - `describe_missing`: explains an unknown ID using
    `downloads/hgnc_withdrawn.txt`.
- **Every `*-hgnc-map.csv` shares these columns:** `ensembl`, `hgnc_id`,
  `hgnc_symbol`, `hgnc_match`, `uniprot`, `gene_description`,
  `needs_review` (`true` or `false`, never empty) and `note`. ASCT+B, OMAP
  and AVR also have `suggested_hgnc_symbol`. A field with several genes
  separates them with `|`. HPA's `ensembl` column is HPA's own key, so it can
  differ from HGNC's Ensembl ID on 10 symbol-fallback rows. That is
  intentional.
- **Source values are kept**, and suspected errors are flagged and explained
  in `note`. Only three corrections are applied automatically:
  - ASCT+B `label_symbol`: the label is exactly another gene's approved
    symbol.
  - AVR `target_symbol`: settled by the PDF filename.
  - OMAP_37: a broken experiment link is repaired.

  Every other fix goes in `suggested_hgnc_symbol`. The overlap report's
  "corrected" numbers apply those suggestions.
- **Adding a map:**
  1. Write `src/build_<source>_hgnc_map.py` following the existing scripts:
     argparse with `-o`, `--hgnc`, and a summary line on stderr.
  2. Add it to `scripts/20-generate-maps.sh`.
  3. Add it to `SOURCES` in `src/build_review_items.py`.
  4. Add it to `SOURCES` in `src/build_overlap_report.py`, choosing the kind:
     HRA, HRA-inferred or related.
  5. Update the README results table.
- **Checking after changes:** rebuild everything, then confirm that every
  `hgnc_id` resolves, that symbols and Ensembl IDs agree with HGNC, and that
  `needs_review` is always true or false. This was done by hand with a short
  inline Python check, so write it again if needed. Also diff gene sets
  against the previous output. Several fixes (the CTann blank-flag bug, the
  redundant SOX10 suggestion) were only caught this way.

## Working with the turtle files

- **Use pyoxigraph:** `Store().bulk_load(path=…, format=RdfFormat.TURTLE)`
  loads the 230 MB `hra.ttl` in about 5 seconds.
- **Avoid one big SPARQL query** with OPTIONAL or property-path joins. One
  such query took 9 minutes. Run several flat queries and join the results
  in Python, which takes about 8 seconds.
- **`hra.ttl`** is the HRA collection: 31 ASCT+B tables, 38 OMAPs and 10
  CTanns, plus the ontologies. Each table appears in only one version.
  - **ASCT+B:** `ccf:AsctbRecord` → `ccf:cell_type`, `ccf:anatomical_structure`,
    `ccf:gene_marker` and `ccf:protein_marker` (also `lipid_`, `metabolite_`
    and `proteoform_marker`) → `ccf:BiomarkerRecord`. Each biomarker has
    `ccf:source_concept` (an `http://identifiers.org/hgnc/N` IRI, or a
    `https://purl.org/ccf/ASCTB-TEMP_*` placeholder), `ccf:ccf_pref_label`
    and `ccf:ccf_biomarker_type`. `ccf:CellMarkerDescriptor` also exists but
    only covers HGNC IRIs, so it isn't used.
  - **OMAP:** `ccf:ExperimentUsedAntibody` → `ccf:detects` (always an HGNC
    IRI), `ccf:based_on` (an instance whose class is the RRID class, a
    subclass of `ccf:Antibody`, which holds manufacturer, clone and so on)
    and `ccf:used_by_experiment` (a
    `ccf:MultiplexedAntibodyBasedImagingExperiment` with organ and method).
    Several OMAPs share an RRID class, so its properties can have several
    values.
  - **CTann:** `ccf:CtAnnMapping` has `owl:annotatedSource`,
    `owl:annotatedTarget` (a `"CL:…"` string) and `sssom:*` labels. It has
    no genes. Genes are inferred by joining each CL term to the ASCT+B
    records whose *most specific* cell type (the highest CT `order_number`)
    is that term.
- **`avr.ttl`:** AVR antibodies are `</schema/avr#Antibody>` subjects with
  `hgnc_id`, `target_symbol`, `uniprot_accession_number`, `rrid`,
  `avr_pdf_filename` and so on. AVR rrids start with `RRID:`; OMAP rrids
  don't, so strip the prefix before matching.
- **Fetching examples:** single digital objects can be downloaded for
  inspection, e.g.
  `https://cdn.humanatlas.io/digital-objects/omap/33-ovary-codes/v1.0/graph.ttl`
  or `…/ctann/frmatch/v1.0/graph.ttl`.

## Known issues in the source data

These were found and flagged. Keep these checks working.

- **ASCT+B:**
  - **Truncated IDs:** SOX10 has HGNC:1119 (BST2) in 209 rows instead of
    HGNC:11190; S100 has HGNC:1050 (BID) in 277 rows instead of S100B,
    HGNC:10500.
  - **Wrong IDs:** CD2 has CD7's ID, and the bone-marrow table gives LYVE1,
    PROX1 and PDPN the IDs of STAU2, PRKY and PANK1.
  - **Withdrawn ID:** HGNC:1676 (CD3W) is used in 262 "CD3" rows.
  - **Ambiguous alias matches:** e.g. CGRP and MCT.
- **OMAP / AVR wrong targets:** CXCL13→RPL9, GFAP→FAP, FOXP3→CD79A,
  CD31→CALB1, CD123→CDC123, CD21→CR1, CD3e→MS4A1, Galectin-3→LGALS1,
  caspase-3→SNAP25, ApoJ→VIM and HLA-DR→CD74. In OMAP-14, TNS1, UMOD, SYNPO
  and PODXL are shuffled among each other. AVR PDFs are often named after
  the wrong gene, so the PDF filename can't be the only check. The curated
  list is `SUSPECTED_TARGETS` in `src/build_omap_hgnc_map.py`.
- **OMAP-37's used antibodies** point to an undefined `#OMAP_37`; the
  experiment is `#OMAP-37`.
- **HPA:** 357 entries have no HGNC record (novel or Ensembl-only genes),
  10 match only by symbol, and 284 use an older gene name.
- **Not HPA genes:** TRA, TRD and TRG are HGNC locus records with no
  Ensembl ID, so they can never match the HPA.

## Decisions the owner made or accepted

- Results go in `output/`, downloads in `downloads/`, and generated output is
  committed (the owner moved files there and added the scripts). The overlap
  report lives in the top-level `overlap-report/` folder as `README.md`, so
  GitHub shows it when the folder is opened. Its links point back to
  `../output/`.
- One row per source record, with multi-gene records separated by `|`. The
  exception is CTann, which has one row per (mapping, gene).
- DeepCell markers are mapped through a hand-curated table. PanCK is left
  unmapped.
- CTann genes are inferred through ASCT+B. This was the owner's suggestion.
- The overlap report leads with the "corrected" numbers and also shows
  "given" and "confident" (flagged rows dropped). The HRA union is
  ASCT+B ∪ OMAP ∪ AVR: CTann is excluded because it is a subset of ASCT+B,
  and DeepCell because it is not an HRA resource.
- The report is GitHub Markdown with Mermaid charts (flowchart, pie,
  `xychart-beta`) and matplotlib PNGs in light and dark versions, embedded
  with `<picture>`. Mermaid rendering on GitHub was never checked, and
  `xychart-beta horizontal` is the most likely to fail. The PNG colors come
  from the dataviz skill and passed its palette validator; re-run the
  validator if you change them.

## Possible next steps (not started)

- Check how the report renders on GitHub, especially the Mermaid xychart.
- Send `output/review-items.csv` to curators. Once the source data is
  fixed, the fixes in `SUSPECTED_TARGETS` and `PROTEIN_NAMES` may no longer
  be needed.
- Raise CTann coverage: 959 of 3,466 mappings get no genes because their CL
  term is never a *most specific* cell type in ASCT+B. Matching on less
  specific CT columns would cover more mappings but give less specific
  genes.
- Organ-level overlap: the HPA has tissue and cell-type expression columns
  in `proteinatlas.tsv.zip`, which could be compared with the organs in
  ASCT+B and OMAP.
- Turn the inline checks into a test script (`scripts/` or `src/`).
