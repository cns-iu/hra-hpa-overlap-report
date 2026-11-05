# HRA × HPA overlap report

This project maps every gene, protein and antibody target referenced by the
[Human Reference Atlas (HRA)](https://humanatlas.io) and related resources to
standard gene identifiers: **HGNC ID**, **approved HGNC symbol**, **Ensembl
gene ID** and **UniProt**. It does the same for the
[Human Protein Atlas (HPA)](https://www.proteinatlas.org).

With every source keyed to the same IDs, you can compare them directly, for
example to measure how much of the HRA's biomarker space HPA covers. Each map
also flags records where the source data looks wrong, so subject matter experts
can review them.

## Results

### 📊 [HPA × HRA overlap report](overlap-report/)

**Start here.** The report opens with an executive summary, then covers:
- how much of each HRA resource the HPA covers
- which HRA genes are missing from the HPA, and why
- how much of the HPA the HRA uses
- which resources share which genes
- HPA antibody evidence for the HRA's markers
- how the data-quality fixes affect the numbers

In short, the HPA covers about 94% of the genes the HRA references. Nearly
all the missing genes are non-coding RNAs, pseudogenes or T-cell receptor
loci. The HRA, in turn, uses about 9% of the HPA. Gene-level detail is in
[`overlap-report/hpa-hra-overlap.csv`](overlap-report/hpa-hra-overlap.csv).

### Gene maps

All gene maps are CSV files in [`output/`](output/). Each source gets one
`*-hgnc-map.csv` file.

| File | Source | One row per | Rows | Distinct HGNC genes | Rows flagged for review |
|---|---|---|---:|---:|---:|
| [`hpa-hgnc-map.csv`](output/hpa-hgnc-map.csv) | Human Protein Atlas gene table | HPA gene | 20,162 | 19,797 | 10 |
| [`asctb-hgnc-map.csv`](output/asctb-hgnc-map.csv) | HRA ASCT+B tables (gene, protein and other biomarkers) | ASCT+B biomarker record | 27,929 | 1,871 | 1,546 |
| [`omap-hgnc-map.csv`](output/omap-hgnc-map.csv) | HRA Organ Mapping Antibody Panels (OMAPs) | antibody used in an OMAP | 1,200 | 374 | 32 |
| [`avr-hgnc-map.csv`](output/avr-hgnc-map.csv) | HRA Antibody Validation Reports (AVRs) | AVR antibody | 1,053 | 319 | 44 |
| [`ctann-hgnc-map.csv`](output/ctann-hgnc-map.csv) | HRA cell type annotation (CTann) crosswalks, genes inferred via ASCT+B | (CTann mapping, gene) pair | 40,107 | 1,157 | 1,729 |
| [`deepcell-types-hgnc-map.csv`](output/deepcell-types-hgnc-map.csv) | [DeepCell Types](https://github.com/vanvalenlab/deepcell-types) channel names | imaging channel | 71 | 77 | 14 |
| [`review-items.csv`](output/review-items.csv) | All of the above | distinct issue for expert review | 291 | – | – |

Counts are from the 2026-11-05 build. They change when the sources are updated.

### Columns shared by every map

| Column | Meaning |
|---|---|
| `ensembl` | Ensembl gene ID(s) |
| `hgnc_id` | HGNC ID(s), e.g. `HGNC:1678` |
| `hgnc_symbol` | Current approved HGNC symbol(s) |
| `hgnc_match` | How the gene was identified (see below) |
| `uniprot` | UniProt accession(s) |
| `gene_description` | HGNC gene name |
| `needs_review` | `true` if a person should check the row |
| `note` | Why the row was flagged, plus any other caveats |

A field holds several values, separated by `|`, when one record maps to
several genes. Examples are pan-cytokeratin antibodies, CD3 (CD3D, CD3E and
CD3G) and ambiguous aliases. The other columns describe the source record,
such as ASCT+B table and record number, OMAP and cycle, antibody RRID, clone
or CTann label.

**`hgnc_match` values**

| Value | Meaning |
|---|---|
| `hgnc_id` | The source gave an HGNC ID that resolves to an approved gene |
| `ensembl` | Matched by Ensembl gene ID |
| `symbol` | The source label is an approved HGNC symbol |
| `prev_symbol` | The source label is a previous HGNC symbol |
| `alias_symbol` | The source label is an HGNC alias. Less reliable: check it |
| `label_symbol` | ASCT+B only. The record's HGNC ID pointed to a different gene than its label, which is an exact approved symbol, so the label's gene was used |
| `target_symbol` | AVR only. Same as `label_symbol`; the AVR PDF filename was used to decide |
| `*_multiple` | The match is ambiguous; all candidate genes are listed |
| *(empty)* | No gene could be identified |

### Expert review

[`output/review-items.csv`](output/review-items.csv) collects every row with
`needs_review = true`. Identical issues are grouped, with a count of how many
rows each affects, so each problem appears only once. Each issue has a
`category`:

| Category | What it means |
|---|---|
| `wrong-id-likely` | The source ID and the source's own label or antibody name point to different genes, or the ID looks truncated. For example, ASCT+B labels "SOX10" with HGNC:1119 (BST2) instead of HGNC:11190 |
| `withdrawn-id` | The source uses an HGNC ID that HGNC has since withdrawn |
| `uniprot-mismatch` | The source's UniProt accession belongs to a different gene |
| `ambiguous-alias` | The label is an alias shared by several genes |
| `label-check` | The label is not a known name of the gene. Often harmless, such as a protein name; sometimes a wrong ID |
| `missing-target` | A primary antibody has no target gene |
| `non-gene-type` | A lipid or metabolite biomarker resolves to a gene |
| `ensembl-mismatch` | HPA's Ensembl ID is not in HGNC, so the gene was matched by symbol |
| `curated-mapping` | A hand-made marker-to-gene mapping that depends on the antibody used (DeepCell) |

`suggested_hgnc_symbol` holds the likely correct gene, where one can be
inferred. It appears in the ASCT+B, OMAP and AVR maps and in the review file.
The overlap report's "corrected" numbers use these suggestions.

Where the source can be fixed, please fix it there. The likely wrong IDs
should go back to the ASCT+B, OMAP and AVR curators.

## How each source is mapped

- **HPA:** each HPA gene is joined to HGNC by Ensembl gene ID. If HGNC doesn't
  know that Ensembl ID, the HPA gene name is tried as an approved symbol, and
  the row is flagged.
- **ASCT+B:** reads every biomarker record in `hra.ttl`, including genes,
  proteins, proteoforms, lipids and metabolites.
  - Records that point to an HGNC ID are checked against their label.
  - Records with only a placeholder ID (`ASCTB-TEMP_*`) are matched by label:
    approved symbol, then Ensembl ID, then previous symbol, then alias.
  - Placeholder records that match no gene are left out. These are mostly
    lipids, metabolites and protein complexes.
- **OMAP:** reads every antibody used in an OMAP experiment and maps the
  HGNC IDs it detects. Each antibody name is checked for a gene it names but
  doesn't detect, such as "GFAP antibody" detecting FAP. A short hand-curated
  list in `SUSPECTED_TARGETS` covers errors that are only visible from
  protein names.
- **AVR:** reads `avr.ttl`. Each record's HGNC ID, target symbol and UniProt
  accession are cross-checked. When the ID and the symbol disagree, the AVR
  PDF filename decides which one is right. AVR records have no antibody
  names, and some PDFs repeat the wrong gene (e.g. `CDC123_OMAP17.pdf` for a
  CD123 antibody). So AVR rows for an antibody that the OMAP map flags inherit
  that flag and its suggested gene.
- **CTann:** CTann crosswalks link an annotation tool's cell type labels to
  Cell Ontology (CL) terms. They contain no genes, so genes are **inferred**:
  each CL term is joined to the ASCT+B records whose most specific cell type
  is that term, and those records' gene markers are collected.
  - `same_organ` shows whether any supporting ASCT+B record covers the
    CTann's organ.
  - `asctb_record_count` and `asctb_tables` show the evidence.
- **DeepCell Types:** channel names are protein/marker names such as `CD20`
  or `PanCK`. They are mapped to genes through a hand-curated table in
  [`src/build_deepcell_types_hgnc_map.py`](src/build_deepcell_types_hgnc_map.py).

## Setup

You need Python 3.10 or later, `curl`, and about 300 MB of disk space for the
downloads.

```bash
git clone https://github.com/cns-iu/hra-hpa-overlap-report.git
cd hra-hpa-overlap-report
pip install -r requirements.txt
```

[pyoxigraph](https://pyoxigraph.readthedocs.io) loads the 230 MB HRA
collection graph in a few seconds. Each map builds in under 10 seconds.

## Running

Run the scripts from the repository root:

```bash
bash scripts/10-download-resources.sh   # fetch source data into downloads/
bash scripts/20-generate-maps.sh        # build every map and the review list into output/
bash scripts/50-overlap-report.sh       # build the HPA × HRA overlap report from those maps
```

| Step | Inputs (in `downloads/`) | Where they come from |
|---|---|---|
| HRA collection graph | `hra.ttl` | https://purl.humanatlas.io/collection/hra (ASCT+B, OMAP and CTann) |
| AVR graph | `avr.ttl` | https://purl.humanatlas.io/graph/avr |
| HPA gene table | `proteinatlas.tsv.zip` | https://www.proteinatlas.org/download/proteinatlas.tsv.zip |
| HGNC | `hgnc_complete_set.txt`, `hgnc_withdrawn.txt` | https://www.genenames.org (HGNC public downloads) |
| DeepCell channels | `deepcell-types-channels.yaml` | the vanvalenlab/deepcell-types repository |

All sources are fetched at their latest version. Re-running the download
script updates the maps to current HRA, HPA and HGNC releases.

You can also build a single map. Every script takes `-o` for the output path
and `--help` for its other options:

```bash
python3 src/build_omap_hgnc_map.py -o output/omap-hgnc-map.csv
```

## Repository layout

```
scripts/
  10-download-resources.sh    download all source data
  20-generate-maps.sh         build all maps and the review list
  50-overlap-report.sh        build the overlap report, figures and gene table
src/
  hgnc_lookup.py              shared HGNC matching (IDs, symbols, aliases, withdrawn IDs, label checks)
  build_hpa_hgnc_map.py
  build_asctb_hgnc_map.py
  build_omap_hgnc_map.py
  build_avr_hgnc_map.py
  build_ctann_hgnc_map.py     uses build_asctb_hgnc_map.py for its gene matching
  build_deepcell_types_hgnc_map.py
  build_review_items.py       gathers flagged rows into output/review-items.csv
  build_overlap_report.py     HPA × HRA overlap report, figures and gene table
output/                       generated gene maps and review-items.csv (committed)
overlap-report/               README.md (the report), figures/ and hpa-hra-overlap.csv (committed)
downloads/                    source data (not committed)
```

## Caveats

- **The source data is not changed, with three exceptions.**
  - In ASCT+B, when an HGNC ID contradicts its own label and the label is an
    exact approved symbol, the label's gene is used (`label_symbol`).
  - In AVR, when the ID and the symbol disagree, the PDF filename decides
    (`target_symbol`).
  - In OMAP, one experiment link that is broken in `hra.ttl` (`#OMAP_37`
    instead of `#OMAP-37`) is repaired.
  - The overlap report's "corrected" numbers apply every
    `suggested_hgnc_symbol`. The report also shows the "given" and
    "confident" numbers.

  In every case the original value is kept and the row is flagged. All other
  suspected errors are reported in `note` and `suggested_hgnc_symbol` only.
- **CTann genes are inferences**, not statements made by the CTann
  crosswalks. Cell types that no ASCT+B record uses as its most specific type
  get no genes; 959 of 3,466 mappings fall in this group.
- **Some matching rules are curated by hand** and should be confirmed by
  experts:
  - the DeepCell marker table
  - the OMAP `SUSPECTED_TARGETS` list
  - the `PROTEIN_NAMES` table in
    [`src/hgnc_lookup.py`](src/hgnc_lookup.py), which treats common protein
    names such as "aSMA" or "E-cadherin" as valid names for their genes
- **Alias matches are the weakest link.** HGNC aliases are not unique, so
  `alias_symbol*` matches are always flagged.

## License

[MIT](LICENSE) © Cyberinfrastructure for Network Science Center at Indiana University.
Source data remain under their providers' terms (HRA: CC BY 4.0; HPA, HGNC and DeepCell
Types: see their sites).
