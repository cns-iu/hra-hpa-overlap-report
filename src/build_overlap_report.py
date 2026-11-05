#!/usr/bin/env python3
"""Measure how the Human Protein Atlas (HPA) overlaps the HRA gene/protein maps.

Reads the *-hgnc-map.csv files in output/ plus the HPA gene table, compares the
sources by HGNC ID, and writes into overlap-report/ (--report-dir):

* README.md               Markdown report (GitHub-flavoured, Mermaid diagrams)
* figures/*.png           UpSet and HPA-reliability charts (light + dark)
* hpa-hra-overlap.csv     one row per gene with membership in every source

Every number is computed three ways:
  given      - genes exactly as each map resolved them
  corrected  - rows with a suggested_hgnc_symbol use the suggested gene instead
  confident  - rows flagged needs_review are dropped
The report leads with "corrected" and shows the others as a sensitivity check.

Usage: python3 build_overlap_report.py [--maps-dir output] [--report-dir overlap-report]
       [--hpa downloads/proteinatlas.tsv.zip]
"""
import argparse
import csv
import io
import itertools
import sys
import zipfile
from collections import Counter, OrderedDict, defaultdict
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

from hgnc_lookup import Hgnc  # noqa: E402

VARIANTS = ("given", "corrected", "confident")
MAIN = "corrected"

# key -> (label, map file, kind). CTann genes are inferred from ASCT+B, so they are
# reported but left out of the HRA union and the UpSet plot (they add nothing new).
SOURCES = OrderedDict([
    ("asctb", ("ASCT+B", "asctb-hgnc-map.csv", "HRA")),
    ("omap", ("OMAP", "omap-hgnc-map.csv", "HRA")),
    ("avr", ("AVR", "avr-hgnc-map.csv", "HRA")),
    ("ctann", ("CTann (inferred)", "ctann-hgnc-map.csv", "HRA-inferred")),
    ("deepcell", ("DeepCell Types", "deepcell-types-hgnc-map.csv", "related")),
])
HRA_KEYS = [k for k, v in SOURCES.items() if v[2] == "HRA"]
UPSET_KEYS = ["asctb", "omap", "avr", "deepcell"]

RELIABILITY = ["Enhanced", "Supported", "Approved", "Uncertain"]
NO_IHC = "No HPA tissue IHC"

# Validated with the dataviz skill's validate_palette.js (ordinal + categorical, both modes).
THEMES = {
    "light": {
        "surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781",
        "grid": "#e1e0d9", "axis": "#c3c2b7",
        "in_hpa": "#2a78d6", "not_hpa": "#eb6834", "none": "#e1e0d9",
        "ordinal": ["#104281", "#256abf", "#5598e7", "#86b6ef"],
        "ordinal_text": ["#ffffff", "#ffffff", "#0b0b0b", "#0b0b0b"],
    },
    "dark": {
        "surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7", "muted": "#898781",
        "grid": "#2c2c2a", "axis": "#383835",
        "in_hpa": "#3987e5", "not_hpa": "#d95926", "none": "#383835",
        "ordinal": ["#9ec5f4", "#5598e7", "#2a78d6", "#184f95"],
        "ordinal_text": ["#0b0b0b", "#0b0b0b", "#ffffff", "#ffffff"],
    },
}


# ---------------------------------------------------------------- data loading

def read_csv(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def ids_of(field):
    return [x for x in field.split("|") if x]


def row_genes(row, variant, hgnc):
    if variant == "confident" and row.get("needs_review") == "true":
        return []
    if variant == "corrected" and row.get("suggested_hgnc_symbol"):
        recs = [hgnc.by_symbol.get(s.upper()) for s in ids_of(row["suggested_hgnc_symbol"])]
        if all(recs):
            return [r["hgnc_id"] for r in recs]
    return ids_of(row["hgnc_id"])


def load_sources(maps_dir, hgnc):
    """-> {variant: {source: Counter(hgnc_id -> rows)}}, plus ASCT+B by biomarker type."""
    sets = {v: {} for v in VARIANTS}
    asctb_types = {v: defaultdict(set) for v in VARIANTS}
    for key, (_, fname, _) in SOURCES.items():
        rows = read_csv(maps_dir / fname)
        for v in VARIANTS:
            c = Counter()
            for r in rows:
                genes = row_genes(r, v, hgnc)
                c.update(genes)
                if key == "asctb":
                    asctb_types[v][r["biomarker_type"]].update(genes)
            sets[v][key] = c
    return sets, asctb_types


def load_hpa(maps_dir, hpa_zip):
    """-> {hgnc_id: HPA row} for HPA genes with an HGNC ID, and the HPA map rows."""
    hpa_map = read_csv(maps_dir / "hpa-hgnc-map.csv")
    hgnc_by_ens = {r["ensembl"]: r["hgnc_id"] for r in hpa_map if r["hgnc_id"]}
    hpa = {}
    with zipfile.ZipFile(hpa_zip) as zf:
        name = next(n for n in zf.namelist() if n.endswith(".tsv"))
        with zf.open(name) as raw:
            for r in csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8"), delimiter="\t"):
                for hid in ids_of(hgnc_by_ens.get(r["Ensembl"], "")):
                    hpa.setdefault(hid, r)
    return hpa, hpa_map


# ---------------------------------------------------------------- metrics

def pct(n, d):
    return 100.0 * n / d if d else 0.0


def fmt_pct(n, d):
    return f"{pct(n, d):.1f}%"


def reliability(hpa_row):
    if hpa_row is None:
        return None
    rel = hpa_row.get("Reliability (IH)", "").strip()
    return rel if rel in RELIABILITY else NO_IHC


def reliability_counts(genes, hpa):
    c = Counter(reliability(hpa[g]) for g in genes if g in hpa)
    return [c.get(k, 0) for k in RELIABILITY + [NO_IHC]]


def protein_classes(hpa_row):
    return [c.strip() for c in hpa_row.get("Protein class", "").split(",") if c.strip()]


# ---------------------------------------------------------------- figures

def style_axes(ax, t):
    ax.set_facecolor(t["surface"])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(t["axis"])
    ax.tick_params(colors=t["ink2"], labelsize=9, length=0)


def rounded_bar(ax, x, y, w, h, color, horizontal=False, radius=0.0):
    """Bar whose free end (top, or right if horizontal) is rounded; base stays square."""
    if w <= 0 or h <= 0:
        return
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={radius}",
                                linewidth=0, facecolor=color, mutation_aspect=None))
    # square off the base end
    if horizontal:
        ax.add_patch(Rectangle((x, y), min(w, w / 2), h, linewidth=0, facecolor=color))
    else:
        ax.add_patch(Rectangle((x, y), w, min(h, h / 2), linewidth=0, facecolor=color))


def draw_upset(path, combos, set_sizes, labels, t):
    """combos: list of (tuple of member keys, in_hpa, not_hpa) sorted by size."""
    n_sets, n_cols = len(labels), len(combos)
    fig = plt.figure(figsize=(max(7.5, 0.55 * n_cols + 3.2), 5.4), facecolor=t["surface"])
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1.6], hspace=0.05,
                          left=0.27, right=0.98, top=0.86, bottom=0.04)
    ax, axm = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    style_axes(ax, t)
    ymax = max(a + b for _, a, b in combos)
    bw = 0.62
    for i, (_, a, b) in enumerate(combos):
        x = i - bw / 2
        gap = ymax * 0.006
        if b and a:
            ax.add_patch(Rectangle((x, 0), bw, a, linewidth=0, facecolor=t["in_hpa"]))
            rounded_bar(ax, x, a + gap, bw, max(b - gap, 0), t["not_hpa"], radius=0.06)
        elif a:
            rounded_bar(ax, x, 0, bw, a, t["in_hpa"], radius=0.06)
        else:
            rounded_bar(ax, x, 0, bw, b, t["not_hpa"], radius=0.06)
        ax.text(i, a + b + ymax * 0.02, f"{a + b:,}", ha="center", va="bottom",
                fontsize=8.5, color=t["ink"])
    ax.set_xlim(-0.6, n_cols - 0.4)
    ax.set_ylim(0, ymax * 1.15)
    ax.set_xticks([])
    ax.yaxis.grid(True, color=t["grid"], linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["bottom"].set_visible(False)
    ax.set_ylabel("Genes in this exact combination", color=t["ink2"], fontsize=9)
    fig.text(0.27, 0.95, "Which resources share which genes", fontsize=12, weight="bold",
             color=t["ink"])
    # legend (always present for 2 series)
    for j, (name, col) in enumerate((("In HPA", t["in_hpa"]), ("Not in HPA", t["not_hpa"]))):
        fig.patches.append(Rectangle((0.27 + j * 0.11, 0.895), 0.014, 0.022, transform=fig.transFigure,
                                     facecolor=col, linewidth=0))
        fig.text(0.27 + j * 0.11 + 0.02, 0.906, name, fontsize=9, color=t["ink2"], va="center")

    # matrix
    axm.set_facecolor(t["surface"])
    for s in axm.spines.values():
        s.set_visible(False)
    axm.set_xlim(-0.6, n_cols - 0.4)
    axm.set_ylim(-0.5, n_sets - 0.5)
    axm.set_xticks([])
    axm.set_yticks(range(n_sets))
    axm.set_yticklabels([f"{labels[i]}  ({set_sizes[i]:,})" for i in range(n_sets)],
                        fontsize=9, color=t["ink"])
    axm.tick_params(length=0)
    axm.invert_yaxis()
    for r in range(n_sets):
        if r % 2 == 0:
            axm.axhspan(r - 0.5, r + 0.5, color=t["grid"], alpha=0.45, linewidth=0)
    for i, (members, _, _) in enumerate(combos):
        rows = [r for r in range(n_sets) if r in members]
        axm.scatter([i] * n_sets, range(n_sets), s=46, color=t["axis"], zorder=2, linewidths=0)
        if len(rows) > 1:
            axm.plot([i, i], [min(rows), max(rows)], color=t["ink"], linewidth=2, zorder=3)
        axm.scatter([i] * len(rows), rows, s=58, color=t["ink"], zorder=4,
                    edgecolors=t["surface"], linewidths=2)
    fig.savefig(path, dpi=180, facecolor=t["surface"])
    plt.close(fig)


def draw_reliability(path, rows, t):
    """rows: list of (label, counts list aligned with RELIABILITY + [NO_IHC])."""
    cats = RELIABILITY + [NO_IHC]
    colors = t["ordinal"] + [t["none"]]
    text_colors = t["ordinal_text"] + [t["ink2"]]
    fig, ax = plt.subplots(figsize=(9, 0.48 * len(rows) + 1.7), facecolor=t["surface"])
    fig.subplots_adjust(left=0.25, right=0.97, top=1 - 1.05 / (0.48 * len(rows) + 1.7), bottom=0.08)
    style_axes(ax, t)
    h = 0.62
    for i, (label, counts) in enumerate(rows):
        total = sum(counts)
        x = 0.0
        for j, c in enumerate(counts):
            w = pct(c, total)
            if w <= 0:
                continue
            ax.add_patch(Rectangle((x, i - h / 2), max(w - 0.25, 0.05), h, linewidth=0,
                                   facecolor=colors[j]))
            if w >= 7:
                ax.text(x + w / 2, i, f"{w:.0f}%", ha="center", va="center", fontsize=8,
                        color=text_colors[j])
            x += w
    ax.set_xlim(0, 100)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{lab}  (n={sum(c):,})" for lab, c in rows], fontsize=9, color=t["ink"])
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.xaxis.grid(True, color=t["grid"], linewidth=0.8)
    ax.set_axisbelow(True)
    fig.text(0.25, 1 - 0.32 / (0.48 * len(rows) + 1.7), "HPA tissue antibody reliability of genes in HPA",
             fontsize=12, weight="bold", color=t["ink"])
    lx = 0.25
    for j, cat in enumerate(cats):
        y = 1 - 0.72 / (0.48 * len(rows) + 1.7)
        fig.patches.append(Rectangle((lx, y - 0.012), 0.014, 0.024, transform=fig.transFigure,
                                     facecolor=colors[j], linewidth=0))
        fig.text(lx + 0.019, y, cat, fontsize=8.5, color=t["ink2"], va="center")
        lx += 0.019 + 0.0105 * len(cat) + 0.03
    fig.savefig(path, dpi=180, facecolor=t["surface"])
    plt.close(fig)


def picture(rel_base, alt):
    return (f'<picture>\n  <source media="(prefers-color-scheme: dark)" srcset="{rel_base}-dark.png">\n'
            f'  <img alt="{alt}" src="{rel_base}-light.png">\n</picture>')


# ---------------------------------------------------------------- report

def md_table(header, rows, align=None):
    align = align or ["---"] * len(header)
    out = ["| " + " | ".join(header) + " |", "| " + " | ".join(align) + " |"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def build(args):
    maps_dir = args.maps_dir
    hgnc = Hgnc(args.hgnc)
    sets, asctb_types = load_sources(maps_dir, hgnc)
    hpa, hpa_map = load_hpa(maps_dir, args.hpa)
    hpa_ids = set(hpa)
    hpa_without_hgnc = sum(1 for r in hpa_map if not r["hgnc_id"])
    label = {k: v[0] for k, v in SOURCES.items()}
    S = {k: set(c) for k, c in sets[MAIN].items()}
    hra = set().union(*(S[k] for k in HRA_KEYS))
    hra_in = hra & hpa_ids
    hra_out = hra - hpa_ids
    every = hra | S["deepcell"] | S["ctann"]

    def locus(g):
        return hgnc.by_id[g]["locus_group"] if g in hgnc.by_id else "unknown"

    # ---- per-source coverage (main variant)
    cov_rows = []
    for k in SOURCES:
        g = S[k]
        cov_rows.append((k, len(g), len(g & hpa_ids)))
    asctb_gene, asctb_prot = asctb_types[MAIN]["gene"], asctb_types[MAIN]["protein"]

    # ---- reliability
    rel_rows = [
        ("ASCT+B gene markers", reliability_counts(asctb_gene, hpa)),
        ("ASCT+B protein markers", reliability_counts(asctb_prot, hpa)),
        ("OMAP antibody targets", reliability_counts(S["omap"], hpa)),
        ("AVR antibody targets", reliability_counts(S["avr"], hpa)),
        ("CTann (inferred)", reliability_counts(S["ctann"], hpa)),
        ("DeepCell Types", reliability_counts(S["deepcell"], hpa)),
        ("All HPA genes", reliability_counts(hpa_ids, hpa)),
    ]

    # ---- UpSet combos (HRA direct sources + DeepCell)
    up_sets = [S[k] for k in UPSET_KEYS]
    combos = []
    for r in range(len(UPSET_KEYS), 0, -1):
        for members in itertools.combinations(range(len(UPSET_KEYS)), r):
            inside = set.intersection(*(up_sets[m] for m in members))
            outside = set().union(*(up_sets[m] for m in range(len(UPSET_KEYS)) if m not in members))
            exact = inside - outside
            if exact:
                combos.append((members, len(exact & hpa_ids), len(exact - hpa_ids), exact))
    combos.sort(key=lambda c: -(c[1] + c[2]))

    report_dir = args.report_dir
    figdir = report_dir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    for name, t in THEMES.items():
        draw_upset(figdir / f"upset-{name}.png", [(m, a, b) for m, a, b, _ in combos],
                   [len(s) for s in up_sets], [label[k] for k in UPSET_KEYS], t)
        draw_reliability(figdir / f"hpa-reliability-{name}.png", rel_rows, t)

    # ---- genes not in HPA
    out_locus = Counter(locus(g) for g in hra_out)
    coding_out = sorted(hgnc.by_id[g]["symbol"] for g in hra_out if locus(g) == "protein-coding gene")
    coding_out_no_ens = [s_ for s_ in coding_out if not hgnc.by_symbol[s_.upper()]["ensembl_gene_id"]]

    def coding_out_text():
        if not coding_out:
            return "None of them are protein-coding."
        txt = f"Only {len(coding_out)} are protein-coding ({', '.join(coding_out)})"
        if coding_out_no_ens == coding_out:
            txt += (". All of these are HGNC locus records with no single Ensembl gene (e.g. "
                    "T-cell receptor loci), so the HPA cannot list them as genes.")
        else:
            txt += (f"; {len(coding_out_no_ens)} of them are HGNC locus records with no single "
                    "Ensembl gene.")
        return txt
    asctb_rows_per_gene = sets[MAIN]["asctb"]
    top_missing = sorted(hra_out, key=lambda g: (-asctb_rows_per_gene.get(g, 0),
                                                hgnc.by_id.get(g, {}).get("symbol", g)))[:15]

    # ---- protein classes (HRA genes in HPA vs all HPA)
    cls_hra = Counter(c for g in hra_in for c in set(protein_classes(hpa[g])))
    cls_all = Counter(c for g in hpa_ids for c in set(protein_classes(hpa[g])))
    cls_rows = []
    for c, n in cls_hra.most_common(14):
        share_hra, share_all = pct(n, len(hra_in)), pct(cls_all[c], len(hpa_ids))
        cls_rows.append((c, f"{n:,}", f"{share_hra:.1f}%", f"{share_all:.1f}%",
                         f"{share_hra / share_all:.1f}×" if share_all else "–"))

    # ---- HPA evidence level
    ev_hra = Counter(hpa[g]["Evidence"] for g in hra_in)
    ev_all = Counter(hpa[g]["Evidence"] for g in hpa_ids)

    # ---- sensitivity
    sens = []
    for k in SOURCES:
        cells = [label[k]]
        for v in VARIANTS:
            g = set(sets[v][k])
            cells.append(f"{len(g & hpa_ids):,} / {len(g):,} ({fmt_pct(len(g & hpa_ids), len(g))})")
        sens.append(cells)
    hra_v = {v: set().union(*(set(sets[v][k]) for k in HRA_KEYS)) for v in VARIANTS}
    sens.append(["**HRA union**"] + [f"**{len(hra_v[v] & hpa_ids):,} / {len(hra_v[v]):,} "
                                      f"({fmt_pct(len(hra_v[v] & hpa_ids), len(hra_v[v]))})**"
                                      for v in VARIANTS])
    gained = (hra_v["corrected"] - hra_v["given"])
    lost = (hra_v["given"] - hra_v["corrected"])

    # ---- antibody resources specifics
    ab = S["omap"] | S["avr"]
    ab_in = ab & hpa_ids
    ab_strong = sum(1 for g in ab_in if reliability(hpa[g]) in ("Enhanced", "Supported"))
    hpa_strong = sum(1 for g in hpa_ids if reliability(hpa[g]) in ("Enhanced", "Supported"))
    ab_has_hpa_ab = sum(1 for g in ab_in if hpa[g].get("Antibody", "").strip())

    # ---- gene-level CSV
    with open(report_dir / "hpa-hra-overlap.csv", "w", encoding="utf-8", newline="") as f:
        cols = ["hgnc_id", "hgnc_symbol", "ensembl", "locus_group", "in_hpa", "hpa_evidence",
                "hpa_reliability_ih", "hpa_antibodies"] + \
               [f"{k}_rows" for k in SOURCES] + ["hra_sources", "n_hra_sources"]
        w = csv.writer(f)
        w.writerow(cols)
        for g in sorted(every | hpa_ids, key=lambda x: int(x.split(":")[1])):
            rec = hgnc.by_id.get(g, {})
            h = hpa.get(g)
            srcs = [label[k] for k in HRA_KEYS if g in S[k]]
            w.writerow([g, rec.get("symbol", ""), rec.get("ensembl_gene_id", ""), locus(g),
                        str(h is not None).lower(), h["Evidence"] if h else "",
                        (h.get("Reliability (IH)", "") if h else ""),
                        (h.get("Antibody", "") if h else "")] +
                       [sets[MAIN][k].get(g, 0) for k in SOURCES] + ["|".join(srcs), len(srcs)])

    # ------------------------------------------------------------ write markdown
    L = []
    add = L.append
    n_hpa = len(hpa_ids)
    add("# HPA × HRA protein/gene overlap report\n")
    add("_Part of [cns-iu/hra-hpa-overlap-report](https://github.com/cns-iu/hra-hpa-overlap-report). See the "
        "[project README](../README.md) for setup and the per-source gene maps._\n")
    add(f"_Generated {date.today().isoformat()} by "
        "[`src/build_overlap_report.py`](../src/build_overlap_report.py) from the maps in the "
        "[`output/`](../output/) folder. Genes are compared by HGNC ID. Unless stated otherwise, numbers use the "
        "**corrected** gene sets (suggested fixes applied; see [Sensitivity](#sensitivity-to-data-quality-fixes))._\n")

    # Executive summary
    add("## Executive summary\n")
    pro_in = len(asctb_prot & hpa_ids)
    gene_in = len(asctb_gene & hpa_ids)
    noncoding_out = sum(v for k, v in out_locus.items() if k != "protein-coding gene")
    add(f"- **The HPA covers {fmt_pct(len(hra_in), len(hra))} of the genes the HRA references.** "
        f"{len(hra_in):,} of the {len(hra):,} distinct genes named in ASCT+B, OMAP and AVR are in the "
        f"HPA gene table.")
    add(f"- **Almost every protein-level marker is in the HPA.** ASCT+B protein markers: "
        f"{fmt_pct(pro_in, len(asctb_prot))} ({pro_in:,}/{len(asctb_prot):,}); OMAP antibody targets: "
        f"{fmt_pct(len(S['omap'] & hpa_ids), len(S['omap']))}; AVR antibody targets: "
        f"{fmt_pct(len(S['avr'] & hpa_ids), len(S['avr']))}. ASCT+B gene (RNA) markers: "
        f"{fmt_pct(gene_in, len(asctb_gene))}.")
    add(f"- **The {len(hra_out):,} HRA genes missing from the HPA are mostly not protein-coding.** "
        f"{noncoding_out:,} of them are non-coding RNAs, pseudogenes or other non-protein loci "
        f"({', '.join(f'{k}: {v:,}' for k, v in out_locus.most_common(3))}). The HPA only covers "
        f"protein-coding genes, so this gap is expected, and it comes almost entirely from ASCT+B "
        f"gene (RNA) markers. {coding_out_text()}")
    add(f"- **The HRA uses {fmt_pct(len(hra_in), n_hpa)} of the HPA.** Only {len(hra_in):,} of the "
        f"{n_hpa:,} HPA genes with an HGNC ID appear in any HRA resource. That leaves a large pool of "
        f"protein-coding genes with HPA evidence that could become HRA biomarkers.")
    add(f"- **The HRA's antibody targets have unusually strong HPA antibody support.** "
        f"{fmt_pct(ab_strong, len(ab_in))} of OMAP/AVR targets in the HPA have *Enhanced* or "
        f"*Supported* HPA tissue staining, against {fmt_pct(hpa_strong, n_hpa)} of all HPA genes.")
    add(f"- **CTann gene sets are inferred, not stated.** CTann crosswalks only reference cell "
        f"types. Their {len(S['ctann']):,} genes come from ASCT+B, of which "
        f"{fmt_pct(len(S['ctann'] & hpa_ids), len(S['ctann']))} are in the HPA.")
    add(f"- **The data-quality fixes barely move the headline numbers** (HRA coverage "
        f"{fmt_pct(len(hra_v['given'] & hpa_ids), len(hra_v['given']))} as given, "
        f"{fmt_pct(len(hra_v['corrected'] & hpa_ids), len(hra_v['corrected']))} corrected, "
        f"{fmt_pct(len(hra_v['confident'] & hpa_ids), len(hra_v['confident']))} without flagged "
        f"rows). They do matter for individual genes: see [`review-items.csv`](../output/review-items.csv).\n")
    add(md_table(["", "Distinct genes", "In HPA", "% in HPA"],
                 [[f"**{label[k]}**", f"{n:,}", f"{i:,}", fmt_pct(i, n)] for k, n, i in cov_rows]
                 + [["**HRA union** (ASCT+B ∪ OMAP ∪ AVR)", f"{len(hra):,}", f"{len(hra_in):,}",
                     fmt_pct(len(hra_in), len(hra))],
                    ["**HPA** (genes with an HGNC ID)", f"{n_hpa:,}", "–", "–"]],
                 ["---", "---:", "---:", "---:"]))
    add("")

    # How
    add("## How the comparison works\n")
    add("```mermaid\nflowchart LR\n"
        "  subgraph SRC[Sources]\n"
        "    HPA[HPA gene table]\n"
        "    ASCTB[ASCT+B biomarkers]\n"
        "    OMAP[OMAP antibodies]\n"
        "    AVR[AVR antibodies]\n"
        "    CTANN[CTann cell types]\n"
        "    DC[DeepCell Types channels]\n"
        "  end\n"
        "  HGNC[(HGNC complete set)]\n"
        "  HPA -- Ensembl ID --> HGNC\n"
        "  ASCTB -- HGNC ID or label --> HGNC\n"
        "  OMAP -- HGNC ID --> HGNC\n"
        "  AVR -- HGNC ID + PDF check --> HGNC\n"
        "  DC -- curated marker table --> HGNC\n"
        "  CTANN -- CL term --> ASCTB\n"
        "  HGNC --> MAPS[*-hgnc-map.csv]\n"
        "  MAPS --> CMP{Compare gene sets<br/>by HGNC ID}\n"
        "  CMP --> REPORT[This report +<br/>hpa-hra-overlap.csv]\n"
        "```\n")
    add("- **A gene** is an HGNC ID. Multi-gene records count every gene they name: pan-cytokeratin, "
        "for example, counts each keratin it detects.")
    add(f"- **The HPA set** is the {n_hpa:,} HPA genes that resolve to an HGNC ID. The "
        f"{hpa_without_hgnc:,} HPA entries without an HGNC record cannot overlap anything.")
    add("- **The HRA union** is ASCT+B ∪ OMAP ∪ AVR. CTann is left out because its genes are a "
        "subset of ASCT+B by construction. DeepCell Types is not an HRA resource and is shown "
        "for comparison.\n")

    # Coverage
    add("## How much of each HRA resource is in the HPA\n")
    add("```mermaid\nxychart-beta horizontal\n"
        "  title \"Genes found in the HPA (%)\"\n"
        f"  x-axis [{', '.join(chr(34) + label[k] + chr(34) for k, _, _ in cov_rows)}, "
        "\"ASCT+B protein\", \"ASCT+B gene\"]\n"
        "  y-axis \"% of the resource's genes\" 0 --> 100\n"
        f"  bar [{', '.join(f'{pct(i, n):.1f}' for _, n, i in cov_rows)}, "
        f"{pct(pro_in, len(asctb_prot)):.1f}, {pct(gene_in, len(asctb_gene)):.1f}]\n"
        "```\n")
    by_type_rows = []
    for t_name in sorted(asctb_types[MAIN], key=lambda x: -len(asctb_types[MAIN][x])):
        g = asctb_types[MAIN][t_name]
        by_type_rows.append([t_name, f"{len(g):,}", f"{len(g & hpa_ids):,}", fmt_pct(len(g & hpa_ids), len(g))])
    add("ASCT+B by biomarker type (one gene can be both a gene and a protein marker):\n")
    add(md_table(["Biomarker type", "Distinct genes", "In HPA", "% in HPA"], by_type_rows,
                 ["---", "---:", "---:", "---:"]))
    add("")

    # Missing
    add("## HRA genes that are not in the HPA\n")
    add(f"{len(hra_out):,} HRA genes are missing from the HPA, grouped below by HGNC locus group. "
        "The HPA profiles protein-coding genes, so non-coding RNAs and pseudogenes are expected "
        f"to be missing. {coding_out_text()}\n")
    add("```mermaid\npie showData title HRA genes missing from the HPA, by locus group\n" +
        "\n".join(f'  "{k}" : {v}' for k, v in out_locus.most_common()) + "\n```\n")
    add(md_table(["Locus group", "Genes"], [[k, f"{v:,}"] for k, v in out_locus.most_common()],
                 ["---", "---:"]))
    add("\nMost-used missing genes (by number of ASCT+B biomarker rows):\n")
    add(md_table(["Gene", "HGNC ID", "Locus type", "ASCT+B rows", "Also in"],
                 [[hgnc.by_id.get(g, {}).get("symbol", g), g, hgnc.by_id.get(g, {}).get("locus_type", ""),
                   f"{asctb_rows_per_gene.get(g, 0):,}",
                   ", ".join(label[k] for k in ("omap", "avr") if g in S[k]) or "–"]
                  for g in top_missing], ["---", "---", "---", "---:", "---"]))
    add("")

    # HPA usage
    add("## How much of the HPA the HRA uses\n")
    add("```mermaid\npie showData title HPA genes referenced by the HRA\n"
        f"  \"In ASCT+B, OMAP or AVR\" : {len(hra_in)}\n"
        f"  \"Not referenced\" : {n_hpa - len(hra_in)}\n```\n")
    add(md_table(["Resource", "HPA genes referenced", "% of HPA"],
                 [[label[k], f"{i:,}", fmt_pct(i, n_hpa)] for k, n, i in cov_rows] +
                 [["**HRA union**", f"**{len(hra_in):,}**", f"**{fmt_pct(len(hra_in), n_hpa)}**"]],
                 ["---", "---:", "---:"]))
    add("")

    # Overlap between resources
    add("## Which resources share which genes\n")
    add("Each column is one exact combination of resources, e.g. genes in ASCT+B and OMAP but in "
        "neither AVR nor DeepCell. Bars are split by HPA membership. CTann is left out because "
        "its genes are a subset of ASCT+B.\n")
    add(picture("figures/upset", "UpSet plot of gene overlap between ASCT+B, OMAP, AVR and DeepCell Types, split by HPA membership"))
    add("")
    add(md_table(["Combination", "Genes", "In HPA", "Not in HPA", "Most-used examples"],
                 [[" + ".join(label[UPSET_KEYS[m]] for m in members), f"{a + b:,}", f"{a:,}", f"{b:,}",
                   ", ".join(hgnc.by_id[g]["symbol"] for g in sorted(
                       exact, key=lambda g: (-sum(sets[MAIN][k].get(g, 0) for k in UPSET_KEYS),
                                             hgnc.by_id[g]["symbol"]))[:6])
                   + (" …" if len(exact) > 6 else "")]
                  for members, a, b, exact in combos], ["---", "---:", "---:", "---:", "---"]))
    add("")
    only_ab = S["omap"] | S["avr"]
    add(f"- {len(only_ab - S['asctb']):,} OMAP/AVR antibody targets are **not** ASCT+B biomarkers. "
        "They may be candidates to add to the ASCT+B tables.")
    add(f"- {len((S['asctb'] & hpa_ids) - only_ab):,} ASCT+B genes in the HPA have no OMAP/AVR antibody. "
        "These could be candidates for new antibody panels, especially those with strong HPA "
        "antibody support (see the gene table).\n")

    # Reliability
    add("## HPA antibody evidence for HRA genes\n")
    add("The HPA rates its tissue-staining (IHC) antibodies as *Enhanced*, *Supported*, *Approved* "
        "or *Uncertain*. The chart compares those ratings for each HRA resource's genes (in the "
        "HPA) against all HPA genes.\n")
    add(picture("figures/hpa-reliability", "Stacked bars of HPA tissue IHC reliability for each resource's genes"))
    add("")
    add(md_table(["Gene set"] + RELIABILITY + [NO_IHC, "n"],
                 [[lab] + [f"{fmt_pct(c, sum(counts))}" for c in counts] + [f"{sum(counts):,}"]
                  for lab, counts in rel_rows], ["---"] + ["---:"] * 6))
    add("")
    add(f"{ab_has_hpa_ab:,} of the {len(ab_in):,} OMAP/AVR targets in the HPA have at least one HPA "
        "antibody listed, so independent HPA staining is available to cross-check most HRA antibody "
        "panels.\n")

    # Evidence
    add("### HPA protein evidence level\n")
    ev_keys = sorted(set(ev_all) | set(ev_hra), key=lambda k: -ev_all.get(k, 0))
    add(md_table(["HPA evidence", "HRA genes", "% of HRA genes in HPA", "% of all HPA"],
                 [[k or "(none)", f"{ev_hra.get(k, 0):,}", fmt_pct(ev_hra.get(k, 0), len(hra_in)),
                   fmt_pct(ev_all.get(k, 0), n_hpa)] for k in ev_keys],
                 ["---", "---:", "---:", "---:"]))
    add("")

    # Protein classes
    add("## What kinds of proteins the HRA uses\n")
    add("These are HPA protein classes for HRA genes in the HPA, compared with all HPA genes. "
        "Enrichment above 1× means the class is over-represented in the HRA. A gene can belong to "
        "several classes.\n")
    add(md_table(["HPA protein class", "HRA genes", "% of HRA genes", "% of all HPA", "Enrichment"],
                 cls_rows, ["---", "---:", "---:", "---:", "---:"]))
    add("")

    # Sensitivity
    add("## Sensitivity to data-quality fixes\n")
    add("The maps flag records that look wrong (see [`review-items.csv`](../output/review-items.csv)). "
        "The table shows genes in the HPA / distinct genes (% in HPA) for each way of handling them:\n")
    add("- **Given:** genes exactly as each map resolved them.")
    add("- **Corrected:** rows with a `suggested_hgnc_symbol` use the suggested gene. Used everywhere else in this report.")
    add("- **Confident:** every row flagged `needs_review` is dropped.\n")
    add(md_table(["Resource", "Given", "Corrected", "Confident"], sens, ["---", "---:", "---:", "---:"]))
    add("")
    add(f"Corrections add {len(gained):,} genes to the HRA union "
        f"({', '.join(sorted(hgnc.by_id[g]['symbol'] for g in gained)) or 'none'}) and remove "
        f"{len(lost):,} that were only present because of an error "
        f"({', '.join(sorted(hgnc.by_id[g]['symbol'] for g in lost)) or 'none'}).\n")

    # Outputs
    add("## Files\n")
    add(md_table(["File", "Contents"], [
        ["[`hpa-hra-overlap.csv`](hpa-hra-overlap.csv)",
         "One row per gene in the HPA or any source: HPA evidence, IHC reliability and antibodies, "
         "plus the number of rows in each source"],
        ["[`figures/`](figures/)", "Chart images (light and dark versions)"],
        ["[`review-items.csv`](../output/review-items.csv)", "Data-quality issues for expert review"],
        ["[`../output/*-hgnc-map.csv`](../output/)", "Per-source gene maps (see the [project README](../README.md))"],
    ]))
    add("\nRegenerate with `bash scripts/20-generate-maps.sh && bash scripts/50-overlap-report.sh`.")

    (report_dir / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"Wrote {report_dir / 'README.md'}, {report_dir / 'hpa-hra-overlap.csv'} and figures: "
          f"HRA union {len(hra):,} genes, {len(hra_in):,} in HPA "
          f"({fmt_pct(len(hra_in), len(hra))})", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--maps-dir", type=Path, default=Path("output"))
    ap.add_argument("--report-dir", type=Path, default=Path("overlap-report"))
    ap.add_argument("--hpa", type=Path, default=Path("downloads/proteinatlas.tsv.zip"))
    ap.add_argument("--hgnc", type=Path, default=Path("downloads/hgnc_complete_set.txt"))
    build(ap.parse_args())


if __name__ == "__main__":
    main()
