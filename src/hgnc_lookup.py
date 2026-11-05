"""Shared HGNC lookups used by the *_hgnc_map build scripts."""
import csv
import re
from collections import defaultdict
from pathlib import Path

ENSEMBL_RE = re.compile(r"ENSG\d{11}", re.IGNORECASE)
CD_RE = re.compile(r"CD\d+[A-Z]?")
# Expression-level qualifiers appended to marker names, e.g. "CD11c+/lo", "CD33(lo)", "CD64hi".
QUALIFIER_RE = re.compile(r"(\s*(\(?(hi|lo|high|low|dim|bright|int)\)?|[+\-/]))+$", re.IGNORECASE)
DEFAULT_WITHDRAWN = Path("downloads/hgnc_withdrawn.txt")
# Words in antibody names that are also HGNC symbols but name reagents, dyes or clones.
NOT_GENE_WORDS = {"KIT", "APC", "SP6", "SP7", "SET", "CAT", "MAX", "REST", "IMPACT", "CLOCK"}


# Common protein/marker names that HGNC does not list as aliases, with the genes
# they may legitimately refer to (curated; keys are normalize_name() forms).
PROTEIN_NAMES = {
    "CD3": {"CD3D", "CD3E", "CD3G"},
    "IGM": {"IGHM"}, "IGD": {"IGHD"}, "IGG1": {"IGHG1"}, "IGG3": {"IGHG3"},
    "HLADR": {"HLA-DRA", "HLA-DRB1"},
    "SMA": {"ACTA2"}, "ASMA": {"ACTA2"}, "SMOOTHMUSCLEACTIN": {"ACTA2"},
    "ALPHASMOOTHMUSCLEACTIN": {"ACTA2"}, "ACTINSMOOTHMUSCLE": {"ACTA2"},
    "CD45RA": {"PTPRC"}, "CD45RO": {"PTPRC"},
    "ECADHERIN": {"CDH1"}, "NCADHERIN": {"CDH2"}, "VECADHERIN": {"CDH5"},
    "PDGFRALPHA": {"PDGFRA"}, "PDGFRBETA": {"PDGFRB"},
    "CALRETININ": {"CALB2"}, "PERFORIN": {"PRF1"}, "TRYPTASE": {"TPSAB1", "TPSB2"},
    "CALDESMON": {"CALD1"}, "BETACATENIN": {"CTNNB1"}, "D240": {"PDPN"},
    "COLLAGENI": {"COL1A1", "COL1A2"}, "COLLAGENIII": {"COL3A1"},
    "COLLAGENIV": {"COL4A1", "COL4A2"}, "COLIV": {"COL4A1", "COL4A2"},
}


def normalize_name(name):
    """Upper-case, drop a trailing +/- marker sign and all non-alphanumerics."""
    return re.sub(r"[^A-Z0-9]", "", name.strip().upper().rstrip("+-"))


def split_list(value):
    return [v for v in value.split("|") if v]


class Hgnc:
    def __init__(self, path, withdrawn_path=DEFAULT_WITHDRAWN):
        self.by_id, self.by_symbol, self.by_ensembl = {}, {}, {}
        self.by_prev, self.by_alias = defaultdict(list), defaultdict(list)
        with open(path, encoding="utf-8", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["status"] != "Approved":
                    continue
                self.by_id[r["hgnc_id"]] = r
                self.by_symbol[r["symbol"].upper()] = r
                if r["ensembl_gene_id"]:
                    self.by_ensembl[r["ensembl_gene_id"]] = r
                for s in split_list(r["prev_symbol"]):
                    self.by_prev[s.upper()].append(r)
                for s in split_list(r["alias_symbol"]):
                    self.by_alias[s.upper()].append(r)
        self.withdrawn = {}
        if withdrawn_path and Path(withdrawn_path).exists():
            with open(withdrawn_path, encoding="utf-8", newline="") as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    self.withdrawn[r["HGNC_ID"]] = r

    def match_label(self, label):
        """Return (records, match type) for a free-text marker label."""
        key = label.strip().upper()
        if key in self.by_symbol:
            return [self.by_symbol[key]], "symbol"
        m = ENSEMBL_RE.search(label)
        if m and m.group(0).upper() in self.by_ensembl:
            return [self.by_ensembl[m.group(0).upper()]], "ensembl"
        for table, name in ((self.by_prev, "prev_symbol"), (self.by_alias, "alias_symbol")):
            hits = table.get(key)
            if hits:
                return hits, name if len(hits) == 1 else name + "_multiple"
        return [], ""

    def describe_missing(self, hgnc_id):
        """Explain why an HGNC ID is not an approved record (withdrawn, merged, unknown)."""
        w = self.withdrawn.get(hgnc_id)
        if w is None:
            return f"{hgnc_id} is not a known HGNC ID"
        merged = w["MERGED_INTO_REPORT(S) (i.e HGNC_ID|SYMBOL|STATUS)"].strip()
        if merged:
            return f"{hgnc_id} ({w['WITHDRAWN_SYMBOL']}) was withdrawn and merged into {merged}"
        return f"{hgnc_id} ({w['WITHDRAWN_SYMBOL']}) was withdrawn by HGNC"

    def is_known_name(self, label, record):
        """True if label is the record's symbol, previous/alias symbol or full name.

        Comparison ignores case, punctuation and a trailing +/- (e.g. "CD4+", "MAP-2").
        """
        names = [record["symbol"], record["name"]]
        names += split_list(record["prev_symbol"]) + split_list(record["alias_symbol"])
        names = {normalize_name(n) for n in names}
        keys = {normalize_name(label), normalize_name(QUALIFIER_RE.sub("", label.strip()))}
        if any(k and k in names for k in keys):
            return True
        if any(record["symbol"] in PROTEIN_NAMES.get(k, ()) for k in keys):
            return True
        return record["hgnc_id"] in self.genes_named_in(label).values()

    def truncation_candidates(self, hgnc_id, label):
        """Approved genes whose ID is hgnc_id plus one digit and whose symbol starts with label."""
        key = normalize_name(label)
        out = []
        for d in "0123456789":
            r = self.by_id.get(hgnc_id + d)
            if r and key and normalize_name(r["symbol"]).startswith(key):
                out.append(r)
        return out

    def genes_named_in(self, text):
        """Genes an antibody or marker name mentions by approved symbol or CD name.

        Text in brackets or parentheses (usually clone IDs) is ignored, and only
        exact approved symbols or unambiguous CD designations count, to keep false
        positives low. Returns {token: hgnc_id}.
        """
        text = re.sub(r"\[[^\]]*\]|\([^)]*\)", " ", text)
        # "EGF Receptor", "MARCKS like protein": the named gene is not the target.
        text = re.sub(r"[A-Za-z0-9-]+\s+(receptor|like)\b", " ", text, flags=re.IGNORECASE)
        found = {}
        for tok in re.split(r"[^A-Za-z0-9\-]+", text):
            for t in {tok, *tok.split("-")}:
                u = t.upper()
                if len(u) < 3 or u.isdigit() or u in NOT_GENE_WORDS:
                    continue
                if u in self.by_symbol:
                    found[u] = self.by_symbol[u]["hgnc_id"]
                elif CD_RE.fullmatch(u):
                    hits = self.by_alias.get(u, []) + self.by_prev.get(u, [])
                    if len({h["hgnc_id"] for h in hits}) == 1:
                        found[u] = hits[0]["hgnc_id"]
        return found
