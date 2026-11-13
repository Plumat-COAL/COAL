"""What to cite for what was run (`coal cite`, the `cite` list of the pipeline's report.json, docs/kaks.md).

COAL's codon alignment follows PAL2NAL's method; its Ka/Ks engine is derived from KaKs_Calculator 2.0, PAML's yn00,
libsequence and seqinr and implements published estimators; `coal pipeline` runs DIAMOND and Clustal Omega / MAFFT.
A paper that reports numbers from them should cite the methods and programs those numbers come from, not only COAL.

REFERENCES: key -> the reference as text (authors, year, title, journal, volume and pages). DOIS: key -> the DOI of
that work (each checked against the publisher's or the institution's page of the paper; COAL itself has none).
keys_for(...) gives the keys for a run; the keys are stable (they are written into reports): a key that is replaced
by a later publication stays in REFERENCES and is named in SUPERSEDED (buchfink2023, the DeepClust preprint ->
buchfink2026, the paper), so a report written earlier still resolves.
"""
from __future__ import annotations

__all__ = ["REFERENCES", "DOIS", "SUPERSEDED", "LABELS", "METHOD_KEYS", "TABLE", "keys_for", "reference", "doi_url",
           "format_references", "table", "main"]

COAL_KEY = "coal"
REFERENCES: dict[str, str] = {
    "coal": "Newton C, Murati K (2026). COAL: codon alignments from protein alignments, version {version}. "
            "https://github.com/plumat-coal/COAL",
    "suyama2006": "Suyama M, Torrents D, Bork P (2006). PAL2NAL: robust conversion of protein sequence alignments "
                  "into the corresponding codon alignments. Nucleic Acids Research 34:W609-W612.",
    "nei1986": "Nei M, Gojobori T (1986). Simple methods for estimating the numbers of synonymous and "
               "nonsynonymous nucleotide substitutions. Molecular Biology and Evolution 3:418-426.",
    "li1985": "Li W-H, Wu C-I, Luo C-C (1985). A new method for estimating synonymous and nonsynonymous rates of "
              "nucleotide substitution considering the relative likelihood of nucleotide and codon changes. "
              "Molecular Biology and Evolution 2:150-174.",
    "li1993": "Li W-H (1993). Unbiased estimation of the rates of synonymous and nonsynonymous substitution. "
              "Journal of Molecular Evolution 36:96-99.",
    "pamilo1993": "Pamilo P, Bianchi NO (1993). Evolution of the Zfx and Zfy genes: rates and interdependence "
                  "between the genes. Molecular Biology and Evolution 10:271-281.",
    "tzeng2004": "Tzeng Y-H, Pan R, Li W-H (2004). Comparison of three methods for estimating rates of synonymous "
                 "and nonsynonymous nucleotide substitutions. Molecular Biology and Evolution 21:2290-2298.",
    "yang2000": "Yang Z, Nielsen R (2000). Estimating synonymous and nonsynonymous substitution rates under "
                "realistic evolutionary models. Molecular Biology and Evolution 17:32-43.",
    "zhang2006": "Zhang Z, Li J, Yu J (2006). Computing Ka and Ks with a consideration of unequal transitional "
                 "substitutions. BMC Evolutionary Biology 6:44.",
    "comeron1995": "Comeron JM (1995). A method for estimating the numbers of synonymous and nonsynonymous "
                   "substitutions per site. Journal of Molecular Evolution 41:1152-1159.",
    "wang2010": "Wang D, Zhang Y, Zhang Z, Zhu J, Yu J (2010). KaKs_Calculator 2.0: a toolkit incorporating "
                "gamma-series methods and sliding window strategies. Genomics, Proteomics & Bioinformatics "
                "8:77-80.",
    "yang2007": "Yang Z (2007). PAML 4: phylogenetic analysis by maximum likelihood. Molecular Biology and "
                "Evolution 24:1586-1591.",
    "thornton2003": "Thornton K (2003). libsequence: a C++ class library for evolutionary genetic analysis. "
                    "Bioinformatics 19:2325-2327.",
    "charif2007": "Charif D, Lobry JR (2007). SeqinR 1.0-2: a contributed package to the R project for "
                  "statistical computing devoted to biological sequences retrieval and analysis. In: Bastolla U, "
                  "Porto M, Roman HE, Vendruscolo M (eds), Structural Approaches to Sequence Evolution: Molecules, "
                  "Networks, Populations. Springer, pp. 207-232.",
    "buchfink2021": "Buchfink B, Reuter K, Drost H-G (2021). Sensitive protein alignments at tree-of-life scale "
                    "using DIAMOND. Nature Methods 18:366-368.",
    "buchfink2026": "Buchfink B, Barb\u00e9 \u00c9, Ashkenazy H, Reuter K, Kennedy JA, Drost H-G (2026). Clustering the "
                    "protein universe of life using DIAMOND DeepClust. Nature Methods 23:724-727.",
    # the preprint of buchfink2026: not in TABLE; kept so that reports written before the paper still resolve
    "buchfink2023": "Buchfink B, Ashkenazy H, Reuter K, Kennedy JA, Drost H-G (2023). Sensitive clustering of "
                    "protein sequences at tree-of-life scale using DIAMOND DeepClust. bioRxiv (the preprint of "
                    "buchfink2026).",
    "sievers2011": "Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, "
                   "Soding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein "
                   "multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539.",
    "katoh2013": "Katoh K, Standley DM (2013). MAFFT multiple sequence alignment software version 7: "
                 "improvements in performance and usability. Molecular Biology and Evolution 30:772-780.",
    "drost2015": "Drost H-G, Gabel A, Grosse I, Quint M (2015). Evidence for active maintenance of "
                 "phylotranscriptomic hourglass patterns in animal and plant embryogenesis. Molecular Biology "
                 "and Evolution 32:1221-1231. (the R package orthologr)",
}

# The DOI of each cited work (none for COAL itself). Printed as https://doi.org/ links.
DOIS: dict[str, str] = {
    "suyama2006": "10.1093/nar/gkl315",
    "nei1986": "10.1093/oxfordjournals.molbev.a040410",
    "li1985": "10.1093/oxfordjournals.molbev.a040343",
    "li1993": "10.1007/BF02407308",
    "pamilo1993": "10.1093/oxfordjournals.molbev.a040003",
    "tzeng2004": "10.1093/molbev/msh242",
    "yang2000": "10.1093/oxfordjournals.molbev.a026236",
    "zhang2006": "10.1186/1471-2148-6-44",
    "comeron1995": "10.1007/BF00173196",
    "wang2010": "10.1016/S1672-0229(10)60008-3",
    "yang2007": "10.1093/molbev/msm088",
    "thornton2003": "10.1093/bioinformatics/btg316",
    "charif2007": "10.1007/978-3-540-35306-5_10",
    "buchfink2021": "10.1038/s41592-021-01101-x",
    "buchfink2026": "10.1038/s41592-026-03030-z",
    "buchfink2023": "10.1101/2023.01.24.525373",
    "sievers2011": "10.1038/msb.2011.75",
    "katoh2013": "10.1093/molbev/mst010",
    "drost2015": "10.1093/molbev/msv012",
}

# A key that a later publication replaces -> the key to cite now
SUPERSEDED: dict[str, str] = {"buchfink2023": "buchfink2026"}


def doi_url(key: str) -> str | None:
    """https://doi.org/<DOI> of a reference, None for COAL itself."""
    return "https://doi.org/" + DOIS[key] if key in DOIS else None


def reference(key: str, version: str | None = None) -> str:
    """The reference as one line of text, with its DOI link. KeyError for an unknown key."""
    text = REFERENCES[key]
    if "{version}" in text:
        if version is None:
            from . import __version__ as version
        text = text.replace("{version}", str(version))
    url = doi_url(key)
    return f"{text} {url}" if url else text


def current(keys) -> list[str]:
    """The keys with each superseded one followed by the key that replaces it; each key once."""
    out: list[str] = []
    for k in keys:
        for one in (k, SUPERSEDED.get(k)):
            if one and one not in out:
                out.append(one)
    return out


# How a reference is named in running text and in the tables of the documentation
LABELS: dict[str, str] = {
    "coal": "COAL", "suyama2006": "Suyama et al. 2006 (PAL2NAL)", "nei1986": "Nei & Gojobori 1986",
    "li1985": "Li, Wu & Luo 1985", "li1993": "Li 1993", "pamilo1993": "Pamilo & Bianchi 1993",
    "tzeng2004": "Tzeng, Pan & Li 2004", "yang2000": "Yang & Nielsen 2000", "zhang2006": "Zhang, Li & Yu 2006",
    "comeron1995": "Comeron 1995", "wang2010": "Wang et al. 2010 (KaKs_Calculator 2.0)",
    "yang2007": "Yang 2007 (PAML 4)", "thornton2003": "Thornton 2003 (libsequence)",
    "charif2007": "Charif & Lobry 2007 (seqinr)", "buchfink2021": "Buchfink et al. 2021 (DIAMOND)",
    "buchfink2026": "Buchfink et al. 2026 (DIAMOND DeepClust)",
    "buchfink2023": "Buchfink et al. 2023 (DIAMOND DeepClust, preprint)",
    "sievers2011": "Sievers et al. 2011 (Clustal Omega)", "katoh2013": "Katoh & Standley 2013 (MAFFT)",
    "drost2015": "Drost et al. 2015 (orthologr)",
}


def markdown_table() -> str:
    """The table for README.md and docs/kaks.md (tests/test_citations.py checks that they hold it)."""
    rows = ["| What you ran | Cite |", "|---|---|"]
    rows += [f"| {what} | {'; '.join(LABELS[k] for k in keys)} |" for what, keys in TABLE]
    return "\n".join(rows)


def listed_keys() -> list[str]:
    """The cited works shown in the documents: every reference but COAL's own and the superseded ones."""
    return [k for k in REFERENCES if k != COAL_KEY and k not in SUPERSEDED]


def markdown_references() -> str:
    return "\n".join(f"- {REFERENCES[k]} <{doi_url(k)}>" for k in listed_keys())


# The estimator of each Ka/Ks method (docs/kaks.md "Methods"), then the program whose code COAL's implementation
# of it is derived from.
METHOD_KEYS: dict[str, tuple[str, ...]] = {
    "NG": ("nei1986", "wang2010"),
    "LWL": ("li1985", "wang2010"),
    "MLWL": ("li1985", "tzeng2004", "wang2010"),
    "LPB": ("li1993", "pamilo1993", "wang2010"),
    "MLPB": ("li1993", "pamilo1993", "tzeng2004", "wang2010"),
    "YN": ("yang2000", "wang2010", "yang2007"),
    "MYN": ("zhang2006", "yang2000", "wang2010", "yang2007"),
    "Comeron": ("comeron1995", "thornton2003"),
    "Li": ("li1993", "charif2007"),
}

# (what was run, keys): the table of `coal cite`, README.md, docs/kaks.md and site/cite.html
TABLE: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("Anything", ("coal",)),
    ("A codon alignment (`coal PEP NUC`, `pal2nal`, every Ka/Ks route of the pipeline)", ("suyama2006",)),
    ("Ka/Ks, method NG", METHOD_KEYS["NG"]),
    ("Ka/Ks, method LWL", METHOD_KEYS["LWL"]),
    ("Ka/Ks, method MLWL", METHOD_KEYS["MLWL"]),
    ("Ka/Ks, method LPB", METHOD_KEYS["LPB"]),
    ("Ka/Ks, method MLPB", METHOD_KEYS["MLPB"]),
    ("Ka/Ks, method YN (the default)", METHOD_KEYS["YN"]),
    ("Ka/Ks, method MYN", METHOD_KEYS["MYN"]),
    ("Ka/Ks, method Comeron", METHOD_KEYS["Comeron"]),
    ("Ka/Ks, method Li", METHOD_KEYS["Li"]),
    ("`coal pipeline pairs` with DIAMOND (orthologs)", ("buchfink2021",)),
    ("`coal pipeline clusters` with DIAMOND DeepClust", ("buchfink2021", "buchfink2026")),
    ("Protein alignments by Clustal Omega (`--aligner clustalo`, the cluster default)", ("sievers2011",)),
    ("Protein alignments by MAFFT (`--aligner mafft`)", ("katoh2013",)),
    ("Through the R package orthologr", ("drost2015",)),
)


def keys_for(methods=(), codon_alignment: bool = True, diamond: bool = False, deepclust: bool = False,
             clustalo: bool = False, mafft: bool = False) -> list[str]:
    """The citation keys of a run, COAL first, each once, in a fixed order."""
    keys = [COAL_KEY]
    if codon_alignment:
        keys.append("suyama2006")
    for m in METHOD_KEYS:   # the order of the methods' rows
        if m in methods:
            keys += METHOD_KEYS[m]
    if diamond or deepclust:
        keys.append("buchfink2021")
    if deepclust:
        keys.append("buchfink2026")
    if clustalo:
        keys.append("sievers2011")
    if mafft:
        keys.append("katoh2013")
    seen: set[str] = set()
    return [k for k in keys if not (k in seen or seen.add(k))]


def format_references(keys, version: str | None = None) -> str:
    """One reference per line, `[key] text https://doi.org/...`."""
    return "\n".join(f"[{k}] " + reference(k, version) for k in keys)


def table() -> str:
    """The what-was-run -> what-to-cite table as plain text."""
    lines = ["What you ran -> what to cite (keys; the references follow)", ""]
    for what, keys in TABLE:
        lines.append(f"  {what.replace('`', '')}\n      {', '.join(keys)}")
    return "\n".join(lines)


USAGE = """usage: coal cite [REPORT.json | KEY ... | METHOD ...]

What to cite for what you ran. Without arguments: the table of what was run -> references, and every reference.
With the report.json of a `coal pipeline` run: the references of that run (its "cite" list).
With citation keys (suyama2006, yang2000, ...) or Ka/Ks method names (YN, NG, ...): those references.
Each cited work is printed with its DOI as a https://doi.org/ link. A report written before the DeepClust paper
appeared names its preprint (buchfink2023): the paper (buchfink2026) is printed after it.
"""


def main(argv) -> int:
    """`coal cite ...`."""
    import json
    import os
    import sys
    args = list(argv)
    if args and args[0] in ("-h", "--help"):
        print(USAGE, end="")
        return 0
    if not args:
        print(table())
        print()
        _print(format_references([COAL_KEY, *listed_keys()]))
        return 0
    keys: list[str] = []
    by_upper = {m.upper(): m for m in METHOD_KEYS} | {"C95": "Comeron"}
    for a in args:
        if a in REFERENCES:
            keys.append(a)
        elif a.upper() in by_upper:
            keys += ["coal", "suyama2006", *METHOD_KEYS[by_upper[a.upper()]]]
        elif os.path.isfile(a):
            try:
                with open(a, encoding="utf-8", errors="backslashreplace") as fh:
                    rep = json.load(fh)
                got = rep["cite"]
                if not isinstance(got, list) or not all(k in REFERENCES for k in got):
                    raise KeyError("cite")
            except (OSError, ValueError, KeyError, TypeError):
                print(f"coal cite: error: {a} is not a `coal pipeline` report with a \"cite\" list (reports of "
                      "COAL 0.6 and later have one)", file=sys.stderr)
                return 1
            keys += got
        else:
            print(f"coal cite: error: {a!r} is not a citation key, a Ka/Ks method or a report file. Keys: "
                  + ", ".join(REFERENCES) + ". Methods: " + ", ".join(METHOD_KEYS), file=sys.stderr)
            return 2
    _print(format_references(current(keys)))
    return 0


def _print(text: str) -> None:
    """print, on a terminal whose encoding lacks a letter of an author's name too."""
    import sys
    try:
        print(text)
    except UnicodeEncodeError:
        enc = getattr(sys.stdout, "encoding", None) or "ascii"
        print(text.encode(enc, "replace").decode(enc))
