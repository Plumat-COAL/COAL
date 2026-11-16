"""Whose files are in a `coal pipeline` output folder: the record of what COAL runs wrote there (work/outputs.json),
and the rules built on it (docs/pipeline.md, "One run per folder, and one run's files only").

A run replaces, and removes when it starts, only what an earlier COAL run is known to have written into DIR, never
a file of the user's:

  * work/outputs.json is written when a run starts, before its first output (so what a killed run left is known):
    "files", the names the run may write (temporary NAME.tmp ones included); "read_in_place", earlier outputs it
    reads where they are (`--pairs DIR/orthologs.tsv`, `--clusters DIR/clusters.tsv`: never removed at a start,
    still COAL's to write again); for the per-cluster folders codon_alignments/ and protein_alignments/,
    "<folder>_dir" (the run writes it) and "<folder>_files" (the .fasta names it wrote there).
  * The folder list is [] when the run starts, dropped just before the engine writes there and written again when
    the engine ends in any way but SIGKILL ([] when the folder is not there). No list (0.6.0, or a run killed while
    its engine wrote) means every .fasta with the engine's lower-case extension is the earlier run's.
  * Without a readable record (a version without it, work/ deleted) the outputs that a COAL report.json of this
    folder lists are the earlier run's; a record that cannot be read but begins as COAL writes it, or a work/ that
    holds only an empty logs/, is still COAL's mark on work/.
  * Refused before anything changes: an alignment folder the run writes that holds what no COAL run wrote there
    (--overwrite does not apply to a folder); an output name that exists and that no run recorded; a --report
    FILE that is not a COAL report; a work/ that no COAL run made (those three unless --overwrite); a work that is
    not a folder and a --report that is not a regular file (always).

OutputsRecord holds the record and these rules for one output folder and needs no engine or tool (unit tests:
tests/test_outputs_record.py); coal.dnds_pipeline (_Run) says what a run plans to write and when the engine starts
and ends. Internal: not part of the supported Python API.
"""
from __future__ import annotations

import json
import os
import re
import stat
from collections.abc import Callable, Iterable

from . import __version__

# Files a run writes into the output folder (besides work/): whatever an earlier run left of them is removed when
# a run starts, so that the folder never mixes two runs (a dnds.kaks or codon_alignments.fasta of a run with other
# options, the dnds.tsv.tmp of a killed run).
# Every table and alignment file is written as NAME.tmp and renamed when it is complete (_atomic, and the engine's
# --new-files), so the .tmp names are a run's files too: a run that was killed leaves them, the next one removes them.
RUN_TABLES = ("dnds.tsv", "dnds.kaks", "codon_alignments.fasta", "protein_alignments.fasta", "cluster_summary.tsv",
              "genes.tsv", "orthologs.tsv", "unpaired.tsv", "clusters.tsv", "realign.tsv", "failed.tsv",
              "report.json")
RUN_OUTPUTS = RUN_TABLES + tuple(n + ".tmp" for n in RUN_TABLES) + ("realign.tsv.body",)

# the folders of per-cluster alignment files (cluster mode): name -> the option that writes it
ALIGNMENT_DIRS = {"codon_alignments": "--keep-alignments", "protein_alignments": "--keep-protein-alignments"}


class OutputsRefused(Exception):
    """A run that would replace or mix with files that are not COAL's: nothing in the folder was changed."""


def starts_like_coal_record(path: str) -> bool:
    """True when `path` is a regular file that begins as COAL writes work/outputs.json ({"coal_version":). Used when
    json could not read that file: COAL writes it through a temporary file and renames it, so only a damaged disk or
    an edit cuts it short, and such a file is still COAL's mark on work/. Any other file of that name is not (#98)."""
    try:
        if not os.path.isfile(path):
            return False
        with open(path, "rb") as fh:
            return fh.read(16) == b'{"coal_version":'
    except OSError:
        return False


def is_coal_report(path: str) -> bool:
    """True when `path` is a JSON report that COAL wrote (a `coal pipeline` report, or `coal --report`'s: its
    "schema" is coal.<name>/<n>); what `--report FILE` may replace without --overwrite. A marker, not a
    validation: a JSON object with such a "schema" counts, whatever else it holds (the question is only whose file
    it is, and a report that COAL cut short, {"schema": ..., "report_truncated": true}, is COAL's too). A file
    larger than any report (256 MiB) is not read."""
    try:
        if not os.path.isfile(path) or os.path.getsize(path) > 256 << 20:
            return False
        with open(path, encoding="utf-8", errors="backslashreplace") as fh:
            rep = json.load(fh)
    except (OSError, ValueError):
        return False
    schema = rep.get("schema") if isinstance(rep, dict) else None
    return isinstance(schema, str) and re.fullmatch(r"coal\.[a-z_]+/\d+", schema) is not None


def check_report_path(report) -> None:
    """--report naming something that is not a regular file (also through a link): a folder, a named pipe, a
    device such as /dev/null. Refused before anything opens it (a named pipe would block the run: #98) and whatever
    --overwrite says (the report is written through FILE.tmp and renamed onto FILE)."""
    if not report:
        return
    try:
        st = os.stat(report)
    except OSError:
        return          # (missing, or a broken link: written as a new file)
    if stat.S_ISREG(st.st_mode):
        return
    if stat.S_ISDIR(st.st_mode):
        raise OutputsRefused(f"--report {report} is a folder: the report is written as a file of that name. Nothing "
                             f"was changed: name a file, e.g. {os.path.join(report, 'report.json')}")
    raise OutputsRefused(f"--report {report} is not a regular file (a device or a named pipe): the report is "
                         "written as a new file and renamed onto that name. Nothing was changed: name a regular file "
                         "(without --report it is DIR/report.json)")


def _and(items) -> str:
    """"a", "a and b", "a, b and c" (as dnds_pipeline._and)."""
    items = [str(x) for x in items]
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


class OutputsRecord:
    """The ownership record of one output folder `out_dir`: what the earlier run wrote there (read_earlier), the
    checks made before a run changes anything (check_foreign), the clean-up at its start (clear), and this run's
    record (new_record, engine_starts, list_alignment_files), which write() puts into work/outputs.json.

    The methods that change the record change it in memory and return whether it changed; the caller writes it
    (and decides what a failure to write means). `schema` is the "schema" of the reports whose outputs count when
    there is no record (coal.pipeline/1)."""

    def __init__(self, out_dir: str, schema: str):
        self.out = out_dir
        self.work = os.path.join(out_dir, "work")
        self.schema = schema
        # What read_earlier found:
        self.known: set = set()       # names (RUN_OUTPUTS) an earlier run is known to have written here
        self.dirs: set = set()        # alignment folders whose .fasta files are the earlier run's
        self.read_in_place: set = set()   # earlier outputs that a later run read in place (record: read_in_place)
        self.listed: set = set()      # the outputs a COAL report.json of this folder lists
        # the .fasta names an earlier run recorded for its codon_alignments/ / protein_alignments/ (absent: no list,
        # a run of 0.6.0 or one killed while it wrote them; every .fasta there is then that run's)
        self.dir_files: dict = {}
        self.work_is_coal = False     # DIR/work exists and an earlier COAL run made it
        self.record: dict = {}        # this run's record (new_record), as write() puts it into work/outputs.json

    def path(self, *parts) -> str:
        return os.path.join(self.out, *parts)

    def wpath(self, *parts) -> str:
        return os.path.join(self.work, *parts)

    @property
    def earlier(self) -> set:
        """Every name an earlier COAL run wrote here or read in place, as far as the record or a report says."""
        return self.known | self.read_in_place | self.listed

    # -- what an earlier run wrote
    def read_earlier(self, report=None, report_is_default: bool = False, writes: Iterable[str] = ()) -> None:
        """What an earlier run is known to have written here: the names in work/outputs.json (written when a run
        starts, so a killed run counts), or, for a run of a version without that file, the outputs its report.json
        lists; and the alignment folders whose .fasta files are that run's (a cluster run with --keep-alignments /
        --keep-protein-alignments). `report`: this run's --report (a COAL report of this very folder there counts
        too, unless it is DIR/report.json itself: report_is_default). `writes`: the alignment folders this run
        writes; one that holds anything no COAL run wrote there is refused (OutputsRefused)."""
        names: set = set()
        dirs: set = set()
        codon_dir = False
        # The reports of earlier runs into this folder: DIR/report.json, and the --report file of this run when it
        # is a COAL report of this very folder (an earlier run with the same --report elsewhere). What they list
        # is COAL's for check_foreign when work/outputs.json is gone (work/ deleted, an older version) or has no
        # read_in_place (a development build's record).
        reports = []
        for path, here_only in ((self.path("report.json"), False), (report, True)):
            if not path or (here_only and report_is_default) or not os.path.isfile(path):
                continue   # (not a regular file: a named pipe would block open(); #98)
            try:
                with open(path, encoding="utf-8", errors="backslashreplace") as fh:
                    rep = json.load(fh)
            except (OSError, ValueError):
                continue
            if not (isinstance(rep, dict) and rep.get("schema") == self.schema
                    and isinstance(rep.get("outputs"), dict)):
                continue
            if here_only and not (isinstance(rep.get("out_dir"), str)
                                  and os.path.realpath(rep["out_dir"]) == os.path.realpath(self.out)):
                continue
            reports.append((rep, not here_only))
        self.listed = {v for rep, _ in reports for v in rep["outputs"].values()
                       if isinstance(v, str) and v in RUN_TABLES}
        try:
            if not os.path.isfile(self.wpath("outputs.json")):
                raise FileNotFoundError(self.wpath("outputs.json"))
            with open(self.wpath("outputs.json"), encoding="utf-8", errors="backslashreplace") as fh:
                rec = json.load(fh)
            names = {n for n in rec["files"] if isinstance(n, str)}
            read = rec.get("read_in_place")
            if isinstance(read, list):
                self.read_in_place = {n for n in read if isinstance(n, str)}
                # this record says what was read in place: a report's list is not needed (and would make a
                # user's own file read in place, which its report lists, look like COAL's)
                self.listed = set()
            codon_dir = rec.get("codon_alignments_dir") is True
            if rec.get("protein_alignments_dir") is True:
                dirs.add("protein_alignments")
            for name in ALIGNMENT_DIRS:
                listed = rec.get(name + "_files")
                if isinstance(listed, list):
                    self.dir_files[name] = {n for n in listed if isinstance(n, str)}
            self.work_is_coal = True
        except (OSError, ValueError, KeyError, TypeError):
            # A work/outputs.json that cannot be read (cut short, edited) is still COAL's mark on work/ (COAL writes
            # it through a temporary file; a user's work/ has none), and so is a work/ that holds only an empty
            # logs/ (a run killed between making work/logs and writing the record): #98
            w = self.work
            try:
                self.work_is_coal = (starts_like_coal_record(self.wpath("outputs.json"))
                                     or (os.listdir(w) == ["logs"] and os.path.isdir(self.wpath("logs"))
                                         and not os.path.islink(self.wpath("logs"))
                                         and not os.listdir(self.wpath("logs"))))
            except OSError:
                pass
            # no record (a version without it, or work/ deleted): the outputs the first of those reports lists
            if reports:
                rep, in_folder = reports[0]
                names = {v for v in rep["outputs"].values() if isinstance(v, str)}
                names |= {"report.json", "failed.tsv"} if in_folder else {"failed.tsv"}
                codon_dir = rep.get("mode") == "clusters" and "codon_alignments" in names
                names |= {n + ".tmp" for n in list(names)}
                # work/ of a version without work/outputs.json: the folder's own report.json says a COAL run was here
                self.work_is_coal = self.work_is_coal or in_folder
        if codon_dir:
            dirs.add("codon_alignments")
        writes = set(writes)
        for name, opt in ALIGNMENT_DIRS.items():
            d = self.path(name)
            if not (os.path.isdir(d) and not os.path.islink(d)):
                dirs.discard(name)
            if not (name in writes and os.path.lexists(d)):
                continue
            if name not in dirs:
                # not COAL's: the engine would write its files among (or over) what the folder holds
                if not os.path.isdir(d) or os.path.islink(d) or os.listdir(d):
                    raise OutputsRefused(f"{d} exists and was not written by a COAL run: {opt} writes one file per "
                                         "cluster there and would mix with (or replace) what it holds. Nothing was "
                                         "changed: move it away or use another output folder (--overwrite does not "
                                         "apply to a folder)")
                continue
            mine = self.foreign_alignment_files(name)
            if mine:
                # COAL's folder, with .fasta files of the user's among the earlier run's: kept at the start, but
                # this run writes a file per cluster there and one may have the same name (#96)
                listed = ", ".join(mine[:3]) + (f" and {len(mine) - 3:,} more" if len(mine) > 3 else "")
                raise OutputsRefused(f"{d} holds {listed}, which no COAL run wrote there: {opt} writes one file per "
                                     f"cluster there and could replace {'it' if len(mine) == 1 else 'them'}. Nothing "
                                     f"was changed: move {'it' if len(mine) == 1 else 'them'} out of {d}, or use "
                                     "another output folder (--overwrite does not apply to a folder)")
        self.known = {n for n in names if n in RUN_OUTPUTS}
        self.dirs = dirs

    def foreign_alignment_files(self, name: str) -> list[str]:
        """The .fasta files in the alignment folder `name` (COAL's) that the earlier run's record does not list: the
        user's own. A folder recorded without a list (0.6.0, or work/ deleted: report.json names the folder only)
        has none but those whose .fasta is not in the engine's lower case."""
        listed = self.dir_files.get(name)
        d = self.path(name)
        if not os.path.isdir(d) or os.path.islink(d):
            return []
        # .fasta in any case (X.FASTA): on a case-insensitive file system the engine's X.fasta would replace it
        # (#98). Without a list the earlier run's are the files with the engine's own .fasta.
        return sorted(n for n in os.listdir(d) if n.lower().endswith(".fasta")
                      and (n not in listed if listed is not None else not n.endswith(".fasta")))

    # -- before anything changes
    def check_foreign(self, planned: Iterable[str], *, overwrite: bool = False, pairs_in_place: bool = False,
                      report=None, report_is_default: bool = False) -> None:
        """Refuse (OutputsRefused), before anything in the folder is removed or written, a run whose outputs
        (`planned`, names in DIR) would replace files that no COAL run wrote there: a genes.tsv, dnds.tsv,
        report.json or clusters.tsv of the user's own (they were replaced without a word, also by
        `-o examples/dnds`, whose orthologs.tsv is part of the repository: #87). A file an earlier run recorded
        (read_earlier) is that run's and is replaced as before; so is everything with --overwrite. Temporary names
        (NAME.tmp) are COAL's own and not checked. Then, the same way (#96): a --report FILE that exists and is not
        a COAL report, and a DIR/work/ that no COAL run made (not empty, no work/outputs.json, no COAL report.json
        in DIR); a `work` that is not a folder is refused always."""
        w = self.work
        if os.path.lexists(w) and not os.path.isdir(w):
            # (a file or a broken link: COAL never removes it, so --overwrite cannot help)
            raise OutputsRefused(f"{w} exists and is not a folder: the run keeps its working files in DIR/work. "
                                 "Nothing was changed: move it away or use another output folder")
        if overwrite:
            return
        known = self.earlier
        foreign = []
        for name in planned:
            if name.endswith((".tmp", ".body")) or name in known:
                continue
            if name == "orthologs.tsv" and pairs_in_place:
                continue    # --pairs DIR/orthologs.tsv: read and kept as it is
            if os.path.lexists(self.path(name)):
                foreign.append(name)
        if foreign:
            raise OutputsRefused(
                f"{self.out} already holds {_and(foreign)}, which this run would replace and which no COAL run "
                f"wrote there (work/outputs.json does not list {'it' if len(foreign) == 1 else 'them'}). Nothing "
                "was changed. Use a new -o folder, or add --overwrite to replace "
                f"{'it' if len(foreign) == 1 else 'them'}")
        if report and not report_is_default and os.path.lexists(report) and not is_coal_report(report):
            # --report onto a file of the user's, in DIR or anywhere (#96)
            raise OutputsRefused(f"--report {report} already exists and is not a COAL report: this run would replace "
                                 "it. Nothing was changed. Name another --report file, or add --overwrite to "
                                 "replace it")
        if os.path.isdir(w) and not self.work_is_coal and os.listdir(w):
            # a work/ of the user's (#96): the run would write its logs, state and copies there
            raise OutputsRefused(f"{w} exists and no COAL run made it (it has no work/outputs.json): this run keeps "
                                 "its working files there (logs, state.json, decompressed inputs) and would mix "
                                 "them with what it holds. Nothing was changed. Use a new -o folder, move it away, "
                                 "or add --overwrite to write into it (files there with COAL's own names are "
                                 "replaced; nothing else is touched)")

    def clear(self, keep: Iterable[str] = (), is_input: Callable[[str], object] = lambda p: None) -> None:
        """Remove the files an earlier COAL run wrote into this folder (read_earlier), so the folder holds this run's
        results only: never one in `keep` (what --resume reuses), never one for which is_input(path) is true (this
        run's inputs), and never a file COAL has no record of (a `clusters.tsv` or a .fasta of the user's in
        codon_alignments/ or protein_alignments/ stays). An alignment folder left empty is removed. OSError on the
        first file that cannot be removed."""
        keep = set(keep)
        for name in sorted(self.known):
            p = self.path(name)
            if name not in keep and (os.path.isfile(p) or os.path.islink(p)) and not is_input(p):
                os.remove(p)
        for d in (self.path(name) for name in sorted(self.dirs)):
            if not os.path.isdir(d) or os.path.islink(d):
                continue
            listed = self.dir_files.get(os.path.basename(d))
            for n in os.listdir(d):
                p = os.path.join(d, n)
                # a .fasta the earlier run recorded (every .fasta when it recorded none), and the .fasta.tmp
                # files a run that was killed left; a user's own .fasta in the folder stays (#96)
                theirs = n.endswith(".fasta.tmp") or (n.endswith(".fasta") and (listed is None or n in listed))
                if theirs and not is_input(p):
                    os.remove(p)
            if not os.listdir(d):
                os.rmdir(d)

    # -- this run's record
    def new_record(self, planned: Iterable[str], writes: Iterable[str] = ()) -> dict:
        """This run's record: `planned`, the names it may write (without an output it reads in place), and
        `writes`, the alignment folders it writes. Call after clear()."""
        planned = list(planned)
        writes = set(writes)
        # An earlier COAL run's output that this run reads where it is (`--pairs DIR/orthologs.tsv`, `--clusters
        # DIR/clusters.tsv`): not in "files", so that no later run removes it at its start (it may have been edited),
        # but still COAL's for check_foreign, so that a later run may write it again.
        # Carried forward from run to run while the file is there and no run writes it (a run that neither reads
        # nor writes it keeps the entry: #87). After clear(), an earlier output that this run does not write is
        # still there only when it is one of this run's inputs, or was read in place.
        read = sorted(n for n in self.earlier
                      if n in RUN_TABLES and n not in planned and os.path.isfile(self.path(n)))
        self.record = {"coal_version": __version__, "files": sorted(planned), "read_in_place": read,
                       "codon_alignments_dir": "codon_alignments" in writes,
                       "protein_alignments_dir": "protein_alignments" in writes}
        # The alignment folders this run writes hold no .fasta now (the earlier run's were removed, and one of
        # anyone else's refused the run): an empty list, so that a run that fails before its engine makes no file
        # there the run's (#98). engine_starts() drops the list just before the engine writes there (a run killed
        # while it writes then left files COAL cannot name: all of them are the run's, as in 0.6.0), and
        # list_alignment_files() lists what the engine wrote when it ends in any way but SIGKILL.
        for name in ALIGNMENT_DIRS:
            if self.record[name + "_dir"]:
                self.record[name + "_files"] = []
        return self.record

    def write(self) -> None:
        """Put the record into work/outputs.json (through outputs.json.tmp, renamed: never cut short)."""
        tmp = self.wpath("outputs.json.tmp")
        with open(tmp, "w", encoding="utf-8", errors="backslashreplace", newline="\n") as fh:
            json.dump(self.record, fh)
            fh.write("\n")
        os.replace(tmp, self.wpath("outputs.json"))

    def engine_starts(self) -> bool:
        """Cluster mode, just before coal-cpp starts writing the per-cluster files: the record of the folders it
        writes has no list from now until list_alignment_files() makes one, so that what a run killed meanwhile
        (SIGKILL, out of memory, power loss) wrote counts as the run's, every .fasta of the folder, as for a folder
        of 0.6.0: the next run removes them, and `--resume --keep-alignments` goes on (#98). The folders hold
        nothing of anyone else's at this point (read_earlier refused that). True when the record changed."""
        dropped = False
        for name in ALIGNMENT_DIRS:
            if self.record.pop(name + "_files", None) is not None:
                dropped = True
        return dropped

    def list_alignment_files(self) -> bool:
        """Cluster mode, after the engine: the .fasta files this run wrote into the alignment folders it writes, by
        name, so that the next run removes those and no others (a file of the user's put there later stays: #96).
        The folder held no .fasta of anyone else's when the engine started (read_earlier refuses that), so what it
        holds now is this run's. A folder that is not there (it could not be made) gets an empty list: with none, a
        folder of that name that the user makes later would be taken for the run's, every .fasta in it removed by
        the next run. True when the record changed. OSError when a folder cannot be read."""
        changed = False
        for name in ALIGNMENT_DIRS:
            if not self.record.get(name + "_dir"):
                continue
            d = self.path(name)
            self.record[name + "_files"] = (sorted(n for n in os.listdir(d) if n.endswith(".fasta"))
                                            if os.path.isdir(d) and not os.path.islink(d) else [])
            changed = True
        return changed
