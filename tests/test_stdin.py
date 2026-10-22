"""`-` = standard input, in `coal` (Python engine), `coal-cpp` and the `pal2nal` drop-in.

Every case runs through each engine as a subprocess with the input bytes on stdin; the output must be byte-identical
to the same command line with the file path, and the two engines must give the same exit code, stdout and stderr
(and the same --report JSON). TMPDIR points into the test's own folder, so the check that no spooled copy
(.coal_in_*) is left behind sees only this test's files.

Run:  python -m unittest tests.test_stdin   (C++ checks need COAL_CPP_BIN or ./build/coal-cpp; else skipped)
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIX = ROOT / "tests" / "fixtures"
TIMEOUT = 60
TWICE = "'-' (standard input) can be given only once: as the protein alignment or as one DNA file"


def cpp_binary() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / (p.name + ".exe")):
        if cand.is_file() and (os.name != "nt" or cand.suffix.lower() == ".exe"):
            return cand
    return None


CPP = cpp_binary()


class StdinCase(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory(prefix="coal_stdin_")
        self.d = Path(self._td.name) / "work"
        self.tmp = Path(self._td.name) / "tmp"  # TMPDIR of the runs: spooled copies go here
        self.d.mkdir()
        self.tmp.mkdir()
        for case in ("clean", "clustal_input", "long_mismatch"):
            shutil.copytree(FIX / case, self.d / case)
        p = (FIX / "clean" / "pep.fasta").read_bytes()
        n = (FIX / "clean" / "nuc.fasta").read_bytes()
        # CRLF line ends and a UTF-8 BOM: read from `-` exactly as from a file
        (self.d / "crlf_pep.fa").write_bytes(b"\xef\xbb\xbf" + p.replace(b"\n", b"\r\n"))
        (self.d / "crlf_nuc.fa").write_bytes(b"\xef\xbb\xbf" + n.replace(b"\n", b"\r\n"))
        recs = n.split(b">")[1:]  # the DNA split over two files
        (self.d / "nuc_a.fa").write_bytes(b">" + recs[0])
        (self.d / "nuc_b.fa").write_bytes(b"".join(b">" + r for r in recs[1:]))
        (self.d / "dup.fa").write_bytes(b">s1\nMKV-L\n>s1\nMKVAL\n")

    def tearDown(self):
        self._td.cleanup()

    def env(self):
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), COAL_ENGINE="python", TMPDIR=str(self.tmp),
                   TEMP=str(self.tmp), TMP=str(self.tmp))
        env.pop("COAL_IMPL", None)
        return env

    def engines(self):
        yield "python", [sys.executable, "-m", "coal"]
        if CPP is not None:
            yield "cpp", [str(CPP)]

    def run_cmd(self, cmd, args, stdin: bytes | None):
        r = subprocess.run([*cmd, *args], cwd=self.d, env=self.env(), input=stdin if stdin is not None else b"",
                           capture_output=True, timeout=TIMEOUT)
        return r.returncode, r.stdout, r.stderr.decode("utf-8", "replace").replace("\r\n", "\n")

    def leftovers(self):
        return sorted(p.name for p in self.tmp.iterdir()) + sorted(
            p.name for p in self.d.rglob(".coal*"))

    def both(self, args, stdin):
        """(rc, stdout, stderr) of every engine; asserts they agree, returns the Python one."""
        res = {name: self.run_cmd(cmd, args, stdin) for name, cmd in self.engines()}
        if "cpp" in res:
            self.assertEqual(res["cpp"], res["python"], f"engines differ for {args}")
        self.assertEqual(self.leftovers(), [], f"temp files left by {args}")
        return res["python"]


class TestStdinMatchesFile(StdinCase):
    CASES = [  # (arguments with the file, which argument goes to stdin)
        (["clean/pep.fasta", "clean/nuc.fasta"], 0),
        (["clean/pep.fasta", "clean/nuc.fasta"], 1),
        (["clean/pep.fasta", "clean/nuc.fasta", "-output", "clustal"], 0),
        (["clean/pep.fasta", "clean/nuc.fasta", "-output", "codon"], 1),
        (["clean/pep.fasta", "clean/nuc.fasta", "-output", "paml", "-nogap"], 0),
        (["clustal_input/pep.aln", "clustal_input/nuc.fasta", "-output", "clustal"], 0),  # CLUSTAL from stdin
        (["clustal_input/pep.aln", "clustal_input/nuc.fasta", "-blockonly"], 1),
        (["long_mismatch/pep.fasta", "long_mismatch/nuc.fasta"], 0),  # with per-codon warnings
        (["long_mismatch/pep.fasta", "long_mismatch/nuc.fasta", "-nomismatch"], 1),
        (["crlf_pep.fa", "crlf_nuc.fa"], 0),  # CRLF + BOM as from a file
        (["crlf_pep.fa", "crlf_nuc.fa"], 1),
        (["clean/pep.fasta", "nuc_a.fa", "nuc_b.fa"], 1),  # `-` as one of several DNA files
        (["clean/pep.fasta", "nuc_a.fa", "nuc_b.fa"], 2),
        (["--", "clean/pep.fasta", "clean/nuc.fasta"], 1),  # after --, `-` is still standard input
    ]

    def test_piped_input_gives_the_file_output(self):
        for args, k in self.CASES:
            with self.subTest(args=args, stdin=args[k]):
                want = self.both(args, None)
                self.assertEqual(want[0], 0, want[2])
                data = (self.d / args[k]).read_bytes()
                got = self.both(args[:k] + ["-"] + args[k + 1:], data)
                self.assertEqual(got, want)

    def test_output_file(self):
        for name, cmd in self.engines():
            with self.subTest(engine=name):
                rc, out, err = self.run_cmd(cmd, ["-", "clean/nuc.fasta", "-o", "o.fa"],
                                            (self.d / "clean/pep.fasta").read_bytes())
                self.assertEqual((rc, out), (0, b""), err)
                ref = self.run_cmd(cmd, ["clean/pep.fasta", "clean/nuc.fasta"], None)[1]
                self.assertEqual((self.d / "o.fa").read_bytes(), ref)
                (self.d / "o.fa").unlink()


class TestStdinErrors(StdinCase):
    def test_two_dashes_is_a_usage_error(self):
        for args in (["-", "-"], ["-", "clean/nuc.fasta", "-"], ["clean/pep.fasta", "-", "-"]):
            with self.subTest(args=args):
                rc, out, err = self.both(args, b">s1\nMKV\n")
                self.assertEqual((rc, out), (2, b""))
                self.assertTrue(err.startswith("usage: coal "), err)
                self.assertEqual(err.splitlines()[-1], f"coal: error: {TWICE}")

    def test_errors_name_standard_input(self):
        cases = [
            (["-", "clean/nuc.fasta"], (self.d / "dup.fa").read_bytes(),
             "coal: error: standard input: duplicate sequence ID 's1': every sequence needs its own ID"),
            (["-", "clean/nuc.fasta"], b"",
             "coal: error: standard input: no sequences found: the file is empty or has no '>' header lines "
             "(expected a protein alignment in FASTA or CLUSTAL format)"),
            (["clean/pep.fasta", "-"], b">s1\nATGAAAGTTCTG\n>\nATG\n",
             "coal: error: standard input:3: empty FASTA header: a '>' line without a sequence ID (every record needs a "
             "name after '>')"),
            (["-", "clean/nuc.fasta"], (self.d / "clean/nuc.fasta").read_bytes(), None),  # swapped files
        ]
        for args, data, want in cases:
            with self.subTest(args=args, data=data[:20]):
                rc, out, err = self.both(args, data)
                self.assertEqual((rc, out), (1, b""), err)
                self.assertNotIn(".coal_in_", err)
                self.assertNotIn("'-'", err)
                self.assertIn("standard input", err)
                if want:
                    self.assertIn(want, err.splitlines())

    @unittest.skipIf(os.name == "nt", "POSIX: closing fd 0 in the child")
    def test_closed_stdin(self):
        for name, cmd in self.engines():
            with self.subTest(engine=name):
                r = subprocess.run([*cmd, "-", "clean/nuc.fasta"], cwd=self.d, env=self.env(), capture_output=True,
                                   timeout=TIMEOUT, preexec_fn=lambda: os.close(0))
                self.assertEqual(r.returncode, 1, r.stderr)
                self.assertEqual(r.stderr, b"coal: error: cannot read standard input: it is closed\n")


class TestStdinResumeAndReport(StdinCase):
    def test_resume_is_a_notice_and_a_normal_run(self):
        ref = self.both(["clean/pep.fasta", "clean/nuc.fasta"], None)[1]
        rc, out, err = self.both(["clean/pep.fasta", "-", "-o", "o.fa", "--resume"],
                                 (self.d / "clean/nuc.fasta").read_bytes())
        self.assertEqual(rc, 0, err)
        self.assertIn("coal: notice: o.fa cannot be resumed (an input is not a regular file); starting a normal run",
                      err.splitlines())
        self.assertEqual((self.d / "o.fa").read_bytes(), ref)
        self.assertEqual(sorted(p.name for p in self.d.iterdir() if "coal" in p.name), [])

    def test_report_records_dash_and_the_data_hash(self):
        data = (self.d / "clean/pep.fasta").read_bytes()
        reports = {}
        for name, cmd in self.engines():
            rc, out, err = self.run_cmd(cmd, ["-", "clean/nuc.fasta", "-o", "o.fa", "--report", "r.json"], data)
            self.assertEqual(rc, 0, err)
            rep = json.loads((self.d / "r.json").read_text())
            rep.pop("runtime_seconds")
            reports[name] = rep
        rep = reports["python"]
        prot = rep["inputs"][0]
        self.assertEqual(prot, {"role": "protein", "path": "-", "sha256": hashlib.sha256(data).hexdigest(),
                                "sequences": 3})
        self.assertEqual(rep["inputs"][1]["path"], "clean/nuc.fasta")
        self.assertEqual(rep["command"], "coal - clean/nuc.fasta -o o.fa")
        self.assertEqual(rep["pal2nal_equivalent"], "pal2nal.pl - clean/nuc.fasta -output fasta")
        if "cpp" in reports:
            self.assertEqual(reports["cpp"], rep)

    @unittest.skipUnless(os.path.exists("/dev/stdin"), "no /dev/stdin (Windows)")
    def test_report_hashes_piped_dev_stdin_data(self):
        # the Python engine used to hash an empty re-read of the pipe; both now hash the data read
        data = (self.d / "clean/nuc.fasta").read_bytes()
        reports = {}
        for name, cmd in self.engines():
            rc, out, err = self.run_cmd(cmd, ["clean/pep.fasta", "/dev/stdin", "-o", "o.fa", "--report", "r.json"],
                                        data)
            self.assertEqual(rc, 0, err)
            rep = json.loads((self.d / "r.json").read_text())
            rep.pop("runtime_seconds")
            reports[name] = rep
        self.assertEqual(reports["python"]["inputs"][1]["path"], "/dev/stdin")
        self.assertEqual(reports["python"]["inputs"][1]["sha256"], hashlib.sha256(data).hexdigest())
        if "cpp" in reports:
            self.assertEqual(reports["cpp"], reports["python"])

    def test_report_of_a_failed_stdin_run(self):
        reports = {}
        for name, cmd in self.engines():
            rc, out, err = self.run_cmd(cmd, ["-", "clean/nuc.fasta", "--report", "r.json"], b"")
            self.assertEqual(rc, 1, err)
            rep = json.loads((self.d / "r.json").read_text())
            rep.pop("runtime_seconds")
            reports[name] = rep
        rep = reports["python"]
        self.assertEqual(rep["inputs"][0]["path"], "-")
        self.assertEqual(rep["inputs"][0]["sha256"], hashlib.sha256(b"").hexdigest())
        self.assertEqual(rep["errors"][0]["text"][:len("standard input: ")], "standard input: ")
        if "cpp" in reports:
            self.assertEqual(reports["cpp"], rep)

    def test_report_of_a_usage_error_has_no_report(self):
        for name, cmd in self.engines():
            with self.subTest(engine=name):
                rc, _, _ = self.run_cmd(cmd, ["-", "-", "--report", "r.json"], b"")
                self.assertEqual(rc, 2)
                self.assertFalse((self.d / "r.json").exists())


class TestStdinHelp(unittest.TestCase):
    def test_help_mentions_dash(self):
        env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), COLUMNS="80")
        r = subprocess.run([sys.executable, "-m", "coal", "--help"], env=env, capture_output=True, text=True,
                           timeout=TIMEOUT)
        # coal-cpp's help is generated from this parser (test_cpp_cli.test_help_text compares the whole text)
        helps = [r.stdout]
        if CPP is not None:
            helps.append(subprocess.run([str(CPP), "--help"], capture_output=True, text=True, timeout=TIMEOUT).stdout)
        for h in helps:
            self.assertIn("protein alignment (FASTA or CLUSTAL; - reads standard\n", h)
            self.assertIn("coding DNA (FASTA; one or more files; - reads standard\n", h)


class TestPal2nalDropInStdin(StdinCase):
    def drop_in(self, args, stdin):
        r = subprocess.run([sys.executable, "-m", "coal_compat", *args], cwd=self.d, env=self.env(),
                           input=stdin if stdin is not None else b"", capture_output=True, timeout=TIMEOUT)
        return r.returncode, r.stdout, r.stderr.decode("utf-8", "replace")

    def test_dash_as_protein_or_dna(self):
        for args, k in ((["clean/pep.fasta", "clean/nuc.fasta"], 0),
                        (["clean/pep.fasta", "clean/nuc.fasta", "-output", "paml"], 1),
                        (["clustal_input/pep.aln", "clustal_input/nuc.fasta", "-blockonly"], 0)):
            with self.subTest(args=args, k=k):
                want = self.drop_in(args, None)
                self.assertEqual(want[0], 0, want[2])
                got = self.drop_in(args[:k] + ["-"] + args[k + 1:], (self.d / args[k]).read_bytes())
                self.assertEqual(got, want)
                self.assertEqual(self.leftovers(), [])

    def test_two_dashes(self):
        rc, out, err = self.drop_in(["-", "-"], b"")
        self.assertEqual((rc, out), (1, b""))
        self.assertEqual(err, f"pal2nal (COAL): error: {TWICE}\n")

    def test_error_names_standard_input(self):
        rc, out, err = self.drop_in(["-", "clean/nuc.fasta"], (self.d / "dup.fa").read_bytes())
        self.assertEqual(rc, 1)
        self.assertIn("coal: error: standard input: duplicate sequence ID 's1'", err)


if __name__ == "__main__":
    unittest.main()
