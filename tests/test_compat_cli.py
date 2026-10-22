"""The `pal2nal` drop-in (coal_compat): option translation, the command, and what its docs claim."""
import os
import subprocess
import sys
import unittest

from coal_compat.cli import UsageError, translate

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(ROOT, "examples", "clean")


class TestTranslate(unittest.TestCase):
    def test_options_any_order(self):  # DEC-4: -output always passed, files after "--"
        self.assertEqual(translate(["-nogap", "p.aln", "-codontable", "2", "n.fa", "-output", "fasta"]),
                         ["-output", "fasta", "-codontable", "2", "-nogap", "--", "p.aln", "n.fa"])

    def test_two_dna_files(self):
        self.assertEqual(translate(["p", "a.fa", "b.fa", "-output", "fasta"]),
                         ["-output", "fasta", "--", "p", "a.fa", "b.fa"])

    def test_default_is_clustal_and_every_format_passes(self):  # DEC-4
        self.assertEqual(translate(["p", "n"]), ["-output", "clustal", "--", "p", "n"])
        for f in ("clustal", "paml", "fasta", "codon"):
            self.assertEqual(translate(["p", "n", "-output", f]), ["-output", f, "--", "p", "n"])
        self.assertEqual(translate(["p", "n", "-blockonly"]), ["-output", "clustal", "-blockonly", "--", "p", "n"])

    def test_clear_errors(self):
        for argv, text in ((["p", "n", "-html"], "html"),
                           (["p", "n", "-foo"], "unknown option"),
                           (["p", "n", "-output", "xml"], "unknown -output"),
                           (["p", "n", "-codontable", "x", "-output", "fasta"], "number"),
                           (["p", "-output", "fasta"], "at least one DNA"),
                           (["p", "n"], "cannot write 'clustal'")):  # DEC-4: only if an engine lacks a format
            with self.assertRaises(UsageError) as cm:
                translate(argv, ("fasta",))
            self.assertIn(text, str(cm.exception))


class TestCommand(unittest.TestCase):
    def test_runs_coal_to_stdout(self):
        r = subprocess.run([sys.executable, "-m", "coal_compat", os.path.join(EX, "protein.aln.fasta"),
                            os.path.join(EX, "dna.fasta"), "-output", "fasta"],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(r.stdout.startswith(">gene01"))

    def test_help(self):
        r = subprocess.run([sys.executable, "-m", "coal_compat", "-h"], capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0)
        self.assertIn("-codontable", r.stdout)

    @unittest.skipIf(os.name == "nt", "a closed pipe is EPIPE on POSIX only (Windows reports EINVAL)")
    def test_usage_into_a_closed_pipe_is_141_and_silent(self):  # pal2nal | head
        # stdout is a pipe whose read end is already closed, so the first write fails with EPIPE (Python ignores
        # SIGPIPE, so it is a BrokenPipeError): no traceback, exit 141 as for `coal` (docs/usage.md).
        for argv in ([], ["-h"]):
            with self.subTest(argv=argv):
                rfd, wfd = os.pipe()
                os.close(rfd)
                try:
                    r = subprocess.run([sys.executable, "-m", "coal_compat", *argv], stdout=wfd,
                                       stderr=subprocess.PIPE, text=True, timeout=60)
                finally:
                    os.close(wfd)
                self.assertEqual(r.returncode, 141, r.stderr)
                self.assertEqual(r.stderr, "")


class TestCompatClaims(unittest.TestCase):
    def test_docstring_and_html_error_match_what_works(self):  # clustal is supported
        from coal_compat import cli
        doc = " ".join(cli.__doc__.split())
        self.assertIn("An unchanged PAL2NAL command line works", doc)
        self.assertIn("default output format is clustal", doc)
        with self.assertRaises(cli.UsageError) as cm:
            cli.translate(["p", "n", "-html"])
        self.assertIn("not supported", str(cm.exception))

    def test_nostderr_is_accepted_as_quiet(self):
        from coal_compat import cli
        self.assertEqual(cli.translate(["p", "n", "-output", "fasta", "-nostderr"]),
                         ["-output", "fasta", "-q", "--", "p", "n"])


if __name__ == "__main__":
    unittest.main()
