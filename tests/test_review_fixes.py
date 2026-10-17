"""Regression tests for the independent engine review (overnight/engine-fixes/CHANGES.md, ENG-a .. ENG-i).

Every input is built inside a temporary directory created by the test; -o is only
ever pointed at paths inside that directory.
"""
from __future__ import annotations

import os
import random
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
TIMEOUT = 60
WINDOWS = os.name == "nt"  # R2-D1

PEP = ">s1\nMKV-AL\n>s2\nMKVSAL\n"
NUC = ">s1\nATGAAAGTTGCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n"
EXPECTED = ">s1\nATGAAAGTT---GCTCTG\n>s2\nATGAAGGTCTCAGCACTA\n"


def coal(args, cwd, stdin=None, timeout=TIMEOUT, **kw):
    env = dict(os.environ, PYTHONPATH=str(SRC))
    return subprocess.run([sys.executable, "-m", "coal", *args], cwd=cwd, env=env, input=stdin,
                          capture_output=True, text=True, timeout=timeout, **kw)


class TempDirCase(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory(prefix="coal_review_")
        self.d = Path(self._td.name)
        (self.d / "p.fa").write_text(PEP)
        (self.d / "n.fa").write_text(NUC)

    def tearDown(self):
        self._td.cleanup()

    def leftovers(self):
        return sorted(p.name for p in self.d.iterdir() if p.name.startswith(".coal_"))


# ---------------------------------------------------------------- ENG-a
class TestEngA_OutputPath(TempDirCase):
    def test_symlink_is_written_through(self):
        target = self.d / "target.fa"
        target.write_text("old\n")
        link = self.d / "link.fa"
        try:
            link.symlink_to(target)
        except (OSError, NotImplementedError):  # R2-D1
            self.skipTest("cannot create symlinks here")
        r = coal(["p.fa", "n.fa", "-o", "link.fa"], self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(link.is_symlink())
        self.assertEqual(target.read_text(), EXPECTED)
        self.assertEqual(self.leftovers(), [])

    @unittest.skipIf(WINDOWS, "R2-D1: no POSIX permission bits / umask on Windows")
    def test_new_file_mode_follows_umask(self):
        old = os.umask(0o022)
        try:
            r = coal(["p.fa", "n.fa", "-o", "out.fa"], self.d)
        finally:
            os.umask(old)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(stat.S_IMODE((self.d / "out.fa").stat().st_mode), 0o644)

    @unittest.skipUnless(hasattr(os, "mkfifo"), "needs FIFOs")
    def test_fifo_is_written_not_replaced(self):
        fifo = self.d / "out.fifo"
        os.mkfifo(fifo)
        reader = subprocess.Popen(["cat", str(fifo)], stdout=subprocess.PIPE, text=True)
        try:
            r = coal(["p.fa", "n.fa", "-o", str(fifo)], self.d)
            out, _ = reader.communicate(timeout=20)
        finally:
            if reader.poll() is None:
                reader.kill()
                reader.wait()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(out, EXPECTED)
        self.assertTrue(stat.S_ISFIFO(fifo.stat().st_mode))
        self.assertEqual(self.leftovers(), [])


@unittest.skipIf(WINDOWS, "R2-E5: POSIX permissions / links")
class TestR2E5_ExistingOutput(TempDirCase):  # R2-E5
    def test_existing_file_keeps_its_permissions(self):
        out = self.d / "out.fa"
        out.write_text("old\n")
        os.chmod(out, 0o600)
        old = os.umask(0o022)
        try:
            r = coal(["p.fa", "n.fa", "-o", "out.fa"], self.d)
        finally:
            os.umask(old)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(out.read_text(), EXPECTED)
        self.assertEqual(stat.S_IMODE(out.stat().st_mode), 0o600)

    def test_hard_link_is_written_in_place(self):
        out, other = self.d / "out.fa", self.d / "hl2.fa"
        out.write_text("old\n")
        os.link(out, other)
        r = coal(["p.fa", "n.fa", "-o", "out.fa"], self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(other.read_text(), EXPECTED)
        self.assertEqual(out.stat().st_ino, other.stat().st_ino)
        self.assertEqual(self.leftovers(), [])

    def test_symlink_loop_is_a_clear_error(self):
        a, b = self.d / "a.fa", self.d / "b.fa"
        a.symlink_to(b)
        b.symlink_to(a)
        r = coal(["p.fa", "n.fa", "-o", "a.fa"], self.d)
        self.assertEqual(r.returncode, 1, r.stdout)
        self.assertIn("symbolic link", r.stderr)
        self.assertTrue(a.is_symlink() and b.is_symlink())
        self.assertEqual(self.leftovers(), [])


class TestR2E6_BadByteLine(TempDirCase):  # R2-E6
    def test_latin1_byte_reports_its_line(self):
        lines = [">s1"] + ["MKVAL"] * 399 + ["MK\xe9AL"]
        (self.d / "lat.fa").write_bytes(("\n".join(lines) + "\n").encode("latin-1"))
        r = coal(["lat.fa", "n.fa"], self.d)
        self.assertEqual(r.returncode, 1)
        self.assertIn("lat.fa:401", r.stderr)
        self.assertNotIn("Traceback", r.stderr)


@unittest.skipUnless(os.path.exists("/dev/stdin"), "R2-E7: no /dev/stdin (Windows)")
class TestR2E7_PipeNames(TempDirCase):  # R2-E7
    def test_error_in_piped_input_names_the_argument(self):
        r = coal(["/dev/stdin", "n.fa"], self.d, stdin=">s1\nMKV-AL\n>s1\nMKVSAL\n")
        self.assertEqual(r.returncode, 1)
        self.assertIn("/dev/stdin", r.stderr)
        self.assertNotIn(".coal_in_", r.stderr)


class TestR3E2_OutputFolderErrors(TempDirCase):  # R3-E2
    def test_missing_folder_names_the_output_path(self):
        r = coal(["p.fa", "n.fa", "-o", os.path.join("nodir", "out.fa")], self.d)
        self.assertEqual(r.returncode, 1, r.stderr)
        self.assertNotIn(".coal_", r.stderr)
        self.assertIn(os.path.join("nodir", "out.fa"), r.stderr)
        self.assertIn("folder", r.stderr)
        self.assertIn("does not exist", r.stderr)
        self.assertFalse((self.d / "nodir").exists())

    def test_folder_symlink_loop_names_the_output_path(self):
        try:
            (self.d / "loopA").symlink_to(self.d / "loopB")
            (self.d / "loopB").symlink_to(self.d / "loopA")
        except (OSError, NotImplementedError):
            self.skipTest("cannot create symlinks here")
        r = coal(["p.fa", "n.fa", "-o", os.path.join("loopA", "out.fa")], self.d)
        self.assertEqual(r.returncode, 1, r.stderr)
        self.assertNotIn(".coal_", r.stderr)
        self.assertIn(os.path.join("loopA", "out.fa"), r.stderr)
        self.assertEqual(self.leftovers(), [])


class TestR3E3_DnaErrorLineNumbers(TempDirCase):  # R3-E3
    def test_data_before_header_has_line_number(self):
        (self.d / "bn.fa").write_text("# comment\n\nATGAAA\n>s1\nATGAAAGTTGCTCTG\n")
        r = coal(["p.fa", "bn.fa"], self.d)
        self.assertEqual(r.returncode, 1)
        self.assertIn("bn.fa:3: sequence data before the first '>' header", r.stderr)  # CHANGED(ENG-UX-11)

    def test_non_ascii_id_has_line_number(self):
        (self.d / "un.fa").write_bytes(b">s1\nATGAAAGTTGCTCTG\n>s\xc3\xa92\nATG\n")
        r = coal(["p.fa", "un.fa"], self.d)
        self.assertEqual(r.returncode, 1)
        self.assertIn("un.fa:3:", r.stderr)

    def test_crlf_and_cr_line_numbers(self):
        for name, nl in (("crlf.fa", b"\r\n"), ("cr.fa", b"\r")):
            (self.d / name).write_bytes(nl.join([b">s1", b"ATGAAAGTTGCTCTG", b">", b"ATG", b""]))
            r = coal(["p.fa", name], self.d)
            self.assertEqual(r.returncode, 1)
            self.assertIn(f"{name}:3: empty FASTA header", r.stderr)


# ---------------------------------------------------------------- ENG-b
@unittest.skipUnless(os.path.exists("/dev/stdin"), "R2-D1: no /dev/stdin (Windows)")
class TestEngB_StreamInputs(TempDirCase):
    def test_protein_from_pipe(self):
        r = coal(["/dev/stdin", "n.fa"], self.d, stdin=PEP)
        self.assertEqual((r.returncode, r.stdout), (0, EXPECTED), r.stderr)

    def test_dna_from_pipe(self):
        tmp = Path(tempfile.gettempdir())
        before = set(tmp.glob(".coal_in_*"))
        r = coal(["p.fa", "/dev/stdin"], self.d, stdin=NUC)
        self.assertEqual((r.returncode, r.stdout), (0, EXPECTED), r.stderr)
        self.assertEqual(set(tmp.glob(".coal_in_*")) - before, set())  # spooled copy removed

    @unittest.skipUnless(hasattr(os, "mkfifo"), "needs FIFOs")
    def test_fifo_inputs(self):  # what <(cat p.fa) gives
        os.mkfifo(self.d / "pf")
        os.mkfifo(self.d / "nf")
        writer = subprocess.Popen(["sh", "-c", "cat p.fa > pf & cat n.fa > nf & wait"], cwd=self.d)
        try:
            r = coal(["pf", "nf"], self.d)
        finally:
            try:
                writer.wait(timeout=5)
            except subprocess.TimeoutExpired:
                writer.kill()
                writer.wait()
        self.assertEqual((r.returncode, r.stdout), (0, EXPECTED), r.stderr)


# ---------------------------------------------------------------- ENG-c
def _random_pair(seed=3, n=60):
    rng = random.Random(seed)
    prot = "M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n - 1))
    dna = "".join(rng.choice("ACGT") for _ in range(3 * n))
    return f">g1\n{prot}\n>g2\n{prot}\n", f">g1\n{dna}\n>g2\n{dna}\n"


class TestEngC_Unrelated(TempDirCase):
    def test_unrelated_protein_and_dna_is_an_error(self):
        pep, nuc = _random_pair()
        (self.d / "rp.fa").write_text(pep)
        (self.d / "rn.fa").write_text(nuc)
        for extra in ([], ["-q"]):
            r = coal(["rp.fa", "rn.fa", "--strict", *extra], self.d)  # DEC-1
            self.assertNotEqual(r.returncode, 0, r.stdout)
            self.assertEqual(r.stdout, "")
            self.assertIn("g1: protein and DNA are inconsistent", r.stderr)

    def test_no_anchor_matched_warns_and_q_silences_it(self):  # DEC-2
        # Golden table3/fasta: PAL2NAL succeeds although no 10-residue anchor matches
        # (7 of 25 codons differ under table 1), so this stays a success, with a notice.
        fx = ROOT / "tests" / "fixtures" / "table3"
        gold = (ROOT / "tests" / "golden" / "table3" / "fasta.stdout").read_text()
        r = coal([str(fx / "pep.fasta"), str(fx / "nuc.fasta")], self.d)
        self.assertEqual((r.returncode, r.stdout), (0, gold), r.stderr)
        self.assertIn("s0: no 10-residue segment", r.stderr)
        self.assertIn("7 of 25", r.stderr)
        r = coal([str(fx / "pep.fasta"), str(fx / "nuc.fasta"), "-q"], self.d)
        self.assertEqual((r.returncode, r.stdout, r.stderr), (0, gold, ""))


class TestR2E3_AmbiguousCodons(TempDirCase):  # R2-E3
    def masked_pair(self, frac=0.55, n=300, seed=7):
        sys.path.insert(0, str(SRC))
        from coal.tables import get_code
        code = get_code(1).code
        sense = sorted(c for c, a in code.items() if a != "*" and c != "ATG")
        rng = random.Random(seed)
        cod = ["ATG"] + [rng.choice(sense) for _ in range(n - 1)]
        prot = "".join(code[c] for c in cod)
        masked = list(cod)
        for i in rng.sample(range(1, n), int(frac * n)):
            masked[i] = "NNN"
        return prot, "".join(masked)

    def test_mostly_masked_gene_is_aligned(self):
        prot, dna = self.masked_pair()
        clean_p, clean_d = "MKVAL", "ATGAAAGTTGCTCTG"
        pad = len(prot) - len(clean_p)
        (self.d / "mp.fa").write_text(f">masked\n{prot}\n>clean\n{clean_p}{'-' * pad}\n")
        (self.d / "mn.fa").write_text(f">masked\n{dna}\n>clean\n{clean_d}\n")
        r = coal(["mp.fa", "mn.fa", "-q"], self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        sys.path.insert(0, str(SRC))
        from coal.pipeline import fasta_record
        self.assertEqual(r.stdout, fasta_record("masked", dna) + fasta_record("clean", clean_d + "---" * pad))

    def test_all_ambiguous_codons_do_not_trigger_the_rule(self):
        (self.d / "ap.fa").write_text(">g\nMKVALWYFHIKL\n")
        (self.d / "an.fa").write_text(">g\nATG" + "NNN" * 11 + "\n")
        r = coal(["ap.fa", "an.fa", "-q"], self.d)
        self.assertEqual((r.returncode, r.stdout), (0, ">g\nATG" + "NNN" * 11 + "\n"), r.stderr)

    def test_definite_mismatches_still_rejected_with_checkable_wording(self):
        pep, nuc = _random_pair()
        (self.d / "rp.fa").write_text(pep)
        (self.d / "rn.fa").write_text(nuc)
        r = coal(["rp.fa", "rn.fa", "--strict"], self.d)  # DEC-1
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("checkable codons do not match", r.stderr)


class TestR3E1_MaskedGeneWithUTR(TempDirCase):  # R3-E1
    PROT = "MTVDCLAHRFISDLVKSPSFGSVLACYVISPLDRNVCGQH"
    DNA = ("GCATGNNNNNNNNNNNNTTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNATCAGTNNNNNNNNNNNN"
           "NNNNNNTGTNNNNNNNNN")

    def test_reviewer_repro_is_not_placed_shifted(self):
        (self.d / "ep.fa").write_text(f">g\n{self.PROT}\n")
        (self.d / "en.fa").write_text(f">g\n{self.DNA}\n")
        r = coal(["ep.fa", "en.fa", "--strict"], self.d)  # DEC-1
        self.assertNotEqual(r.returncode, 0, r.stdout)
        self.assertEqual(r.stdout, "")
        self.assertIn("g: protein and DNA are inconsistent", r.stderr)
        self.assertIn("too few", r.stderr)

    def test_randomized_heavy_mask_with_utr_never_shifted_at_exit_0(self):
        sys.path.insert(0, str(SRC))
        from coal.core import CoalError, thread_one
        from coal.tables import get_code
        code = get_code(1).code
        sense = sorted(c for c, a in code.items() if a != "*" and c != "ATG")
        rng = random.Random(2)  # at 3f9ee66 this seed gave 12 shifted outputs at exit 0
        for _ in range(3000):
            n = rng.randint(20, 120)
            cod = ["ATG"] + [rng.choice(sense) for _ in range(n - 1)]
            prot = "".join(code[c] for c in cod)
            for i in rng.sample(range(1, n), int(rng.uniform(0.70, 0.98) * (n - 1))):
                cod[i] = "NNN"
            utr = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 30)))
            cds = "".join(cod)
            try:
                got = thread_one("g", prot, utr + cds, 1, strict=True).dna  # DEC-1
            except CoalError:
                continue
            self.assertEqual(got, cds, (prot, utr, cds))
            # without the UTR the placement is forced: always accepted, unchanged
            self.assertEqual(thread_one("g", prot, cds, 1).dna, cds)

    def test_randomized_heavy_mask_exact_length_is_accepted(self):
        sys.path.insert(0, str(SRC))
        from coal.core import thread_one
        from coal.tables import get_code
        code = get_code(1).code
        sense = sorted(c for c, a in code.items() if a != "*" and c != "ATG")
        rng = random.Random(5)
        for _ in range(500):
            n = rng.randint(20, 120)
            cod = ["ATG"] + [rng.choice(sense) for _ in range(n - 1)]
            prot = "".join(code[c] for c in cod)
            for i in rng.sample(range(1, n), int(rng.uniform(0.70, 0.98) * (n - 1))):
                cod[i] = "NNN"
            self.assertEqual(thread_one("g", prot, "".join(cod), 1).dna, "".join(cod))

    def test_masked_gene_without_utr_and_trailing_stop_still_errors_or_is_right(self):
        sys.path.insert(0, str(SRC))
        from coal.core import CoalError, thread_one
        cds = "ATG" + "NNN" * 11
        try:
            self.assertEqual(thread_one("g", "MKVALWYFHIKL", cds + "TAA", 1).dna, cds)
        except CoalError:
            pass

    def test_notice_counts_checkable_codons(self):
        # golden table3: no anchor, 7 of 25 checkable codons differ, none ambiguous
        fx = ROOT / "tests" / "fixtures" / "table3"
        r = coal([str(fx / "pep.fasta"), str(fx / "nuc.fasta")], self.d)  # DEC-2: was with -q
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("(7 of 25 checkable codons do not match)", r.stderr)


# ---------------------------------------------------------------- ENG-d
class TestEngD_PartialIdMatch(TempDirCase):
    def test_some_ids_match_warns_and_pairs_by_order(self):
        (self.d / "ip.fa").write_text(">geneA\nMKVAL\n>geneB\nMKVAL\n")
        (self.d / "in.fa").write_text(">geneB\nATGAAAGTTGCTCTG\n>geneX\nATGAAGGTCGCACTA\n")
        r = coal(["ip.fa", "in.fa"], self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stdout, ">geneA\nATGAAAGTTGCTCTG\n>geneB\nATGAAGGTCGCACTA\n")
        self.assertIn("1 of 2 protein IDs", r.stderr)
        self.assertIn("'geneA' gets DNA 'geneB'", r.stderr)
        r = coal(["ip.fa", "in.fa", "-q"], self.d)  # DEC-2
        self.assertEqual((r.returncode, r.stdout, r.stderr), (0, ">geneA\nATGAAAGTTGCTCTG\n>geneB\nATGAAGGTCGCACTA\n", ""))

    def test_no_ids_match_stays_silent(self):  # golden id_mismatch_same_order: PAL2NAL pairs by order
        (self.d / "ip.fa").write_text(">a\nMKVAL\n>b\nMKVAL\n")
        (self.d / "in.fa").write_text(">x\nATGAAAGTTGCTCTG\n>y\nATGAAGGTCGCACTA\n")
        r = coal(["ip.fa", "in.fa"], self.d)
        self.assertEqual((r.returncode, r.stderr), (0, ""))


# ---------------------------------------------------------------- ENG-e
class TestEngE_LongFallback(unittest.TestCase):
    def test_100k_codons_with_one_n_is_fast(self):
        sys.path.insert(0, str(SRC))
        import time
        from coal.core import thread_one
        from coal.tables import get_code
        code = get_code(1).code
        sense = sorted(c for c, a in code.items() if a != "*")
        rng = random.Random(5)
        cod = ["ATG"] + [rng.choice(sense) for _ in range(99_999)]
        prot = "".join(code[c] for c in cod)
        dna = "".join(cod)
        k = 3 * 50_000 + 1
        dna = dna[:k] + "N" + dna[k + 1:]
        t = time.perf_counter()
        r = thread_one("g", prot, dna, 1)
        dt = time.perf_counter() - t
        self.assertEqual(r.dna, dna)
        self.assertEqual(r.mismatch_cols, {50_001})
        self.assertLess(dt, 3.0)  # R2-E10: was 1.0 (flaky on slow CI); ~0.1 s here, quadratic was ~10 s

    def test_many_failed_anchors_with_n_are_fast(self):  # R2-E9
        sys.path.insert(0, str(SRC))
        import time
        from coal.core import thread_one
        from coal.tables import get_code
        code = get_code(1).code
        sense = sorted(c for c, a in code.items() if a != "*")
        rng = random.Random(5)
        cod = ["ATG"] + [rng.choice(sense) for _ in range(99_999)]
        prot = "".join(code[c] for c in cod)
        for i in range(50, len(cod), 100):
            cod[i] = "NNN"
        dna = "".join(cod)
        t = time.perf_counter()
        r = thread_one("g", prot, dna, 1)
        dt = time.perf_counter() - t
        self.assertEqual(r.dna, dna)
        self.assertEqual(len(r.mismatch_cols), 1000)
        self.assertLess(dt, 2.5)  # ~0.5 s here; was ~5 s


# ---------------------------------------------------------------- ENG-f
class TestEngF_Interrupts(TempDirCase):
    """coal is stopped while it copies its finished temp file to a stdout pipe that the
    test does not drain (so it is blocked mid-copy, with the temp file present)."""

    def setUp(self):
        super().setUp()
        sys.path.insert(0, str(ROOT / "benchmarks"))
        from gen_data import generate
        generate(200, 1000, self.d / "big", seed=3)   # ~600 kB of output, > a pipe buffer
        self.tmpdir = self.d / "tmp"
        self.tmpdir.mkdir()

    def start(self):
        env = dict(os.environ, PYTHONPATH=str(SRC), TMPDIR=str(self.tmpdir))
        p = subprocess.Popen([sys.executable, "-m", "coal", "big/pep.fasta", "big/nuc.fasta", "-j", "1"],
                             cwd=self.d, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertTrue(p.stdout.readline().startswith(b">seq0"))  # now in the copy phase
        return p

    def stop(self, p, sig):
        p.send_signal(sig)
        try:
            p.stdout.close()  # only after the signal, so no broken pipe first
            err = p.stderr.read().decode()
            p.wait(timeout=TIMEOUT)
        finally:
            if p.poll() is None:
                p.kill()
                p.wait()
            p.stderr.close()  # R2-E10
        return p.returncode, err

    @unittest.skipUnless(hasattr(__import__("signal"), "SIGHUP"), "R2-E8: no SIGHUP here")
    def test_sighup_cleans_up(self):  # R2-E8
        import signal
        p = self.start()
        rc, err = self.stop(p, signal.SIGHUP)
        self.assertEqual(rc, 128 + signal.SIGHUP, err)
        self.assertNotIn("Traceback", err)
        self.assertEqual(list(self.tmpdir.iterdir()), [])

    @unittest.skipIf(WINDOWS, "R2-D1: Windows has no SIGTERM handler (send_signal = TerminateProcess)")
    def test_sigterm_cleans_up(self):
        import signal
        p = self.start()
        rc, err = self.stop(p, signal.SIGTERM)
        self.assertEqual(rc, 128 + signal.SIGTERM, err)
        self.assertNotIn("Traceback", err)
        self.assertEqual(list(self.tmpdir.iterdir()), [])

    @unittest.skipIf(WINDOWS, "R2-D1: send_signal(SIGINT) is not supported on Windows")
    def test_sigint_message_and_130(self):
        import signal
        p = self.start()
        rc, err = self.stop(p, signal.SIGINT)
        self.assertEqual(rc, 130, err)
        self.assertNotIn("Traceback", err)
        self.assertIn("coal: interrupted", err)
        self.assertEqual(list(self.tmpdir.iterdir()), [])

    def test_broken_pipe_is_quiet(self):  # coal ... | head -1
        from gen_data import generate
        generate(2000, 1000, self.d / "huge", seed=3)   # ~6 MB of output, > any pipe buffer
        env = dict(os.environ, PYTHONPATH=str(SRC), TMPDIR=str(self.tmpdir))
        p = subprocess.Popen([sys.executable, "-m", "coal", "huge/pep.fasta", "huge/nuc.fasta", "-j", "1"],
                             cwd=self.d, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            self.assertTrue(p.stdout.readline().startswith(b">seq0"))
            p.stdout.close()
            err = p.stderr.read().decode()
            p.wait(timeout=TIMEOUT)
        finally:
            if p.poll() is None:
                p.kill()
                p.wait()
            p.stderr.close()  # R2-E10
        self.assertEqual(err, "")
        self.assertEqual(p.returncode, 141)
        self.assertEqual(list(self.tmpdir.iterdir()), [])


class TestR2D1_WindowsClosedStdout(TempDirCase):  # R2-D1
    class _ClosedPipe:
        def __init__(self, err):
            self.err = err

        def write(self, s):
            raise self.err

        def flush(self):
            pass

    def run_to(self, err, windows):
        from unittest import mock
        sys.path.insert(0, str(SRC))
        from coal import pipeline
        with mock.patch.object(pipeline, "_WINDOWS", windows), mock.patch.object(sys, "stdout", self._ClosedPipe(err)):
            pipeline.run(str(self.d / "p.fa"), [str(self.d / "n.fa")])

    def test_einval_on_windows_is_a_broken_pipe(self):
        import errno
        with self.assertRaises(BrokenPipeError):
            self.run_to(OSError(errno.EINVAL, "Invalid argument"), True)

    def test_einval_elsewhere_is_a_real_error(self):
        import errno
        with self.assertRaises(OSError) as cm:
            self.run_to(OSError(errno.EINVAL, "Invalid argument"), False)
        self.assertNotIsInstance(cm.exception, BrokenPipeError)


# ---------------------------------------------------------------- ENG-g
class TestEngG_Table3Starts(TempDirCase):
    def test_gtg_start_in_table3(self):
        (self.d / "gp.fa").write_text(">g\nMKVAS\n")
        (self.d / "gn.fa").write_text(">g\nGTGAAAGTTGCTTCT\n")
        r = coal(["gp.fa", "gn.fa", "-codontable", "3"], self.d)
        self.assertEqual((r.returncode, r.stdout, r.stderr), (0, ">g\nGTGAAAGTTGCTTCT\n", ""))

    def test_table3_start_set(self):
        sys.path.insert(0, str(SRC))
        from coal.tables import get_code
        self.assertEqual(set(get_code(3).starts), {"ATA", "ATG", "GTG"})


# ---------------------------------------------------------------- ENG-h
class TestEngH_InputEncodings(TempDirCase):
    def check(self, pep: bytes, nuc: bytes):
        (self.d / "hp.fa").write_bytes(pep)
        (self.d / "hn.fa").write_bytes(nuc)
        r = coal(["hp.fa", "hn.fa"], self.d)
        self.assertEqual((r.returncode, r.stdout), (0, EXPECTED), r.stderr)

    def test_utf8_bom(self):
        bom = b"\xef\xbb\xbf"
        self.check(bom + PEP.encode(), NUC.encode())
        self.check(PEP.encode(), bom + NUC.encode())

    def test_cr_only_line_endings(self):
        self.check(PEP.replace("\n", "\r").encode(), NUC.encode())
        self.check(PEP.encode(), NUC.replace("\n", "\r").encode())
        self.check(PEP.replace("\n", "\r\n").encode(), NUC.replace("\n", "\r\n").encode())

    def test_hash_line_before_first_dna_header_is_skipped(self):
        self.check(PEP.encode(), b"# made by hand\n" + NUC.encode())

    def test_non_ascii_error_names_the_file(self):
        (self.d / "bad.fa").write_bytes(PEP.encode().replace(b"MKV-AL", b"MKV\xe9AL"))
        r = coal(["bad.fa", "n.fa"], self.d)
        self.assertEqual(r.returncode, 1)
        self.assertIn("bad.fa", r.stderr)
        self.assertNotIn("Traceback", r.stderr)


# ---------------------------------------------------------------- ENG-i2
class TestEngI_NomismatchWarnings(TempDirCase):
    def test_no_mismatch_warnings_with_nomismatch(self):
        fx = ROOT / "tests" / "fixtures" / "long_mismatch"
        args = [str(fx / "pep.fasta"), str(fx / "nuc.fasta")]
        self.assertIn("does not correspond", coal(args, self.d).stderr)
        r = coal(args + ["-nomismatch"], self.d)
        gold = (ROOT / "tests" / "golden" / "long_mismatch" / "fasta_nomismatch.stdout").read_text()
        self.assertEqual((r.returncode, r.stdout, r.stderr), (0, gold, ""))


if __name__ == "__main__":
    unittest.main()
