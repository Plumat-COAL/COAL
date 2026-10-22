"""packaging/orthologr/pal2nal.pl: the Perl stand-in for orthologr's bundled pal2nal.pl runs COAL with the
unchanged PAL2NAL command line (`<aln> <nuc> -output fasta`), prefers coal-cpp, and gives the same bytes as
calling COAL directly. Needs perl; skipped otherwise."""
import os
import shutil
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STANDIN = ROOT / "packaging" / "orthologr" / "pal2nal.pl"
EX = ROOT / "examples" / "clean"
PERL = shutil.which("perl")


def _bin_dir():
    return str(Path(sys.executable).parent)


@unittest.skipUnless(PERL, "perl is not installed")
class TestStandIn(unittest.TestCase):
    def run_standin(self, args, env_extra=None):
        env = dict(os.environ, PATH=_bin_dir() + os.pathsep + os.environ.get("PATH", ""))
        env.pop("COAL_PAL2NAL", None)
        env.update(env_extra or {})
        return subprocess.run([PERL, str(STANDIN), *args], capture_output=True, cwd=EX, env=env, timeout=120,
                              stdin=subprocess.DEVNULL)

    def test_same_output_as_coal_pal2nal(self):
        want = subprocess.run([sys.executable, "-m", "coal_compat", "protein.aln.fasta", "dna.fasta", "-output",
                               "fasta"], capture_output=True, cwd=EX, timeout=120, stdin=subprocess.DEVNULL)
        self.assertEqual(want.returncode, 0, want.stderr)
        pal2nal = Path(_bin_dir()) / ("pal2nal.exe" if os.name == "nt" else "pal2nal")
        if not pal2nal.exists():
            self.skipTest("the pal2nal console script is not installed next to this Python")
        got = self.run_standin(["protein.aln.fasta", "dna.fasta", "-output", "fasta"],
                               {"COAL_PAL2NAL": str(pal2nal)})
        self.assertEqual(got.returncode, 0, got.stderr)
        self.assertEqual(got.stdout, want.stdout)

    def test_prefers_coal_cpp_when_available(self):
        cpp = os.environ.get("COAL_CPP_BIN") or shutil.which("coal-cpp")
        if not cpp or not Path(cpp).exists():
            self.skipTest("coal-cpp not built (set COAL_CPP_BIN)")
        env = {"PATH": str(Path(cpp).parent) + os.pathsep + _bin_dir() + os.pathsep + os.environ.get("PATH", "")}
        got = self.run_standin(["--version"], env)
        self.assertEqual(got.returncode, 0, got.stderr)
        self.assertIn(b"(C++)", got.stdout)  # coal-cpp --version: "coal <v> (C++)"
        aln = self.run_standin(["protein.aln.fasta", "dna.fasta", "-output", "fasta"], env)
        direct = subprocess.run([cpp, "protein.aln.fasta", "dna.fasta", "-output", "fasta"], capture_output=True,
                                cwd=EX, timeout=120, stdin=subprocess.DEVNULL)
        self.assertEqual((aln.returncode, aln.stdout), (0, direct.stdout))

    def test_missing_coal_is_a_clear_error(self):
        got = self.run_standin(["a", "b", "-output", "fasta"], {"COAL_PAL2NAL": str(EX / "no-such-program")})
        self.assertEqual(got.returncode, 1)
        self.assertIn(b"cannot run", got.stderr)
        self.assertIn(b"pip install coal-align", got.stderr)
        self.assertNotIn(b"Can't exec", got.stderr)  # Perl's own warning is silenced


if __name__ == "__main__":
    unittest.main()
