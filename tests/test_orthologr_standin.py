"""packaging/orthologr/pal2nal.pl: the Perl stand-in for orthologr's bundled pal2nal.pl runs COAL with the
unchanged PAL2NAL command line (`<aln> <nuc> -output fasta`), prefers coal-cpp for the options it understands,
hands PAL2NAL-only options (-nostderr, -h) to the `pal2nal` drop-in, and gives the same bytes as calling COAL
directly. Needs perl; skipped otherwise."""
import os
import shutil
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests"))
from subprocess_env import absolute_pythonpath, python_env  # noqa: E402

STANDIN = ROOT / "packaging" / "orthologr" / "pal2nal.pl"
EX = ROOT / "examples" / "clean"
PERL = shutil.which("perl")


def _bin_dir():
    """Where pip put the console scripts: next to python on Linux and macOS, Scripts\\ on Windows (sysconfig)."""
    import sysconfig
    name = "pal2nal.exe" if os.name == "nt" else "pal2nal"
    for d in (sysconfig.get_path("scripts"), str(Path(sys.executable).parent)):
        if d and (Path(d) / name).exists():
            return d
    return str(Path(sys.executable).parent)


@unittest.skipUnless(PERL, "perl is not installed")
class TestStandIn(unittest.TestCase):
    def run_standin(self, args, env_extra=None):
        env = python_env({"PATH": _bin_dir() + os.pathsep + os.environ.get("PATH", "")})   # runs in examples/clean
        env.pop("COAL_PAL2NAL", None)
        env.update(env_extra or {})
        return subprocess.run([PERL, str(STANDIN), *args], capture_output=True, cwd=EX, env=env, timeout=120,
                              stdin=subprocess.DEVNULL)

    def test_same_output_as_coal_pal2nal(self):
        want = subprocess.run([sys.executable, "-m", "coal_compat", "protein.aln.fasta", "dna.fasta", "-output",
                               "fasta"], capture_output=True, cwd=EX, timeout=120, stdin=subprocess.DEVNULL,
                              env=python_env())
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

    def _cpp_env(self):
        cpp = os.environ.get("COAL_CPP_BIN") or shutil.which("coal-cpp")
        if not cpp or not Path(cpp).exists():
            self.skipTest("coal-cpp not built (set COAL_CPP_BIN)")
        pal2nal = Path(_bin_dir()) / ("pal2nal.exe" if os.name == "nt" else "pal2nal")
        if not pal2nal.exists():
            self.skipTest("the pal2nal console script is not installed next to this Python")
        return cpp, {"PATH": str(Path(cpp).parent) + os.pathsep + _bin_dir() + os.pathsep
                     + os.environ.get("PATH", "")}

    def test_pal2nal_only_options_go_to_the_drop_in(self):
        # coal-cpp rejects -nostderr and prints coal's help for -h; with coal-cpp first on PATH these must still
        # behave as PAL2NAL's (the drop-in accepts them)
        cpp, env = self._cpp_env()
        want = subprocess.run([sys.executable, "-m", "coal_compat", "protein.aln.fasta", "dna.fasta", "-output",
                               "fasta", "-nostderr"], capture_output=True, cwd=EX, timeout=120,
                              stdin=subprocess.DEVNULL)
        self.assertEqual(want.returncode, 0, want.stderr)
        got = self.run_standin(["protein.aln.fasta", "dna.fasta", "-output", "fasta", "-nostderr"], env)
        self.assertEqual(got.returncode, 0, got.stderr)
        self.assertEqual(got.stdout, want.stdout)
        self.assertNotIn(b"unrecognized arguments", got.stderr)
        for args in (["-h"], ["-help"], []):
            with self.subTest(args=args):
                got = self.run_standin(args, env)
                direct = subprocess.run([sys.executable, "-m", "coal_compat", *args], capture_output=True, cwd=EX,
                                        timeout=120, stdin=subprocess.DEVNULL)
                self.assertEqual(got.returncode, direct.returncode, got.stderr)
                self.assertEqual(got.stdout + got.stderr, direct.stdout + direct.stderr)  # PAL2NAL-style help
                self.assertNotIn(b"usage: coal [-h]", got.stdout + got.stderr)

    def test_coal_cpp_still_takes_the_orthologr_command_lines(self):
        cpp, env = self._cpp_env()
        for args in (["protein.aln.fasta", "dna.fasta", "-output", "fasta", "-codontable", "1", "-nogap",
                      "-nomismatch", "-blockonly"],
                     ["protein.aln.fasta", "dna.fasta", "-output", "paml"]):
            with self.subTest(args=args):
                got = self.run_standin(args, env)
                direct = subprocess.run([cpp, *args], capture_output=True, cwd=EX, timeout=120,
                                        stdin=subprocess.DEVNULL)
                self.assertEqual((got.returncode, got.stdout), (direct.returncode, direct.stdout), got.stderr)
        # COAL_PAL2NAL stays the override, also for a PAL2NAL-only option the chosen program may not take
        got = self.run_standin(["-h"], {**env, "COAL_PAL2NAL": cpp})
        self.assertEqual(got.returncode, 0, got.stderr)
        self.assertIn(b"usage: coal [-h]", got.stdout)  # coal-cpp's help, as asked

    def test_missing_coal_is_a_clear_error(self):
        got = self.run_standin(["a", "b", "-output", "fasta"], {"COAL_PAL2NAL": str(EX / "no-such-program")})
        self.assertEqual(got.returncode, 1)
        self.assertIn(b"cannot run", got.stderr)
        self.assertIn(b"pip install coal-align", got.stderr)
        self.assertNotIn(b"Can't exec", got.stderr)  # Perl's own warning is silenced


class TestSubprocessEnv(unittest.TestCase):
    """tests/subprocess_env.py: a relative PYTHONPATH (PYTHONPATH=src) is made absolute, src/ first, so a
    subprocess in another folder imports this checkout (REVIEW-COAL-3 L4)."""

    def test_relative_entries_become_absolute(self):
        sep = os.pathsep
        self.assertEqual(absolute_pythonpath(None), str(ROOT / "src"))
        self.assertEqual(absolute_pythonpath("src", base=str(ROOT)), str(ROOT / "src"))   # no second copy
        got = absolute_pythonpath(sep.join(["", "lib", str(ROOT / "x")]), base=str(ROOT))
        self.assertEqual(got.split(sep), [str(ROOT / "src"), str(ROOT / "lib"), str(ROOT / "x")])
        env = python_env({"A": "1"}, base={"PYTHONPATH": "src", "PATH": "/bin"})
        self.assertEqual(env["A"], "1")
        self.assertEqual(env["PATH"], "/bin")
        self.assertTrue(all(os.path.isabs(p) for p in env["PYTHONPATH"].split(sep)), env["PYTHONPATH"])

    def test_a_subprocess_in_another_folder_imports_this_checkout(self):
        r = subprocess.run([sys.executable, "-c", "import coal, coal_compat; print(coal.__file__)"],
                           capture_output=True, text=True, cwd=EX, timeout=60, stdin=subprocess.DEVNULL,
                           env=python_env(base={**os.environ, "PYTHONPATH": "src"}))
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(Path(r.stdout.strip()).resolve().is_relative_to((ROOT / "src").resolve()), r.stdout)


if __name__ == "__main__":
    unittest.main()
