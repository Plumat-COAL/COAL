"""The platform wheels that carry coal-cpp (packaging/wheels/build_platform_wheel.py, issue #16).

Builds the pure wheel (`python -m build --wheel`), turns it into a manylinux_2_28_x86_64 wheel with the local
build/coal-cpp and build/coal-cpp-align (the tag is nominal here: the binary is the host's, so the wheel is for this test's venv only),
installs it into a fresh virtual environment (uv when it is around, else `python -m venv`; pip verifies RECORD's
hashes on the way) and checks that `coal --cpp-path` is the bundled binary, that it is executable and that
`COAL_ENGINE=cpp coal examples/clean/...` is byte-identical to the Python engine, and runs
packaging/wheels/check_wheel_binaries.py there (both binaries found inside the wheel, one alignment through
coal-cpp-align, one `coal pipeline pairs` run with an empty PATH). Skipped without `build` and
`wheel` (test-only dependencies: `pip install build wheel`) or without build/coal-cpp (COAL_CPP_BIN).

Run:  python -m unittest tests.test_wheels
"""
from __future__ import annotations

import importlib.util
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "packaging" / "wheels" / "build_platform_wheel.py"
EXAMPLE = (ROOT / "examples" / "clean" / "protein.aln.fasta", ROOT / "examples" / "clean" / "dna.fasta")
TAG = "manylinux_2_28_x86_64" if sys.platform.startswith("linux") else (
    "win_amd64" if os.name == "nt" else "macosx_11_0_universal2")


def cpp_binary() -> Path | None:
    p = Path(os.environ.get("COAL_CPP_BIN", ROOT / "build" / "coal-cpp"))
    if not p.is_absolute():
        p = ROOT / p
    for cand in (p, p.with_suffix(".exe"), p.parent / "Release" / "coal-cpp.exe"):
        if cand.is_file():
            return cand
    return None


def align_binary() -> Path | None:
    """coal-cpp-align of the same build (COAL_CPP_ALIGN, else next to coal-cpp)."""
    cpp = cpp_binary()
    env = os.environ.get("COAL_CPP_ALIGN")
    cands = [Path(env)] if env else []
    if cpp is not None:
        cands.append(cpp.parent / ("coal-cpp-align" + cpp.suffix))
    for cand in cands:
        if cand.is_file():
            return cand
    return None


def _have(mod: str) -> bool:
    return importlib.util.find_spec(mod) is not None


def _load_script():
    spec = importlib.util.spec_from_file_location("build_platform_wheel", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def run(cmd, **kw) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, capture_output=True, text=True, timeout=600, **kw)


class TestAuditParse(unittest.TestCase):
    """manylinux_ok reads `auditwheel show` (no auditwheel needed for this)."""

    def test_reads_the_tag_auditwheel_names(self):
        m = _load_script()
        ok28 = 'x.whl is consistent with\nthe following platform tag: "manylinux_2_28_x86_64".\n'
        self.assertTrue(m.manylinux_ok(ok28, "manylinux_2_28_x86_64")[0])
        # auditwheel 6.4 wraps the sentence at ~70 columns, so the break can fall anywhere in it (this exact text
        # failed the v0.5.0rc2 Linux wheel jobs)
        real = ('\ncoal_align-0.5.0rc2-py3-none-manylinux_2_28_x86_64.whl is consistent\nwith the following platform tag: '
                '"manylinux_2_28_x86_64".\n\nThe wheel references external versioned symbols in these\n')
        self.assertTrue(m.manylinux_ok(real, "manylinux_2_28_x86_64")[0])
        self.assertFalse(m.manylinux_ok(real.replace("2_28_x86_64\".", "2_38_x86_64\"."), "manylinux_2_28_x86_64")[0])
        self.assertTrue(m.manylinux_ok(ok28.replace("2_28", "2_17"), "manylinux_2_28_x86_64")[0])  # lower is fine
        ok, why = m.manylinux_ok(ok28.replace("2_28", "2_39"), "manylinux_2_28_x86_64")
        self.assertFalse(ok)
        self.assertIn("needs manylinux_2_39_x86_64", why)
        self.assertFalse(m.manylinux_ok("nothing useful", "manylinux_2_28_aarch64")[0])
        self.assertTrue(m.manylinux_ok("", "win_amd64")[0])  # not a manylinux tag: nothing to check


@unittest.skipIf(cpp_binary() is None, "coal-cpp not built (set COAL_CPP_BIN)")
@unittest.skipIf(align_binary() is None, "coal-cpp-align not built (set COAL_CPP_ALIGN)")
@unittest.skipIf(not (_have("build") and _have("wheel")), "needs the build and wheel packages")
class TestPlatformWheel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="coal_wheel_"))
        dist = cls.tmp / "dist"
        r = run([sys.executable, "-m", "build", "--wheel", "--outdir", str(dist), str(ROOT)])
        if r.returncode != 0:
            raise unittest.SkipTest(f"python -m build failed here:\n{r.stdout[-2000:]}{r.stderr[-2000:]}")
        pure = [p for p in dist.iterdir() if p.name.endswith("-py3-none-any.whl")]
        assert len(pure) == 1, list(dist.iterdir())
        cls.pure = pure[0]
        r = run([sys.executable, str(SCRIPT), "--binary", str(cpp_binary()), "--align-binary", str(align_binary()),
                 "--wheel", str(cls.pure), "--platform-tag", TAG, "--outdir", str(dist)])
        assert r.returncode == 0, r.stdout + r.stderr
        cls.wheel = Path(r.stdout.strip())
        cls.exe = "coal-cpp.exe" if os.name == "nt" else "coal-cpp"
        cls.align_exe = "coal-cpp-align.exe" if os.name == "nt" else "coal-cpp-align"

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_name_and_contents(self):
        self.assertEqual(self.wheel.name, self.pure.name.replace("-py3-none-any.whl", f"-py3-none-{TAG}.whl"))
        self.assertTrue(self.wheel.is_file())
        with zipfile.ZipFile(self.wheel) as z:
            info = z.getinfo(f"coal/bin/{self.exe}")
            self.assertTrue(stat.S_ISREG(info.external_attr >> 16))
            self.assertEqual((info.external_attr >> 16) & 0o777, 0o755)  # the executable bit travels in the zip
            self.assertEqual(z.read(info), cpp_binary().read_bytes())
            ainfo = z.getinfo(f"coal/bin/{self.align_exe}")
            self.assertTrue(stat.S_ISREG(ainfo.external_attr >> 16))
            self.assertEqual((ainfo.external_attr >> 16) & 0o777, 0o755)
            self.assertEqual(z.read(ainfo), align_binary().read_bytes())
            self.assertEqual(sorted(n for n in z.namelist() if n.startswith("coal/bin/")),
                             sorted([f"coal/bin/{self.exe}", f"coal/bin/{self.align_exe}"]))
            dist_info = [n for n in z.namelist() if n.endswith(".dist-info/WHEEL")][0].rsplit("/", 1)[0]
            meta = z.read(f"{dist_info}/WHEEL").decode()
            self.assertIn(f"Tag: py3-none-{TAG}\n", meta)
            self.assertIn("Root-Is-Purelib: false\n", meta)
            self.assertNotIn("Root-Is-Purelib: true", meta)
            self.assertNotIn("\n\n", meta.strip())
            record = z.read(f"{dist_info}/RECORD").decode()
            self.assertEqual(record.count(f"coal/bin/{self.exe},sha256="), 1)
            self.assertEqual(record.count(f"coal/bin/{self.align_exe},sha256="), 1)
            self.assertEqual(record.count(f"{dist_info}/WHEEL,sha256="), 1)
            # every file in the zip is in RECORD and vice versa
            in_zip = {n for n in z.namelist() if not n.endswith("/")}
            in_record = {ln.split(",")[0] for ln in record.splitlines() if ln}
            self.assertEqual(in_zip, in_record)
        # `wheel unpack` verifies every hash in RECORD
        r = run([sys.executable, "-m", "wheel", "unpack", "--dest", str(self.tmp / "unpacked"), str(self.wheel)])
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        # the pure wheel is untouched and carries no binary
        with zipfile.ZipFile(self.pure) as z:
            self.assertFalse(any(n.startswith("coal/bin/") for n in z.namelist()))

    def test_aligner_is_optional_only_on_request(self):
        """Without --align-binary the aligner next to --binary is taken; --no-align leaves it out; a named aligner
        that is not there is an error (release.yml names it, so a build without it cannot make a wheel)."""
        out = self.tmp / "variants"
        r = run([sys.executable, str(SCRIPT), "--binary", str(cpp_binary()), "--wheel", str(self.pure),
                 "--platform-tag", TAG, "--outdir", str(out / "sibling")])
        self.assertEqual(r.returncode, 0, r.stderr)
        with zipfile.ZipFile(r.stdout.strip()) as z:
            self.assertEqual(f"coal/bin/{self.align_exe}" in z.namelist(),
                             (cpp_binary().parent / self.align_exe).is_file())
        r = run([sys.executable, str(SCRIPT), "--binary", str(cpp_binary()), "--no-align", "--wheel", str(self.pure),
                 "--platform-tag", TAG, "--outdir", str(out / "none")])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("no coal-cpp-align in this wheel", r.stderr)
        with zipfile.ZipFile(r.stdout.strip()) as z:
            self.assertEqual([n for n in z.namelist() if n.startswith("coal/bin/")], [f"coal/bin/{self.exe}"])
            record = z.read([n for n in z.namelist() if n.endswith(".dist-info/RECORD")][0]).decode()
            self.assertNotIn("coal-cpp-align", record)
        r = run([sys.executable, str(SCRIPT), "--binary", str(cpp_binary()), "--align-binary", str(self.tmp / "none"),
                 "--wheel", str(self.pure), "--platform-tag", TAG, "--outdir", str(out / "x")])
        self.assertEqual(r.returncode, 1)
        self.assertIn("no such binary", r.stderr)
        r = run([sys.executable, str(SCRIPT), "--binary", str(cpp_binary()), "--align-binary", str(align_binary()),
                 "--no-align", "--wheel", str(self.pure), "--platform-tag", TAG, "--outdir", str(out / "x")])
        self.assertEqual(r.returncode, 1)

    def test_refuses_bad_inputs(self):
        r = run([sys.executable, str(SCRIPT), "--binary", str(cpp_binary()), "--wheel", str(self.wheel),
                 "--platform-tag", TAG, "--outdir", str(self.tmp / "x")])
        self.assertEqual(r.returncode, 1)
        self.assertIn("expected the pure wheel", r.stderr)
        r = run([sys.executable, str(SCRIPT), "--binary", str(cpp_binary()), "--wheel", str(self.pure),
                 "--platform-tag", "linux_x86_64", "--outdir", str(self.tmp / "x")])
        self.assertEqual(r.returncode, 1)
        self.assertIn("not one of", r.stderr)
        r = run([sys.executable, str(SCRIPT), "--binary", str(self.tmp / "none"), "--wheel", str(self.pure),
                 "--platform-tag", TAG, "--outdir", str(self.tmp / "x")])
        self.assertEqual(r.returncode, 1)
        self.assertIn("no such binary", r.stderr)

    def test_installs_and_coal_uses_the_bundled_engine(self):
        venv = self.tmp / "venv"
        uv = shutil.which("uv") or ("/root/.local/bin/uv" if os.path.exists("/root/.local/bin/uv") else None)
        if uv:
            r = run([uv, "venv", "-q", "-p", sys.executable, str(venv)])
            self.assertEqual(r.returncode, 0, r.stderr)
            r = run([uv, "pip", "install", "-q", "-p", str(venv), str(self.wheel)])
        else:
            r = run([sys.executable, "-m", "venv", str(venv)])
            self.assertEqual(r.returncode, 0, r.stderr)
            py = venv / ("Scripts" if os.name == "nt" else "bin") / "python"
            r = run([str(py), "-m", "pip", "install", "-q", str(self.wheel)])
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)  # pip/uv verify RECORD's hashes here
        bindir = venv / ("Scripts" if os.name == "nt" else "bin")
        coal = bindir / ("coal.exe" if os.name == "nt" else "coal")
        env = {k: v for k, v in os.environ.items() if k not in ("COAL_CPP", "COAL_ENGINE", "PYTHONPATH")}
        r = run([str(coal), "--cpp-path"], env=env, cwd=str(self.tmp))
        self.assertEqual(r.returncode, 0, r.stderr)
        path = Path(r.stdout.strip())
        self.assertEqual(path.name, self.exe)
        self.assertEqual(path.parent.name, "bin")
        self.assertEqual(path.parent.parent.name, "coal")
        self.assertIn("site-packages", path.parts)
        self.assertTrue(path.is_file())
        self.assertTrue(os.access(path, os.X_OK))
        self.assertEqual(path.read_bytes(), cpp_binary().read_bytes())
        self.assertTrue(path.resolve().is_relative_to(venv.resolve()), (path, venv))  # macOS: /var -> /private/var
        with_cpp = run([str(coal), *map(str, EXAMPLE)], env={**env, "COAL_ENGINE": "cpp"}, cwd=str(self.tmp))
        self.assertEqual(with_cpp.returncode, 0, with_cpp.stderr)
        with_py = run([str(coal), *map(str, EXAMPLE)], env={**env, "COAL_ENGINE": "python"}, cwd=str(self.tmp))
        self.assertEqual(with_py.returncode, 0, with_py.stderr)
        self.assertEqual(with_cpp.stdout, with_py.stdout)
        self.assertTrue(with_cpp.stdout.startswith(">"))
        # the aligner: found inside the wheel by coal.align, runs, and the pipeline runs on the wheel alone
        py = bindir / ("python.exe" if os.name == "nt" else "python")
        r = run([str(py), "-c", "import coal.align as a; print(a.find_align_binary())"],
                env={k: v for k, v in env.items() if k != "COAL_CPP_ALIGN"}, cwd=str(self.tmp))
        self.assertEqual(r.returncode, 0, r.stderr)
        apath = Path(r.stdout.strip())
        self.assertEqual((apath.name, apath.parent.name, apath.parent.parent.name), (self.align_exe, "bin", "coal"))
        self.assertTrue(apath.resolve().is_relative_to(venv.resolve()), (apath, venv))
        self.assertEqual(apath.read_bytes(), align_binary().read_bytes())
        r = run([str(py), str(SCRIPT.parent / "check_wheel_binaries.py")], env=env, cwd=str(self.tmp))
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("coal pipeline pairs: 6 rows", r.stdout)
        self.assertIn("coal.align.align_pairs: 3 pairs by coal-cpp-align", r.stdout)
        # the bundled binary is this version, so a large input goes to it under COAL_ENGINE=auto too
        r = run([str(path), "--version"])
        self.assertEqual(r.stdout.split()[:2], ["coal", run([str(coal), "--version"]).stdout.split()[1]])


if __name__ == "__main__":
    unittest.main()
