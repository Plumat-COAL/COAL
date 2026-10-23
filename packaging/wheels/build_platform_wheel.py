#!/usr/bin/env python3
"""Make a platform wheel of coal-align that carries the C++ engine.

    python packaging/wheels/build_platform_wheel.py --binary build/coal-cpp \
        --wheel dist/coal_align-0.4.2-py3-none-any.whl --platform-tag manylinux_2_28_x86_64 --outdir dist

Takes the pure-Python wheel (`python -m build --wheel`) and a built coal-cpp, copies the binary into the wheel as
`coal/bin/coal-cpp` (`coal-cpp.exe` for a win_* tag) with the executable bit set in the zip entry (pip, installer
and uv restore it), rewrites RECORD, marks the wheel as platform-specific (Root-Is-Purelib: false) and retags it
with `python -m wheel tags --python-tag py3 --abi-tag none --platform-tag <tag>`. The package code is unchanged:
`coal.dispatch.find_cpp()` looks in coal/bin after COAL_CPP, so `pip install coal-align` on a platform with such
a wheel gives both engines, and the pure wheel stays the fallback everywhere else. release.yml runs this once per
platform (Linux x86_64 and aarch64 in manylinux_2_28 containers, macOS, Windows); tests/test_wheels.py runs it on
the local build. No cibuildwheel, no compiled Python extension: the binary is a data file the Python engine execs.

Prints the path of the wheel it made. Exit 1 with the reason on stderr for a bad input.
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

PLATFORM_TAGS = ("manylinux_2_28_x86_64", "manylinux_2_28_aarch64", "macosx_11_0_arm64", "macosx_10_15_x86_64",
                 "macosx_11_0_universal2", "win_amd64")  # what release.yml builds; others are allowed (--force)


def _record_line(name: str, data: bytes) -> str:
    digest = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
    return f"{name},sha256={digest},{len(data)}"


def inject(wheel: str, binary: str, platform_tag: str, out_wheel: str) -> str:
    """Write out_wheel: `wheel` plus the binary at coal/bin/, RECORD and WHEEL updated. Returns the package's
    dist-info folder name (coal_align-<v>.dist-info)."""
    exe = "coal-cpp.exe" if platform_tag.startswith("win") else "coal-cpp"
    dest = f"coal/bin/{exe}"
    with open(binary, "rb") as fh:
        blob = fh.read()
    with zipfile.ZipFile(wheel) as zin:
        names = zin.namelist()
        dist_info = sorted({n.split("/")[0] for n in names if n.endswith(".dist-info/RECORD")})
        if len(dist_info) != 1:
            raise SystemExit(f"{wheel}: expected one .dist-info/RECORD, found {dist_info}")
        info = dist_info[0]
        if not info.startswith("coal_align-"):
            raise SystemExit(f"{wheel}: not a coal-align wheel ({info})")
        if dest in names:
            raise SystemExit(f"{wheel} already carries {dest}")
        if not any(n.startswith("coal/") for n in names):
            raise SystemExit(f"{wheel}: no coal/ package inside")
        record_name = f"{info}/RECORD"
        wheel_name = f"{info}/WHEEL"
        record = zin.read(record_name).decode("utf-8")
        wheel_meta = zin.read(wheel_name).decode("utf-8")
        with zipfile.ZipFile(out_wheel, "w", compression=zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                if item.filename in (record_name, wheel_name):
                    continue
                zout.writestr(item, zin.read(item.filename))
            # the binary, executable in the zip (0755, regular file) so installers restore the bit
            bi = zipfile.ZipInfo(dest, date_time=(1980, 1, 1, 0, 0, 0))
            bi.compress_type = zipfile.ZIP_DEFLATED
            bi.external_attr = (0o100755 & 0xFFFF) << 16
            bi.create_system = 3  # Unix: the mode in external_attr is honoured
            zout.writestr(bi, blob)
            # platform-specific wheels install to platlib (PEP 427); `wheel tags` sets the Tag line below
            lines = [ln for ln in wheel_meta.splitlines() if ln.strip() and not ln.startswith("Root-Is-Purelib:")]
            lines.append("Root-Is-Purelib: false")
            new_wheel_meta = "\n".join(lines) + "\n"
            wi = zipfile.ZipInfo(wheel_name, date_time=(1980, 1, 1, 0, 0, 0))
            wi.compress_type = zipfile.ZIP_DEFLATED
            wi.external_attr = (0o100644 & 0xFFFF) << 16
            wi.create_system = 3
            zout.writestr(wi, new_wheel_meta.encode("utf-8"))
            rec_lines = [ln for ln in record.splitlines()
                         if ln and not ln.startswith((wheel_name + ",", record_name + ","))]
            rec_lines.append(_record_line(dest, blob))
            rec_lines.append(_record_line(wheel_name, new_wheel_meta.encode("utf-8")))
            rec_lines.append(f"{record_name},,")
            ri = zipfile.ZipInfo(record_name, date_time=(1980, 1, 1, 0, 0, 0))
            ri.compress_type = zipfile.ZIP_DEFLATED
            ri.external_attr = (0o100644 & 0xFFFF) << 16
            ri.create_system = 3
            zout.writestr(ri, ("\n".join(rec_lines) + "\n").encode("utf-8"))
    return info


def build(binary: str, wheel: str, platform_tag: str, outdir: str) -> str:
    if not os.path.isfile(binary):
        raise SystemExit(f"no such binary: {binary}")
    if not os.path.isfile(wheel):
        raise SystemExit(f"no such wheel: {wheel}")
    base = os.path.basename(wheel)
    if not base.endswith("-py3-none-any.whl"):
        raise SystemExit(f"{base}: expected the pure wheel (…-py3-none-any.whl) as the input")
    os.makedirs(outdir, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        staged = os.path.join(tmp, base)
        inject(wheel, binary, platform_tag, staged)
        # `wheel tags` reads the wheel back (verifying RECORD), renames it after the new tags and rewrites WHEEL's Tag
        r = subprocess.run([sys.executable, "-m", "wheel", "tags", "--python-tag", "py3", "--abi-tag", "none",
                            "--platform-tag", platform_tag, "--remove", staged],
                           capture_output=True, text=True, cwd=tmp)
        if r.returncode != 0:
            raise SystemExit(f"wheel tags failed:\n{r.stdout}{r.stderr}")
        made = [f for f in os.listdir(tmp) if f.endswith(".whl")]
        if len(made) != 1:
            raise SystemExit(f"wheel tags left {made} in {tmp}")
        want = base.replace("-py3-none-any.whl", f"-py3-none-{platform_tag}.whl")
        if made[0] != want:
            raise SystemExit(f"wheel tags made {made[0]}, expected {want}")
        out = os.path.join(outdir, want)
        shutil.move(os.path.join(tmp, made[0]), out)
    return out


def _glibc(tag: str) -> tuple[int, int] | None:
    m = re.match(r"manylinux_(\d+)_(\d+)_", tag)
    return (int(m.group(1)), int(m.group(2))) if m else None


def manylinux_ok(audit_text: str, platform_tag: str) -> tuple[bool, str]:
    """Read `auditwheel show` output: is the wheel really installable under platform_tag? auditwheel names the
    lowest manylinux tag the binary's glibc symbols allow; the wheel's tag must not be lower than that (anything
    dynamic besides glibc, or a newer glibc than the manylinux_2_28 image's, raises the tag and fails here)."""
    want = _glibc(platform_tag)
    if want is None:
        return True, "not a manylinux tag; nothing to check"
    # auditwheel wraps its text at ~70 columns, so any of these spaces may be a line break
    m = re.search(r'consistent\s+with\s+the\s+following\s+platform\s+tag:\s+"(manylinux_\d+_\d+_\w+)"', audit_text)
    if not m:
        return False, "auditwheel show named no manylinux tag:\n" + audit_text
    got = _glibc(m.group(1))
    if got > want:
        return False, f"auditwheel: the binary needs {m.group(1)}, so it cannot be tagged {platform_tag}"
    return True, f"auditwheel: {m.group(1)} <= {platform_tag}"


def audit(wheel: str, platform_tag: str) -> str:
    r = subprocess.run([sys.executable, "-m", "auditwheel", "show", wheel], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"auditwheel show failed:\n{r.stdout}{r.stderr}")
    ok, why = manylinux_ok(r.stdout, platform_tag)
    if not ok:
        raise SystemExit(why)
    return why


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], prog="build_platform_wheel.py")
    ap.add_argument("--binary", required=True, help="the built coal-cpp (coal-cpp.exe on Windows)")
    ap.add_argument("--wheel", required=True, help="the pure wheel, dist/coal_align-<v>-py3-none-any.whl")
    ap.add_argument("--platform-tag", required=True, help="e.g. " + ", ".join(PLATFORM_TAGS))
    ap.add_argument("--outdir", default="dist", help="where the platform wheel goes (default: dist)")
    ap.add_argument("--force", action="store_true", help="allow a platform tag outside the known list")
    ap.add_argument("--audit", action="store_true",
                    help="manylinux tags: run `auditwheel show` on the result and fail if the tag is too low for it")
    a = ap.parse_args(argv)
    if a.platform_tag not in PLATFORM_TAGS and not a.force:
        print(f"{a.platform_tag}: not one of {', '.join(PLATFORM_TAGS)} (use --force for another)", file=sys.stderr)
        return 1
    try:
        out = build(a.binary, a.wheel, a.platform_tag, a.outdir)
        if a.audit:
            print(audit(out, a.platform_tag), file=sys.stderr)
        print(out)
    except SystemExit as e:
        if e.code and not isinstance(e.code, int):
            print(f"build_platform_wheel: {e.code}", file=sys.stderr)
            return 1
        raise
    return 0


if __name__ == "__main__":
    sys.exit(main())
