"""Check a release tag against the version in the source and print the version (used by .github/workflows/release.yml).

    python scripts/release_version.py v0.4.0      -> prints 0.4.0, exit 0
    python scripts/release_version.py v0.5.0      -> error (source says 0.4.0), exit 1

The tag must be "v" + the one version of both packages: src/coal/__init__.py and pyproject.toml (the engine,
coal-align), apps/coal-app/pyproject.toml (the app, coal-align-app) together with the coal-align version that the app
pins (coal-align==<version>), and the C++ engine's VERSION in cpp/src/main.cpp and cpp/src/report.cpp (`coal`
hands work to a `coal-cpp` of its own version only, so a stale one would silently never be used). So a release can
never ship files whose names, metadata, dependency, `coal --version` and `coal-cpp --version` disagree. Prints `prerelease=true|false` too with --github-output.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAG_RE = re.compile(r"^v(\d+(?:\.\d+){1,2}(?:(?:a|b|rc)\d+)?)$")


APP_PYPROJECT = os.path.join("apps", "coal-app", "pyproject.toml")


def source_versions(root=ROOT):
    """(src/coal/__init__.py, pyproject.toml): the engine's version in its two places."""
    with open(os.path.join(root, "src", "coal", "__init__.py"), encoding="utf-8") as fh:
        init = re.search(r'__version__\s*=\s*"([^"]+)"', fh.read()).group(1)
    with open(os.path.join(root, "pyproject.toml"), encoding="utf-8") as fh:
        proj = re.search(r'(?m)^version\s*=\s*"([^"]+)"', fh.read()).group(1)
    return init, proj


def app_versions(root=ROOT):
    """(the app's version, the coal-align version it pins) from apps/coal-app/pyproject.toml."""
    with open(os.path.join(root, APP_PYPROJECT), encoding="utf-8") as fh:
        text = fh.read()
    app = re.search(r'(?m)^version\s*=\s*"([^"]+)"', text).group(1)
    pin = re.search(r'"coal-align==([^"\s;]+)"', text)
    return app, (pin.group(1) if pin else None)


CPP_SOURCES = (os.path.join("cpp", "src", "main.cpp"), os.path.join("cpp", "src", "report.cpp"))


def cpp_versions(root=ROOT):
    """The C++ engine's VERSION string in each of CPP_SOURCES (what `coal-cpp --version` and its reports say)."""
    out = []
    for rel in CPP_SOURCES:
        with open(os.path.join(root, rel), encoding="utf-8") as fh:
            m = re.search(r'static const char\* VERSION = "([^"]+)";', fh.read())
        out.append(m.group(1) if m else None)
    return tuple(out)


def check(tag, root=ROOT):
    """Return (version, prerelease) for a valid tag, else raise ValueError with a readable reason."""
    m = TAG_RE.match(tag or "")
    if not m:
        raise ValueError(f"tag {tag!r} is not v<major>.<minor>[.<patch>][a|b|rc<N>], e.g. v0.3.0")
    version = m.group(1)
    init, proj = source_versions(root)
    app, pin = app_versions(root)
    cpp = cpp_versions(root)
    if not (version == init == proj == app == pin) or any(c != version for c in cpp):
        cpp_says = ", ".join(f"{rel.replace(os.sep, '/')} says {c}" for rel, c in zip(CPP_SOURCES, cpp))
        raise ValueError(f"tag {tag} says {version}, but src/coal/__init__.py says {init}, pyproject.toml says "
                         f"{proj}, {APP_PYPROJECT.replace(os.sep, '/')} says {app} and pins coal-align=={pin}, "
                         f"{cpp_says}: bump all six, commit, then tag that commit")
    return version, bool(re.search(r"(a|b|rc)\d+$", version))


def main(argv):
    args = [a for a in argv if a != "--github-output"]
    if len(args) != 1:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    try:
        version, pre = check(args[0])
    except ValueError as e:
        print(f"release_version: {e}", file=sys.stderr)
        return 1
    print(version)
    if "--github-output" in argv and os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as fh:
            fh.write(f"version={version}\nprerelease={'true' if pre else 'false'}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
