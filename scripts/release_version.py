"""Check a release tag against the version in the source and print the version (used by .github/workflows/release.yml).

    python scripts/release_version.py v0.4.0      -> prints 0.4.0, exit 0
    python scripts/release_version.py v0.5.0      -> error (source says 0.4.0), exit 1
    python scripts/release_version.py --dry-run main   -> prints the source's version (e.g. 0.6.1.dev0), exit 0 when
                                                     the six places agree; for a build that publishes nothing
    ... --github-output   also appends version=, prerelease=true|false and r_version= to $GITHUB_OUTPUT

A release (not a dry run) is also refused while CHANGELOG.md has entries under "## Unreleased" (the release commit
moves them into "## X.Y.Z"), while coalr's DESCRIPTION does not carry r_version(X.Y.Z) or while coalr's NEWS.md
does not start with "# coalr X.Y.Z". r_version is empty for a
prerelease (a1, b1, rc1: not a valid R package version), whose release has no coalr asset; a dry run is reported
as a prerelease.

The tag must be "v" + the one version of both packages: src/coal/__init__.py and pyproject.toml (the engine,
coal-align), apps/coal-app/pyproject.toml (the app, coal-align-app) together with the coal-align version that the app
pins (coal-align==<version>), and the C++ engine's VERSION in cpp/src/main.cpp and cpp/src/report.cpp (`coal`
hands work to a `coal-cpp` of its own version only, so a stale one would silently never be used). So a release can
never ship files whose names, metadata, dependency, `coal --version` and `coal-cpp --version` disagree. Prints `prerelease=true|false` too with --github-output.

Between releases the source carries a development version, X.Y.Z.devN (PEP 440: 0.6.0.dev0 is the tree that becomes
0.6.0), in those six places, so that no build of the development tree reports the version of a published release.
A development version cannot be tagged (TAG_RE has no .dev): the release commit sets X.Y.Z everywhere first. The
places that take digits only or describe the release carry the base version X.Y.Z (base_version): CMake's
project(VERSION), CITATION.cff, the website's version spans; coalr's DESCRIPTION has R's form of it
(r_version: 0.6.0.dev0 -> 0.6.0.9000). tests/test_feature_combinations.py checks every place.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAG_RE = re.compile(r"^v(\d+(?:\.\d+){1,2}(?:(?:a|b|rc)\d+)?)$")


APP_PYPROJECT = os.path.join("apps", "coal-app", "pyproject.toml")
DEV_RE = re.compile(r"^(\d+(?:\.\d+){1,2})\.dev(\d+)$")


def base_version(version):
    """The release a version belongs to, digits only: 0.6.0.dev0 -> 0.6.0, 0.5.0rc1 -> 0.5.0, 0.5.0 -> 0.5.0."""
    m = DEV_RE.match(version)
    return m.group(1) if m else re.sub(r"(a|b|rc)\d+$", "", version)


PRE_RE = re.compile(r"(a|b|rc)\d+$")


def r_version(version):
    """The version of the R package coalr for a COAL version: the same, or for a development version R's
    convention for one (a fourth component from 9000): 0.6.0.dev0 -> 0.6.0.9000, 0.6.0.dev3 -> 0.6.0.9003.
    None for a prerelease (0.6.1rc1): R has no form of it that sorts before 0.6.1, so a prerelease ships no coalr
    (release.yml skips the coalr job) and its commit leaves DESCRIPTION at the development form."""
    m = DEV_RE.match(version)
    if m:
        return f"{m.group(1)}.{9000 + int(m.group(2))}"
    return None if PRE_RE.search(version) else version


DESCRIPTION = os.path.join("packaging", "r", "coalr", "DESCRIPTION")


def description_version(root=ROOT):
    """coalr's Version: in packaging/r/coalr/DESCRIPTION (None when the file or the field is missing)."""
    try:
        with open(os.path.join(root, DESCRIPTION), encoding="utf-8") as fh:
            m = re.search(r"(?m)^Version:\s*(\S+)\s*$", fh.read())
    except OSError:
        return None
    return m.group(1) if m else None


NEWS = os.path.join("packaging", "r", "coalr", "NEWS.md")


def news_heading(root=ROOT):
    """The first "# " heading of coalr's NEWS.md ("" when the file is missing or has none)."""
    try:
        with open(os.path.join(root, NEWS), encoding="utf-8") as fh:
            m = re.search(r"(?m)^# (.+?)\s*$", fh.read())
    except OSError:
        return ""
    return m.group(1) if m else ""


def unreleased_text(root=ROOT):
    """The text under CHANGELOG.md's "## Unreleased" heading up to the next "## " heading, stripped ("" when the
    heading is missing or the section is empty). A release has none: the release commit renames the section."""
    try:
        with open(os.path.join(root, "CHANGELOG.md"), encoding="utf-8") as fh:
            text = fh.read()
    except OSError:
        return ""
    m = re.search(r"(?m)^## Unreleased[ \t]*\n(.*?)(?=^## |\Z)", text, re.S)
    return m.group(1).strip() if m else ""


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
        if DEV_RE.match(init or ""):
            raise ValueError(f"the source is the development version {init}: a release is tagged from a commit that "
                             f"sets {base_version(init)} (tag {tag} says {version}) in all six places (bump all six: "
                             f"src/coal/__init__.py, pyproject.toml, {APP_PYPROJECT.replace(os.sep, '/')} twice, "
                             + ", ".join(rel.replace(os.sep, "/") for rel in CPP_SOURCES) + ")")
        cpp_says = ", ".join(f"{rel.replace(os.sep, '/')} says {c}" for rel, c in zip(CPP_SOURCES, cpp))
        raise ValueError(f"tag {tag} says {version}, but src/coal/__init__.py says {init}, pyproject.toml says "
                         f"{proj}, {APP_PYPROJECT.replace(os.sep, '/')} says {app} and pins coal-align=={pin}, "
                         f"{cpp_says}: bump all six, commit, then tag that commit")
    pre = bool(PRE_RE.search(version))
    if unreleased_text(root):
        raise ValueError(f"CHANGELOG.md has entries under \"## Unreleased\": the release commit moves them into "
                         f"\"## {version}\" (docs/development.md \"Releasing COAL\")")
    want_r = r_version(version)
    if want_r is not None and description_version(root) != want_r:
        raise ValueError(f"{DESCRIPTION.replace(os.sep, '/')} says Version: {description_version(root)}, not {want_r} "
                         f"(r_version of {version}): set it, run scripts/sync_r_sources.py, commit, then tag")
    if want_r is not None and news_heading(root) != f"coalr {want_r}":
        raise ValueError(f"{NEWS.replace(os.sep, '/')} starts with \"# {news_heading(root)}\", not \"# coalr {want_r}\": "
                         f"the release's coalr tarball ships it")
    return version, pre


def dry_run(root=ROOT):
    """For a build that publishes nothing (release.yml dispatched with publish off on a branch or commit, e.g. main
    between releases): the source's version, which may be a development version, once the six places agree. Raises
    ValueError when they do not."""
    init, proj = source_versions(root)
    app, pin = app_versions(root)
    cpp = cpp_versions(root)
    if len({init, proj, app, pin, *cpp}) != 1:
        raise ValueError(f"dry run: the versions disagree: {init}, {proj}, app {app} pinning {pin}, coal-cpp {cpp}")
    return init


def main(argv):
    flags = {"--github-output", "--dry-run"}
    args = [a for a in argv if a not in flags]
    if len(args) != 1:
        print(__doc__.strip(), file=sys.stderr)
        return 2
    try:
        if "--dry-run" in argv:
            version, pre = dry_run(), True
            print(f"dry run of {args[0]}: version {version}", file=sys.stderr)
        else:
            version, pre = check(args[0])
    except ValueError as e:
        print(f"release_version: {e}", file=sys.stderr)
        return 1
    print(version)
    if "--github-output" in argv and os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a", encoding="utf-8") as fh:
            fh.write(f"version={version}\nprerelease={'true' if pre else 'false'}\nr_version={r_version(version) or ''}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
