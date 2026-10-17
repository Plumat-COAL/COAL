"""Check a release tag against the version in the source and print the version (used by .github/workflows/release.yml).

    python scripts/release_version.py v0.4.0      -> prints 0.4.0, exit 0
    python scripts/release_version.py v0.5.0      -> error (source says 0.4.0), exit 1

The tag must be "v" + the version in both src/coal/__init__.py and pyproject.toml, so a release can never ship
files whose name, metadata and `coal --version` disagree. Prints `prerelease=true|false` too with --github-output.
"""
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TAG_RE = re.compile(r"^v(\d+(?:\.\d+){1,2}(?:(?:a|b|rc)\d+)?)$")


def source_versions(root=ROOT):
    with open(os.path.join(root, "src", "coal", "__init__.py"), encoding="utf-8") as fh:
        init = re.search(r'__version__\s*=\s*"([^"]+)"', fh.read()).group(1)
    with open(os.path.join(root, "pyproject.toml"), encoding="utf-8") as fh:
        proj = re.search(r'(?m)^version\s*=\s*"([^"]+)"', fh.read()).group(1)
    return init, proj


def check(tag, root=ROOT):
    """Return (version, prerelease) for a valid tag, else raise ValueError with a readable reason."""
    m = TAG_RE.match(tag or "")
    if not m:
        raise ValueError(f"tag {tag!r} is not v<major>.<minor>[.<patch>][a|b|rc<N>], e.g. v0.3.0")
    version = m.group(1)
    init, proj = source_versions(root)
    if not (version == init == proj):
        raise ValueError(f"tag {tag} says {version}, but src/coal/__init__.py says {init} and pyproject.toml says "
                         f"{proj}: bump both, commit, then tag that commit")
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
