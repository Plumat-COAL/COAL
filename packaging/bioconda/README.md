# Bioconda recipe for `coal-align`

This folder is a ready-to-submit [Bioconda](https://bioconda.github.io/) recipe. Once it is merged into
`bioconda-recipes`, pipeline users get both engines with one command:

```
conda install -c conda-forge -c bioconda coal-align
coal --version        # Python engine, the `coal` and `pal2nal` commands, `import coal`
coal-cpp --version    # C++ engine; `coal` hands inputs >= 1 MB to it automatically
```

`coal` finds `coal-cpp` without any setting: `coal.dispatch.find_cpp()` looks at `COAL_CPP`, then inside the installed
package (`coal/bin/coal-cpp`, where the PyPI platform wheels carry it; a conda build from the sdist has no `coal/bin`,
so this step finds nothing here), then next to the Python interpreter running COAL (`$CONDA_PREFIX/bin`, where
`build.sh` installs the binary), then on `PATH`. The recipe test `test_dispatch.sh` checks both the discovery and that
a 1 MB input really runs through `coal-cpp`. The recipe builds from the sdist on purpose (Bioconda builds from source,
and the conda package's own `coal-cpp` in `$PREFIX/bin` is the one `coal` should run), so it does not use the platform
wheels. `coal --cpp-path` prints which binary `coal` found.

Files:

| File | Purpose |
| --- | --- |
| `meta.yaml` | The recipe (version, sha256, requirements, tests, metadata). |
| `build.sh` | Builds `coal-cpp` with CMake (Release, `COAL_BUILD_TOOLS=OFF`) into `$PREFIX/bin`, then `pip install .` of the Python package. |
| `test_protein.aln.fasta`, `test_dna.fasta` | Tiny fixture pair (a copy of `examples/clean/`) for the recipe's tests. |
| `test_dispatch.sh` | Recipe test: `COAL_ENGINE=cpp` works with no `COAL_CPP`; a >= 1 MB input dispatches to `coal-cpp`; both engines' output is identical. |

Bioconda builds Linux (x86_64, optionally aarch64) and macOS (x86_64, arm64). It does not build Windows; the
Windows app and `coal-cpp.exe` come from [GitHub Releases](https://github.com/plumat-coal/COAL/releases).

## Which release the recipe fits

The recipe is written for the **next release, 0.6.0**: the first under GPL-3.0-only and the first that has Ka/Ks
(`coal kaks`, `coal dnds`), `coal pipeline` and the second binary `coal-cpp-align`, all of which the recipe's tests
run. That tag does not exist yet, so `meta.yaml` carries two **placeholders, marked `PLACEHOLDER`**: `version`
(`0.6.0`) and `sha256` (64 zeros, which is not a checksum; `conda build` fails on it by design). At release time set
both (below, "The sha256 step") and nothing else needs to change.

`build.sh` builds the C++ programs and copies `coal-cpp` **and** `coal-cpp-align` into `$PREFIX/bin` (nothing
else of the C++ build: see "What the package is" below), then runs `pip install .` for the Python package. The
licence files travel as `about: license_file:`.

**Dependencies.** `run:` is Python only: codon alignment, `coal kaks` and `coal dnds` need nothing else.
`coal pipeline` runs DIAMOND, and Clustal Omega or MAFFT in cluster mode; conda has no optional dependencies, so they
are listed under `run_constrained:` (Bioconda's convention for optional companions: not installed by this package,
but held to versions COAL works with when the user installs them: `diamond >=2.1.0`, the first with DeepClust;
COAL is tested with 2.2.8). The description and `extra: notes:` say `conda install -c bioconda diamond clustalo`.
If the Bioconda reviewers prefer hard dependencies, move `diamond` and `clustalo` to `run:` (both are small).

**Licence.** `about: license:` is `GPL-3.0-only` (`license_family: GPL3`, `LICENSE` and `NOTICE` as licence files),
which is right for the releases after 0.5.0. COAL 0.5.0 and earlier are MIT and have no Ka/Ks: do not submit this
recipe for one of those tags.

## Testing the recipe locally

Without conda at all, the build steps can be run by hand against a virtual environment, which is what `find_cpp`
sees in a conda environment too (binary next to the interpreter):

```
python3 -m venv .venv-conda && .venv-conda/bin/pip install 'setuptools>=77.0.3'
PREFIX=$PWD/.venv-conda PYTHON=$PWD/.venv-conda/bin/python bash packaging/bioconda/build.sh
cd packaging/bioconda && PATH=$OLDPWD/.venv-conda/bin:$PATH bash test_dispatch.sh
```

With conda (the real thing; [bioconda-utils](https://bioconda.github.io/contributor/building-locally.html) is what
Bioconda's CI runs):

```
# conda-build directly: fastest check that the recipe builds and its tests pass
conda install -c conda-forge conda-build
conda build -c conda-forge -c bioconda packaging/bioconda
conda create -n coal-test -c local -c conda-forge -c bioconda coal-align && conda activate coal-test && coal --version

# bioconda-utils: the same lint + build the PR will get (inside a bioconda-recipes checkout)
conda create -n bioconda -c conda-forge -c bioconda bioconda-utils
conda activate bioconda
bioconda-utils lint recipes config.yml --packages coal-align
bioconda-utils build recipes config.yml --packages coal-align
# or in Bioconda's Docker image: bioconda-utils build --docker --mulled-test recipes config.yml --packages coal-align
```

Status of the current recipe (2026-11-13, second review):

- **Rendered with conda-build 26's `conda-render`** (`-c conda-forge`, `source:` pointed at a copy of the tree, the
  version set to the tree's): valid YAML with the Jinja resolved, 23 test commands, `python` pinned by the build
  variant (`python >=3.12,<3.13.0a0` and `python_abi` in `run:` for the 3.12 variant), `additional-platforms`
  as written. With `--python 3.9` the recipe renders to nothing: the `py<310` selector skips it.
- **`build.sh` run against a virtual environment** as above (`PREFIX` = the environment): `bin/coal`,
  `bin/pal2nal`, `bin/coal-cpp`, `bin/coal-cpp-align` and the Python package; no `lib/libcoal_core.a`, no
  `include/coal/`, no `lib/cmake/coal/`.
- **All 23 rendered test commands run under `set -ex`** in that environment (no `COAL_CPP`, no `COAL_ENGINE`),
  `test_dispatch.sh` included: all pass. The line that failed before (`pal2nal > usage.txt && grep ...`: `pal2nal`
  without arguments exits 1 with the usage on stderr) now tests that behaviour.
- **Not done: a real `conda build` and `bioconda-utils lint`.** The build was started and stopped while it
  unpacked the compiler packages: the machine had too little free disk. So unverified until someone runs the
  commands above: that the package builds with conda's compilers on linux-64, linux-aarch64, osx-64 and
  osx-arm64 (the C++ is plain C++17, built in CI with GCC, Clang, Apple Clang and MSVC), that
  `bioconda-utils lint` has no complaint (the `extra: identifiers: doi:TODO` and `recipe-maintainers` placeholders
  must be filled or removed first), and the run-time behaviour inside a conda environment.

What the package is: **one output, the command line and the Python package** (`coal`, `pal2nal`, `coal-cpp`,
`coal-cpp-align`, `import coal`). The static library `libcoal_core.a`, its headers and CMake package, which
`cmake --install` adds for C++ projects that embed the engine, are not shipped: nothing in Bioconda links against
them. If a C++ consumer appears, add a second output (`libcoal`) rather than putting them back into this one.
The recipe is not `noarch` (it holds compiled programs), so `python` is unpinned in `host:` / `run:` and the
Pythons below 3.10 are skipped with a selector; `python >=3.10` there is the `noarch: python` form and is wrong
for a compiled package.

An earlier version of the recipe (the PAL2NAL scope only):
on 2026-10-22 it (with `source:` pointed at a tarball of the `packaging/bioconda` branch, i.e. the split
layout) built and passed all of its tests with `conda build -c conda-forge` on linux-64 (python 3.15 variant); the
package held `bin/coal`, `bin/coal-cpp`, `bin/pal2nal`, the `coal` and `coal_compat` modules and
`share/coal/examples/`.

To test a not-yet-tagged checkout, point `source:` at a local tarball:
`git archive --format=tar.gz --prefix=COAL-<version>/ -o /tmp/coal.tar.gz HEAD`, then set
`url: file:///tmp/coal.tar.gz` and its sha256 in a copy of `meta.yaml` (do not commit that copy).

## The sha256 step

Every version bump needs the checksum of the GitHub tag tarball (the one in `source: url:`):

```
v=0.6.0   # the new tag, without the v in meta.yaml
curl -sL https://github.com/plumat-coal/COAL/archive/refs/tags/v$v.tar.gz | sha256sum
```

Put the version in `{% set version = "..." %}` and the hash in `{% set sha256 = "..." %}`, and reset
`build: number:` to 0. GitHub's auto-generated tag tarballs are stable for a given tag; do not retag a release.

## Submitting to Bioconda

1. Fork and clone https://github.com/bioconda/bioconda-recipes, create a branch.
2. Copy this folder to `recipes/coal-align/` (`meta.yaml`, `build.sh`, the two fixtures and `test_dispatch.sh`; **not**
   this README).
3. Fill in the TODOs in `meta.yaml`: `extra: recipe-maintainers:` (GitHub handles of the people who will review
   version bumps and answer BiocondaBot), `extra: identifiers:` (a DOI, or drop the line if there is none yet), and
   the release version + sha256.
4. Lint and build locally as above (at least `conda build`).
5. Open a pull request titled `Add coal-align 0.x.y`. Bioconda's CI builds Linux and macOS and runs the tests.
   When CI is green, comment `@BiocondaBot please add label` to request review; a Bioconda core member merges.
   The package then appears on https://anaconda.org/bioconda/coal-align within a few hours, and a Biocontainer
   (`quay.io/biocontainers/coal-align`) is built automatically.
6. For `linux-aarch64` and `osx-arm64` builds, add under `extra:`
   `additional-platforms: [linux-aarch64, osx-arm64]` (the code is plain C++17 and Python, so this should just work).

Name checks before submitting: `coal` is taken on PyPI by another project (hence `coal-align`); check that nothing
called `coal-align` or `coal` exists on anaconda.org (conda-forge, bioconda) at submission time. The `pal2nal` command
installed here does not clash with Bioconda's `pal2nal` package, which installs `pal2nal.pl`; both can be in one
environment.

## New releases

Bioconda's autobump bot opens a PR for a new version by itself when it finds a new **GitHub release tag** matching the
`source: url:` pattern (`v{{ version }}`), recomputes the sha256 and resets the build number. Tags alone are not enough
for every bot path: publish a GitHub Release for each tag (the repo's `release.yml` does this on `v*` tags). If a bump
needs a recipe change (new dependency, new test), edit `recipes/coal-align/meta.yaml` by hand in the bot's PR or in a
new one. A rebuild of the same version (recipe fix) increments `build: number:` instead.
