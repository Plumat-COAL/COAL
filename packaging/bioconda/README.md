# Bioconda recipe for `coal-align`

This folder is a ready-to-submit [Bioconda](https://bioconda.github.io/) recipe. Once it is merged into
`bioconda-recipes`, pipeline users get both engines with one command:

```
conda install -c conda-forge -c bioconda coal-align
coal --version        # Python engine, the `coal` and `pal2nal` commands, `import coal`
coal-cpp --version    # C++ engine; `coal` hands inputs >= 1 MB to it automatically
```

`coal` finds `coal-cpp` without any setting: `coal.dispatch.find_cpp()` looks at `COAL_CPP`, then next to the Python
interpreter running COAL (`$CONDA_PREFIX/bin`, where `build.sh` installs the binary), then on `PATH`. The recipe test
`test_dispatch.sh` checks both the discovery and that a 1 MB input really runs through `coal-cpp`.

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

The recipe installs the root `pyproject.toml` as the `coal-align` package (engine only, no dependencies) and runs
`cmake --install` when `CMakeLists.txt` has install rules (branch `feat/stdin-and-cpp-lib`; without them `build.sh`
copies `build-conda/coal-cpp` by hand). **The `v0.4.2` tag predates both**: its root package is still `coal` with the
app inside. The sha256 in `meta.yaml` is v0.4.2's, so the recipe is syntactically complete, but submit it for the first
release that contains the split layout (bump `version` and recompute the sha256, below).

## Testing the recipe locally

Without conda at all, the build steps can be run by hand against a virtual environment, which is what `find_cpp`
sees in a conda environment too (binary next to the interpreter):

```
python3 -m venv .venv-conda && .venv-conda/bin/pip install 'setuptools>=69.3'
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

Status: on 2026-10-22 this recipe (with `source:` pointed at a tarball of the `packaging/bioconda` branch, i.e. the split
layout) built and passed all of its tests with `conda build -c conda-forge` on linux-64 (python 3.15 variant); the
package held `bin/coal`, `bin/coal-cpp`, `bin/pal2nal`, the `coal` and `coal_compat` modules and
`share/coal/examples/`.

To test a not-yet-tagged checkout, point `source:` at a local tarball:
`git archive --format=tar.gz --prefix=COAL-<version>/ -o /tmp/coal.tar.gz HEAD`, then set
`url: file:///tmp/coal.tar.gz` and its sha256 in a copy of `meta.yaml` (do not commit that copy).

## The sha256 step

Every version bump needs the checksum of the GitHub tag tarball (the one in `source: url:`):

```
v=0.4.3   # the new tag, without the v in meta.yaml
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
