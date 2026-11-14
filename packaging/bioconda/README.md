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
| `build_local.sh` | Not submitted. `conda build` of this recipe against the current checkout (below, "Testing the recipe locally"). |

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

`build_local.sh` (in this folder; **not** part of the submission) builds the recipe with `conda build` against the
committed `HEAD` of the checkout: it makes a `git archive` tarball, copies the recipe to a work folder with `version`
set to the tree's version and `url` / `sha256` pointing at that tarball, and builds and tests the copy. The committed
`meta.yaml` keeps its placeholders. It needs `conda-build` on `PATH` (conda-build 26.9 was used; boa is not needed,
conda-build uses the libmamba solver).

`{{ stdlib('c') }}` only resolves with Bioconda's pinnings, so pass them as Bioconda's CI does: conda-forge's
`conda_build_config.yaml` (from the `conda-forge-pinning` package, the version bioconda-utils pins; it installs the
file at the root of its environment) and bioconda-utils' own `bioconda_utils-conda_build_config.yaml`
(`c_stdlib_version` 2.17 on Linux, 10.13 / 11.0 on macOS):

```
micromamba create -p ./cb -c conda-forge conda-build 'conda-forge-pinning==2026.01.07.13.16.30'
curl -sLO https://raw.githubusercontent.com/bioconda/bioconda-utils/master/bioconda_utils/bioconda_utils-conda_build_config.yaml
PATH=$PWD/cb/bin:$PATH bash packaging/bioconda/build_local.sh /tmp/coal-conda \
  -e cb/conda_build_config.yaml -e bioconda_utils-conda_build_config.yaml
# one Python only (Bioconda builds 3.10-3.13): add  -m py.yaml  with
#   python: ["3.10.* *_cpython"]
#   is_python_min: [true]
```

The full check Bioconda's CI runs is bioconda-utils inside a `bioconda-recipes` checkout (copy the recipe to
`recipes/coal-align/` first):

```
bioconda-utils lint recipes config.yml --packages coal-align
bioconda-utils build recipes config.yml --packages coal-align     # add --docker --mulled-test for the container test
```

Status of the recipe (2026-11-14):

- **`conda build` succeeds and all 24 test commands pass on linux-64** (conda-build 26.9.1, libmamba solver,
  `-c conda-forge -c bioconda`, strict channel priority, source = `git archive` of the tree at 0.6.0.dev0):
  - with Bioconda's pinnings (conda-forge-pinning 2026.01.07.13.16.30 + bioconda-utils' config): Python 3.10
    variant, `gxx_linux-64` 14.4, `sysroot_linux-64` 2.17. Package `coal-align-0.6.0.dev0-py310hd766df8_0.conda`,
    0.99 MB (3.6 MB installed). Build step 75 s, tests 4 s (1 job, `nice -n 10`, 2-core VM; the first run of
    `conda build` takes about 4.5 min more to fetch repodata and the compiler packages).
  - without pinnings (conda-build's defaults, before `stdlib('c')` was added): Python 3.12 variant, GCC 16, 1.13 MB.
  - Package contents: `bin/coal`, `bin/pal2nal`, `bin/coal-cpp`, `bin/coal-cpp-align`, the `coal` and `coal_compat`
    modules, `share/coal/examples/`; `info/licenses/` holds `LICENSE` and `NOTICE`. Depends on `libgcc`, `libstdcxx`,
    `python`, `python_abi`; `run_exports` = `coal-align >=0.6.0,<0.7.0a0` (weak, `max_pin="x.x"`); `coal-cpp` has
    `RPATH $ORIGIN/../lib` and links only `libstdc++`, `libm`, `libgcc_s`, `libc`.
  - Harmless conda-build messages: "Number of parsed outputs does not match detected raw metadata blocks"
    (conda/conda-build#5571, printed for every recipe without `outputs:`) and CMake's unused `CMAKE_*` variables
    from conda's toolchain.
- **`bioconda-utils lint`: "All checks OK"** (bioconda-utils `master` at 5a0a355, 2026-11-14, with the pinnings above). Before
  the fix it reported `compiler_needs_stdlib_c`. The linter does not check `recipe-maintainers`: fill it anyway.
- **Not built here:** Python 3.11 / 3.13 (only the site-packages path differs from 3.10; the binaries are the same),
  linux-aarch64, osx-64, osx-arm64 (Bioconda's CI builds them). macOS: `build.sh` defines
  `_LIBCPP_DISABLE_AVAILABILITY`, because Bioconda's osx-64 deployment target (10.13) predates `std::filesystem` in
  libc++ and `coal-cpp` uses it (conda-forge's documented fix; conda ships a modern libc++).


What the package is: **one output, the command line and the Python package** (`coal`, `pal2nal`, `coal-cpp`,
`coal-cpp-align`, `import coal`). The static library `libcoal_core.a`, its headers and CMake package, which
`cmake --install` adds for C++ projects that embed the engine, are not shipped: nothing in Bioconda links against
them. If a C++ consumer appears, add a second output (`libcoal`) rather than putting them back into this one.
The recipe is not `noarch` (it holds compiled programs), so `python` is unpinned in `host:` / `run:` and the
Pythons below 3.10 are skipped with a selector; `python >=3.10` there is the `noarch: python` form and is wrong
for a compiled package.

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
   this README or `build_local.sh`).
3. Fill in `meta.yaml`: the release version + sha256 (the two `PLACEHOLDER` lines), `extra: recipe-maintainers:`
   (replace `TODO-github-handle` with the GitHub handles of the people who will review version bumps and answer
   BiocondaBot), and, if COAL has a DOI by then, uncomment `extra: identifiers:`.
4. Lint and build locally as above (at least `conda build`).
5. Open a pull request titled `Add coal-align 0.x.y`. Bioconda's CI builds Linux and macOS and runs the tests.
   When CI is green, comment `@BiocondaBot please add label` to request review; a Bioconda core member merges.
   The package then appears on https://anaconda.org/bioconda/coal-align within a few hours, and a Biocontainer
   (`quay.io/biocontainers/coal-align`) is built automatically.
6. `extra: additional-platforms: [linux-aarch64, osx-arm64]` is already set, so the CI also builds ARM Linux and
   Apple silicon; if one of them fails and the fix is not quick, remove that platform rather than holding the PR.

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
