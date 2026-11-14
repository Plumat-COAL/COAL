# Installing COAL

**What each install gives.** COAL does three things (codon alignments; Ka/Ks for pairs: `coal kaks`, `coal dnds`;
the one-step `coal pipeline` from coding sequences to dN/dS), and not every install has all three:

| Install | Codon alignment | Ka/Ks (`coal kaks`, `coal dnds`) | `coal pipeline` |
|---|---|---|---|
| **Platform wheel** (`coal_align-<version>-py3-none-<platform>.whl`; [below](#the-c-engine-comes-with-the-wheel)) | yes, both engines | yes | yes, with DIAMOND (and Clustal Omega or MAFFT for clusters) |
| Source package (`coal_src.tar.gz`, the pip commands on this page), the `py3-none-any` wheel | yes, the Python engine | no | no |
| A checkout with `coal-cpp` built (`cmake`; README) and `COAL_CPP` set | yes, both engines | yes | yes, with DIAMOND |
| The Windows app, the web GUI | yes | no (the command line of the same install has what its wheel has) | no |
| R package `coalr` | yes | yes, in process | pairs and clusters you give, in process; the search and clustering through Python (reticulate) |

An install without `coal-cpp` has the Python engine only: it makes codon alignments, and `coal kaks`, `coal dnds`
and `coal pipeline` stop with a message that says so and names the platform wheel. `coal --cpp-path` tells which
kind an install is (it prints the path of `coal-cpp`, or exits 1). `coal-align` is **not on PyPI yet**
(publication is pending): every command on this page installs from a GitHub release.

**DIAMOND and the aligners** are separate programs that only `coal pipeline` uses: DIAMOND 2.1.0 or later for the
ortholog search and the clustering (tested with 2.2.8), and Clustal Omega (or MAFFT) for clusters of three or more.
COAL looks for them on `PATH`, in `$COAL_DIAMOND` / `$COAL_CLUSTALO` / `$COAL_MAFFT` and in its own cache, and
installs them into the cache, checksum-verified, with

```
coal pipeline install diamond clustalo
```

(or `--install-tools` on a run; `conda install -c bioconda diamond clustalo` works as well). On Windows there is no
automatic Clustal Omega or MAFFT install ([alignment.md](alignment.md)); pair mode needs DIAMOND only.

## Option 1: download the app from the website (recommended)

**https://plumat-coal.github.io/COAL/** has the Windows installer and step-by-step install instructions for
Windows, macOS and Linux. No terminal is needed on Windows.

- **Windows 10/11:** click **Download Windows installer** on the website (it installs for your user account; no
  administrator rights needed)
  ([direct link to the latest installer](https://github.com/plumat-coal/COAL/releases/latest/download/COAL-Setup-windows-x64.exe);
  all files, including a portable zip and `SHA256SUMS.txt`, are on the
  [Releases page](https://github.com/plumat-coal/COAL/releases/latest)). The app is not code-signed yet, so Windows
  shows **"Windows protected your PC"**: click **More info**, then **Run anyway**. The
  [install page](https://plumat-coal.github.io/COAL/install.html) walks through every screen, including what to do
  if Smart App Control blocks it. On the installer's last page, **Connect COAL to Claude Desktop** lets Claude run
  COAL for you ([AI agents](mcp.md)). User guide: `apps/coal-app/desktop/README.md`.
- **macOS and Linux:** follow the [install page](https://plumat-coal.github.io/COAL/install.html#macos) (three
  commands, also below).

## Option 2: macOS and Linux (or Windows with Python)

Needs Python 3.10+. COAL is two packages: the engine `coal-align` (the `coal` and `pal2nal` commands; no
dependencies) and the app `coal-align-app` (the web GUI and the MCP server; it installs the engine of the same version
too). Neither is on PyPI yet, so both are installed from the latest release:

```
python3 -m venv ~/coal-env
~/coal-env/bin/python -m pip install "coal-align @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_src.tar.gz" "coal-align-app[gui] @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_app_src.tar.gz"
~/coal-env/bin/coal-gui          # the app, in your browser at http://127.0.0.1:7860
~/coal-env/bin/coal --help       # the command line
```

For the command line alone (pipelines such as OrthologR or a DIAMOND workflow), install only the engine:
`pip install "coal-align @ https://github.com/plumat-coal/COAL/releases/latest/download/coal_src.tar.gz"`. The package is `coal-align`: `pip install coal` installs a
different project. Check downloads against
[SHA256SUMS.txt](https://github.com/plumat-coal/COAL/releases/latest/download/SHA256SUMS.txt).

### The C++ engine comes with the wheel

Every release also has one wheel per platform with the C++ engine `coal-cpp` inside, next to the pure-Python
wheel and the sdist:

| Wheel | For |
| --- | --- |
| `coal_align-<version>-py3-none-manylinux_2_28_x86_64.whl` | Linux x86_64, glibc 2.28+ (RHEL/Rocky 8+, Debian 10+, Ubuntu 18.10+ and anything newer) |
| `coal_align-<version>-py3-none-manylinux_2_28_aarch64.whl` | Linux aarch64 (ARM servers, Raspberry Pi OS 64-bit), glibc 2.28+ |
| `coal_align-<version>-py3-none-macosx_11_0_universal2.whl` | macOS 11+, Apple silicon and Intel (one file) |
| `coal_align-<version>-py3-none-win_amd64.whl` | Windows 10/11, 64-bit |
| `coal_align-<version>-py3-none-any.whl`, `coal_align-<version>.tar.gz` | everything else: the same `coal`, Python engine only |

From a release, download the wheel for your platform and
`python -m pip install coal_align-<version>-py3-none-<platform>.whl`; the `.tar.gz` command above installs the
Python engine only, which makes codon alignments but has no Ka/Ks and no pipeline. (Once `coal-align` is on PyPI,
`python -m pip install coal-align` will pick the right wheel by itself; publication is pending.) With the wheel
nothing else is needed for codon alignment and Ka/Ks: `coal` hands
inputs of 1 MB or more to the bundled `coal-cpp` ([cpp-engine.md](cpp-engine.md) "Where coal looks for
coal-cpp"), and `coal --cpp-path` prints its path (exit 1 if there is none), e.g. to call it directly from a
pipeline. A `coal-cpp` you built yourself still wins when `COAL_CPP` points at it.

## Conda (Bioconda): prepared, not yet published

A Bioconda recipe is ready under `packaging/bioconda/` (with `packaging/bioconda/README.md` on testing and
submitting it). Once it is merged into `bioconda-recipes`, `conda install -c conda-forge -c bioconda coal-align`
will install the engine (`coal`, `pal2nal`) and the C++ engine `coal-cpp` together on Linux and macOS, and `coal`
will find `coal-cpp` in the same environment with no settings. Until then, use the options above.

## Option 3: from source

```
git clone https://github.com/plumat-coal/COAL && cd COAL
python -m pip install -e .
coal examples/clean/protein.aln.fasta examples/clean/dna.fasta > codon_aln.fasta
```

That installs the engine only. Its optional extra `".[fast]"` adds a numpy speed-up for large inputs. The app is
the package in `apps/coal-app`, with two extras: `"apps/coal-app[gui]"` (web GUI and Windows app, Gradio) and
`"apps/coal-app[mcp]"` (stand-alone MCP server for AI agents). Install it next to the engine:

```
python -m pip install -e . -e "apps/coal-app[gui]"
```

**Install `gui` and `mcp` in separate virtual environments**: they pin conflicting versions of the `mcp` package, so
installing one into the other's environment can break it. The Windows launchers
(`apps\coal-app\desktop\run_from_source.bat`, `apps\coal-app\desktop\setup_mcp.bat`) do this for you. Quote
the extras (`"apps/coal-app[gui]"`): some shells (zsh, PowerShell) treat bare brackets specially.

## Ways to use COAL

- **Windows app:** installer from the [website](https://plumat-coal.github.io/COAL/); a local GUI, nothing is
  uploaded. From a source checkout, double-click `apps\coal-app\desktop\run_from_source.bat` (it creates
  `.venv-gui` on first start, repairs it if an install was interrupted, then opens the window). Building it:
  `apps/coal-app/desktop/BUILD.md`.
- **Web GUI on any OS:** `coal-gui` (or `python -m coal_app`, from `coal-align-app[gui]`) runs on this computer at
  http://127.0.0.1:7860. Its **Align** tab has a per-gene health check, an alignment view that marks where a gene
  breaks, genetic-code suggestions and a ready-to-paste methods sentence; its **Ka/Ks** tab gives dN, dS and dN/dS
  for pairs of a codon alignment or AXT file, and its **dN/dS pipeline** tab goes from the coding sequences of two
  species to dN/dS per ortholog pair (DIAMOND best or reciprocal best hits). Ka/Ks and the pipeline need `coal-cpp`
  (below), and the pipeline DIAMOND (`coal pipeline install diamond`).
- **Command line:** `coal` / `python -m coal`.
- **Drop-in for PAL2NAL pipelines:** `pal2nal pep.aln nuc.fasta [PAL2NAL options]` (or `python -m coal_compat`)
  accepts an unchanged PAL2NAL v14 command line, with PAL2NAL's default clustal output, all four formats,
  `-blockonly` and `-nostderr` (= `-q`); `-html` is a clear error.
- **AI agents (MCP):** four tools: `codon_align` and `explain_report` (codon alignment and its report),
  `dnds_pairs` (Ka/Ks for pairs of a codon alignment or AXT file) and `pipeline_pairs` (two species' coding
  sequences to dN/dS per ortholog pair); the last two need `coal-cpp`; see [mcp.md](mcp.md).
- **C++ engine for the largest jobs:** `coal-cpp`, same options and output ([cpp-engine.md](cpp-engine.md)); it
  comes inside the platform wheels (above) and the Windows app, or build it from source.
- **Two-species preset:** `demo/diamond_pipeline.py` (DIAMOND reciprocal best hits, then COAL on every pair).
- **Hosted demo:** a Hugging Face Space is prepared but not deployed yet: `apps/coal-app/space/README.md`,
  `apps/coal-app/space/DEPLOY.md`.
