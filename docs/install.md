# Installing COAL

<!-- dev-version banner (tests/test_new_user_ux.py): remove it when 0.6.0 is released -->
> **These pages describe the development version, 0.6.0.dev0.** The latest release, 0.5.0, makes codon
> alignments only: `coal dnds`, `coal kaks`, `coal translate` and `coal pipeline` are not in it (its wheels answer
> `unrecognized arguments` or `cannot open 'pipeline'`). Until 0.6.0 is released, install from a checkout (Python
> 3.10+, a C++17 compiler, CMake 3.16+; about a minute):
>
> ```
> git clone https://github.com/plumat-coal/COAL && cd COAL
> python -m pip install -e .
> cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --config Release
> export COAL_CPP=$PWD/build/coal-cpp     # for this shell
> ```
>
> To keep `COAL_CPP` in new shells, and for what each install gives: [Option 3](#option-3-from-source) below.

**What each install gives.** COAL does three things (codon alignments; Ka/Ks for pairs: `coal kaks`, `coal dnds`;
the one-step `coal pipeline` from coding sequences to dN/dS), and not every install has all three:

| Install | Codon alignment | Ka/Ks (`coal kaks`, `coal dnds`) | `coal pipeline` |
|---|---|---|---|
| **Platform wheel** (`coal_align-<version>-py3-none-<platform>.whl`; [below](#the-c-engine-comes-with-the-wheel)) | yes, both engines | yes | yes, with DIAMOND (and Clustal Omega or MAFFT for clusters) |
| Source package (`coal_src.tar.gz`, the pip commands on this page), the `py3-none-any` wheel | yes, the Python engine | no | no |
| A checkout with `coal-cpp` built and `COAL_CPP` set ([Option 3](#option-3-from-source)) | yes, both engines | yes | yes, with DIAMOND |
| The Windows app, the web GUI | yes | no (the command line of the same install has what its wheel has) | no |
| R package `coalr` | yes | yes, in process | pairs and clusters you give, in process; the search and clustering through Python (reticulate) |

An install without `coal-cpp` has the Python engine only: it makes codon alignments, and `coal kaks`, `coal dnds`
and `coal pipeline` stop with a message that says so and names the platform wheel. `coal --cpp-path` tells which
kind an install is (it prints the path of `coal-cpp`, or exits 1). `coal-align` is **not on PyPI yet**
(publication is pending): every command on this page installs from a GitHub release, which until 0.6.0 is
released means 0.5.0, codon alignment only (the note at the top).

**DIAMOND and the aligners** are separate programs that only `coal pipeline` uses: DIAMOND 2.1.0 or later for the
ortholog search and the clustering (tested with 2.2.8), and Clustal Omega (or MAFFT) for clusters of three or more.
COAL looks for them on `PATH`, in `$COAL_DIAMOND` / `$COAL_CLUSTALO` / `$COAL_MAFFT` and in its own cache, and
installs them into the cache, checksum-verified, with

```
coal pipeline install diamond clustalo
```

(or `--install-tools` on a run; `conda install -c bioconda diamond clustalo` works as well). On Windows only
DIAMOND installs that way: for Clustal Omega and MAFFT see
[below](#clustal-omega-and-mafft-on-windows). Pair mode needs DIAMOND only.

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

Needs Python 3.10+, git, a C++17 compiler (gcc, clang, or Visual Studio 2019+) and CMake 3.16+ for `coal-cpp`,
which Ka/Ks, `coal translate` and `coal pipeline` need (the build takes about a minute):

```
git clone https://github.com/plumat-coal/COAL && cd COAL
python -m pip install -e .
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --config Release
export COAL_CPP=$PWD/build/coal-cpp
coal --cpp-path                      # prints the path of coal-cpp (exit 1: not found)
coal examples/clean/protein.aln.fasta examples/clean/dna.fasta > codon_aln.fasta
coal dnds examples/dnds/codon.fasta --pairs consecutive -o dnds.tsv
```

`export COAL_CPP` lasts for the current shell only. To make it permanent, copy the program next to the Python
that runs COAL, where `coal` looks for it ([cpp-engine.md](cpp-engine.md) "Where coal looks for coal-cpp"):
`cp build/coal-cpp "$(dirname "$(command -v python)")/"` (in a virtual environment: its `bin/`), or add the
`export COAL_CPP=...` line, with the full path, to `~/.bashrc` (or `~/.zshrc`). On Windows the program is
`build\Release\coal-cpp.exe`. Without the `cmake` lines, the checkout has the Python engine only: codon
alignments, and `coal dnds`, `coal kaks`, `coal translate` and `coal pipeline` say that they need `coal-cpp`.

`python -m pip install -e .` installs the engine only. Its optional extra `".[fast]"` adds a numpy speed-up for large inputs. The app is
the package in `apps/coal-app`, with two extras: `"apps/coal-app[gui]"` (web GUI and Windows app, Gradio) and
`"apps/coal-app[mcp]"` (stand-alone MCP server for AI agents). Install it next to the engine:

```
python -m pip install -e . -e "apps/coal-app[gui]"
```

**Install `gui` and `mcp` in separate virtual environments**: they pin conflicting versions of the `mcp` package, so
installing one into the other's environment can break it. The Windows launchers
(`apps\coal-app\desktop\run_from_source.bat`, `apps\coal-app\desktop\setup_mcp.bat`) do this for you. Quote
the extras (`"apps/coal-app[gui]"`): some shells (zsh, PowerShell) treat bare brackets specially.

## Clustal Omega and MAFFT on Windows

Cluster mode (`coal pipeline clusters`) aligns clusters of three or more genes with Clustal Omega or MAFFT. Neither
Bioconda nor conda-forge has a Windows package of either program (their package indexes, checked 2026-11-14).

**MAFFT: `coal pipeline install mafft` works on x86-64 Windows** (and `--install-tools` during a run). It installs
MAFFT's own Windows package, the "all-in-one" `mafft-7.526-win64-signed.zip` (MAFFT 7.526, the version COAL pins
elsewhere) from https://mafft.cbrc.jp/alignment/software/windows_without_cygwin.html, into COAL's cache
(`%LOCALAPPDATA%\coal\Cache\tools\mafft\mafft-win\`). MAFFT publishes no checksum, so COAL pins the file
itself: 27,476,513 bytes, SHA-256 `7c3a508c8512af04007d68019c6fbd4ff00d891244b58f87ec9aa0e53fbae8a7` (computed
on 2026-11-15 from two identical downloads). A download of another size or hash is deleted before anything is
unpacked, every name in the zip is checked before a file is written, and the install is recorded with the hash
of each file, as for the other tools. COAL then finds `mafft-win\mafft.bat` there; run with `--aligner mafft`.
COAL's CI downloads and runs it on its Windows runner, including one `coal pipeline clusters --aligner mafft`.

- **Licence:** MAFFT is BSD-licensed; the Windows package also bundles MSYS2 components (bash and its runtime),
  and MAFFT's page gives the GPL for the package as a whole. COAL does not redistribute it: each install
  downloads it from MAFFT's site.
- **How it runs:** COAL does not run `mafft.bat`. That batch file expands paths again, unquoted, inside
  `cmd.exe` (its own folder, the current folder, `TMP`, `PATH`), so a path with `&` or `^` in it would be read as
  a command there. COAL starts what the batch file starts, the package's `usr\bin\bash.exe` with MAFFT's script
  `/usr/bin/mafft`, with no shell and with the settings `mafft.bat` makes (`PATH`, `MAFFT_BINARIES`, a temporary
  folder of the run's own). Paths are then only text: COAL's CI runs it with the cache folder, the current folder
  and `TMP` / `TEMP` all in paths with spaces, `&`, `^`, `(` and `)`, and with an input given by a UNC path.
- **A MAFFT batch file of your own** (`COAL_MAFFT` naming a `.bat` that is not `mafft-win\mafft.bat` of MAFFT's
  package) runs through `cmd.exe`, so COAL refuses to start it, with a message, while its folder, the current
  folder or `TMP` / `TEMP` has one of `& | < > ^ % ! "` in its path. Spaces and `( )` are fine there (tested in
  CI with a stand-in batch file).
- Windows on ARM: the package is for x86-64; there the command says what works (below).

**Clustal Omega: no automatic install on Windows.** The pre-compiled Clustal Omega binaries of its authors' site
are gone ([alignment.md](alignment.md)), and COAL downloads a program only with a SHA-256 pinned in its source and
computed from the file itself. `coal pipeline install clustalo` says what works instead:

1. **WSL** (Windows Subsystem for Linux): run COAL there, and `coal pipeline install clustalo mafft` installs the
   pinned Bioconda / conda-forge builds as on Linux (the route COAL's CI tests, on Linux).
2. **Clustal Omega from MSYS2** (https://www.msys2.org), which packages Clustal Omega 1.2.4, the version COAL pins
   on Linux and macOS (`mingw-w64-ucrt-x86_64-clustal-omega` 1.2.4-4; pacman checks the package's signature). In an
   MSYS2 shell:

   ```
   pacman -S mingw-w64-ucrt-x86_64-clustal-omega
   ```

   then set `COAL_CLUSTALO` to `C:\msys64\ucrt64\bin\clustalo.exe` (MSYS2's default folder), or put
   `C:\msys64\ucrt64\bin` on `PATH`. Any other `clustalo.exe` you trust works the same way.
3. **MAFFT instead** (`--aligner mafft`), installed as above. A MAFFT Windows package you extracted yourself
   (keep its LF line ends; MAFFT's page warns that some unzip tools convert them) works too: set `COAL_MAFFT` to
   the full path of its `mafft-win\mafft.bat`.

Route 2, and a MAFFT of your own, are installs of your own: COAL neither downloads nor checks them, uses them as
they are (no version note, as for any program on `PATH`), and they are not tested with COAL.

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
