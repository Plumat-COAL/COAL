# Installing COAL

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
  http://127.0.0.1:7860, with a
  per-gene health check, an alignment view that marks where a gene breaks, genetic-code suggestions and a
  ready-to-paste methods sentence.
- **Command line:** `coal` / `python -m coal`.
- **Drop-in for PAL2NAL pipelines:** `pal2nal pep.aln nuc.fasta [PAL2NAL options]` (or `python -m coal_compat`)
  accepts an unchanged PAL2NAL v14 command line, with PAL2NAL's default clustal output, all four formats,
  `-blockonly` and `-nostderr` (= `-q`); `-html` is a clear error.
- **AI agents (MCP):** two tools, `codon_align` and `explain_report`; see [mcp.md](mcp.md).
- **C++ engine for the largest jobs:** `coal-cpp`, same options and output ([cpp-engine.md](cpp-engine.md)).
- **Two-species preset:** `demo/diamond_pipeline.py` (DIAMOND reciprocal best hits, then COAL on every pair).
- **Hosted demo:** a Hugging Face Space is prepared but not deployed yet: `apps/coal-app/space/README.md`,
  `apps/coal-app/space/DEPLOY.md`.
