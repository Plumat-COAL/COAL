# COAL validation data

Test sets for measuring COAL's **accuracy** and **speed**, plus the scripts that build and score them.
The data itself lives in `validation/data/` (git-ignored); every set is rebuilt deterministically by
its script.

## Tiers

| Tier | Built by | What it is | Reference answer |
|---|---|---|---|
| `accuracy` | `build_sim.py accuracy` | 412 simulated gene families across 40 scenarios (edge cases, all 17 genetic codes) | exact (simulation truth) |
| `speed` | `build_sim.py speed` | single large alignments: 1k–1M sequences × 400 codons; 100–35,000 codons per gene | exact |
| `families` | `build_sim.py families` | 2,000 realistic families (2–500 seqs, 30–5,000 codons): per-family throughput | exact |
| `compara` | `build_compara.py` | 2,634 real vertebrate gene trees from Ensembl Compara 116 (of 54,308): up to 400 per size bin (2–10, 11–50, 51–200, 201–1000, all 289 trees >1000 seqs) plus trees with selenocysteine, `*`, B/Z and ambiguous bases | Ensembl's own CDS alignment |
| `ncbi` | `build_ncbi.py` | 3,437 real tree-of-life families from 316 RefSeq genomes (250 bacteria, 30 archaea, fungi, protists, plants, animals; genetic codes 11, 1, 4, 6, 12), aligned with FAMSA | none: reference-free checks + PAL2NAL |

### Simulation (`simulate.py`, `cases.py`)

Codons evolve down a coalescent tree under a GY94 codon model (κ, ω, F3x4 frequencies, gamma rate
variation) with whole-codon insertions and deletions, so the true codon alignment is known exactly and
gap patterns are tree-structured. Genetic codes come from Biopython's NCBI tables, not from COAL.
For the 1M-sequence speed set a faster "clades" mode is used (full tree for 512 clade ancestors, then
short independent branches).

Scenarios on top of the simulation (see `cases.py`): UTRs and stop codons to trim, protein `*` columns,
internal stops, alternative start codons, lower/mixed case, RNA, IUPAC ambiguity, frameshift digits
(PAL2NAL convention), wrong codons that must produce warnings (`expected_warnings.tsv`), shuffled IDs,
pairing by position, CDS split over two files, FASTA headers with descriptions, wrapped lines, CRLF line
endings, every NCBI table PAL2NAL supports, divergent/saturated sequences, very gappy alignments,
2,000-sequence families, titin-length genes (34,350 codons), very short genes and pairs.

Each case folder holds `pep.fasta`, `nuc.fasta` (or `nuc1/nuc2`), `truth.fasta` (when known) and
`meta.json` (scenario, options such as `-codontable`, sizes, flags). Each tier has a `manifest.jsonl`.

## Prebuilt data

The built sets are also packed as three archives (`coal_testset_sim.tar.zst`,
`coal_testset_ncbi.tar.zst`, `coal_testset_compara.tar.zst`). Unpack them into `validation/data/`:

```
cd validation/data && for f in /path/to/coal_testset_*.tar.zst; do tar --zstd -xf "$f"; done
```

The 1,000,000-sequence speed set is not in the archive (about 3 GB); rebuild it with
`python validation/build_sim.py speed` (about a minute; same seed, same files).

## Scoring (`check.py`, `run.py`)

`run.py` runs any tool over a tier, times it (median of `--repeats`), records peak memory (needs
`psutil`) and scores every output:

- **reference-free** (any data, any scale): rows have the protein-alignment IDs in order and equal
  length; every residue sits on a codon that encodes it (start codons for the first Met; `X` takes
  any codon; a mismatch is allowed only where the tool warned); gaps are gaps; frameshift digits have
  that many bases; removing the gaps gives a contiguous piece of the input CDS;
- **truth**: rows identical to the true codon alignment, and per-column agreement;
- **warnings**: mismatch warnings must be exactly the expected (ID, column) pairs;
- **PAL2NAL**: `--compare-to DIR` checks byte-identity with stored reference outputs.

```
pip install numpy scipy biopython psutil pyfamsa
python scripts/fetch_pal2nal.py

python validation/build_sim.py accuracy
python validation/build_sim.py families --count 2000
python validation/build_sim.py speed --max-n 100000          # drop --max-n for the 1M set (~3 GB)
python validation/build_compara.py --download && python validation/build_compara.py
python validation/build_ncbi.py

python validation/run.py accuracy --tool coal
python validation/run.py speed --tool coal --repeats 3
python validation/run.py families --tool pal2nal --keep-outputs validation/data/ref/pal2nal_families
python validation/run.py families --tool coal --compare-to validation/data/ref/pal2nal_families
python validation/run.py compara --tool coal
python validation/run.py ncbi --tool coal --compare-to validation/data/ref/pal2nal_ncbi
python validation/run.py accuracy --cmd "./coal-cpp {pep} {nuc} {opts} -o {out}" --name coal_cpp
```

Results go to `validation/results/<name>_<tier>_<time>.csv` (one row per case) with a `.txt` noting
the command, platform and CPU count. Only compare timings from the same machine.

## Notes and known conventions

- `table10` cases: PAL2NAL treats TGA as a stop in code 10 (COAL follows NCBI: Cys), so PAL2NAL's
  warnings fail the warnings check there. This is a documented deliberate difference.
- `iupac`: COAL accepts any ambiguity code whose expansions all give the residue; PAL2NAL accepts a
  fixed list.
- Compara: Ensembl's CDS alignment is a back-translation of its protein alignment, so it is an
  independent check of placement. `meta.json` → `features` flags trees with selenocysteine `U`, `X`,
  B/Z, `*` and ambiguous bases. Checked against its own inputs, Ensembl's answer passes the
  reference-free checks in all but a handful of trees, where Ensembl has a Met on a non-start codon
  or a residue on a stop codon (readthrough). A tool should warn there, not silently agree.
- NCBI: CDS come straight from RefSeq annotation, with stop codons, partial genes, `transl_except`
  (selenocysteine) and programmed frameshifts (`exception=ribosomal slippage`). The members with such
  quirks are listed in `meta.json` → `quirky_members`. The genetic code of each genome is inferred
  from its data (`genomes.json`).
