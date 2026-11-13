# Ka/Ks (dN/dS): `coal kaks` and `coal dnds`

COAL estimates Ka (dN) and Ks (dS) for pairs of aligned coding sequences in process, with a C++ engine
(`cpp/src/kaks.cpp`) derived from KaKs_Calculator 2.0 and PAML's yn00 (see "Where the code comes from" below). No
external program runs. It has two commands:

- **`coal kaks`** does the job of `KaKs_Calculator -i IN.axt -o OUT.kaks -m METHOD` for seven of KaKs_Calculator's
  methods (NG, LWL, MLWL, LPB, MLPB, YN, MYN). It takes KaKs_Calculator's `-i`, `-o`, `-m` and `-c` and writes its
  tab-separated format: the same 22 columns, the same NA rules and 6 significant digits. It is not a drop-in for
  every KaKs_Calculator command line: the other methods are refused, and without `-m` it computes YN, where
  KaKs_Calculator computes MA (model averaging).
- **`coal dnds`** writes a tidy TSV from the same engine, with one row per pair and method, plus standard errors,
  a status and a note.

```
coal kaks -i pairs.axt -o pairs.kaks                     # YN00 (the default)
coal kaks -i pairs.axt -o pairs.kaks -m NG -m YN -m MYN  # several methods; rows come in KaKs_Calculator's order
coal dnds codon_aln.fasta -o dnds.tsv                    # every pair of an aligned FASTA (MSA), YN00
coal dnds pairs.fasta --pairs consecutive -m YN -m NG -j 8 -c 2
```

Both commands are subcommands of `coal-cpp`, and `coal` (Python) hands them to the `coal-cpp` of its own version
(the lookup of `coal --cpp-path`). The platform wheels and the Windows app include `coal-cpp`. A pure-Python
install has no Ka/Ks engine; `coal kaks` then says so and exits 1.

## Methods

| `-m` | Method | Reference | Standard errors (`coal dnds`) |
|---|---|---|---|
| `YN` (default) | Yang & Nielsen (2000): kappa from non-degenerate and 4-fold sites, codon frequencies F3x4, pathways weighted by P(t) and iterated, F84 distances with K80/JC69 fallbacks | Mol Biol Evol 17:32-43 | yes (F84 / K80 / JC69 variance, as PAML) |
| `NG` | Nei & Gojobori (1986) as KaKs_Calculator counts sites (see the note below), Jukes-Cantor correction | Mol Biol Evol 3:418-426 | yes (Jukes-Cantor variance) |
| `LWL` | Li, Wu & Luo (1985): 0-, 2- and 4-fold sites, K2P per class | Mol Biol Evol 2:150-174 | NA |
| `MLWL` | LWL with Tzeng, Pan & Li's (2004) kappa-based sites and Ile/Arg rules | Mol Biol Evol 21:2290-2298 | NA |
| `LPB` | Li (1993); Pamilo & Bianchi (1993) | J Mol Evol 36:96-99; Mol Biol Evol 10:271-281 | NA |
| `MLPB` | LPB with Tzeng et al.'s (2004) Ile/Arg rules | Mol Biol Evol 21:2290-2298 | NA |
| `MYN` | Zhang, Li & Yu (2006): YN00 with separate transition rates for purines and pyrimidines (TN93) | BMC Evol Biol 6:44 | yes (TN93 variance; F84 fallback) |

**`NG` is KaKs_Calculator's NG, not the textbook Nei-Gojobori count.** KaKs_Calculator takes a codon's synonymous
sites from its first and third positions only, leaves the changes that lead to a stop codon out of the
nonsynonymous sites, and then rescales S + N to the compared length; COAL does the same, so that `-m NG` gives
KaKs_Calculator's numbers. The classical NG86 procedure (as in MEGA, or PAML's NG86 output) counts the sites of
all three positions without that rescaling, and its dS differs: by up to about 0.02 on the 999 orthologr pairs
(measured by the independent review of 2026-11-12). Use `YN` (the default) when the absolute values matter, and
do not mix NG values of COAL / KaKs_Calculator with NG86 values of other programs.

The other KaKs_Calculator methods are refused with exit status 2 and the message "COAL does not implement
natively". These are `GY` (maximum likelihood), `MS` and `MA` (model selection and averaging), `NONE`, the gamma
variants `GNG GLWL GLPB GMLWL GMLPB GYN GMYN`, and `ALL`. KaKs_Calculator's default method is MA. COAL's default is
YN, so `coal kaks` without `-m` gives YN rows.

**Where the code comes from.** The engine is not an independent implementation: `cpp/src/kaks.cpp` is derived
from the source code of KaKs_Calculator 2.0 (Copyright (C) 2009 Beijing Institute of Genomics of CAS; Da-Peng Wang,
Yu-Bin Zhang, Zhang Zhang, Kristian K Ullrich) and of PAML's `yn00.c` (Copyright 1998 Ziheng Yang), both under the
GNU GPL version 3. Its functions follow theirs in structure, constants and order of operations, which is why the
numbers are the same; they were rewritten without classes and global state, for any of COAL's genetic codes and
for threads, and changed where listed under "Differences from KaKs_Calculator". The header of `kaks.cpp` and the
`NOTICE` file give the details. This is the reason COAL is distributed under GPL-3.0-only from this release on.
Neither program is included as such, and `coal kaks` / `coal dnds` do not run them.

## Policy

COAL applies the same rules to every pair and every method. When a pair cannot be used, it gets a status
(below), and the rest of the run carries on. One bad pair never stops a run.

| Situation | What COAL does | Status |
|---|---|---|
| Stop codon as the last codon of a sequence (its last codon that is not a gap) | that column is dropped | `ok` (not a warning) |
| Gap (`-` or `.`) anywhere in a codon of either sequence | that column is skipped (pairwise deletion) | `ok` |
| A letter other than A/C/G/T/U in a codon of either sequence (N, IUPAC, X) | that column is skipped | `ok` |
| Stop codon anywhere else (internal) in either sequence | that column is skipped; the note gives the number | `warning` |
| YN or MYN: the iteration reaches its round limit (10 or 100) without converging (short or nearly saturated pairs; "Differences", 12) | no estimate: Ka, Ks, Ka/Ks, the standard errors and the weighted difference counts (Sd, Nd) are NA. The site counts stay. The note gives the last round's Ka and Ks, which are not estimates ("Platform reproducibility") | `not_converged` |
| Multiple-hit correction undefined (log of a non-positive number; for YN00/MYN, F84, K80 and JC69 all fail, or P + Q > 1) | that estimate is NA; the note says which (Ka, Ks) | `saturated` |
| No codon column left | nothing is estimated | `no_data` |
| No synonymous site (e.g. Met/Trp codons only), or no nonsynonymous site (e.g. `GCT` against `GCA`; 200 × `AAA` against `AAG`) | that estimate (Ks, or Ka) and Ka/Ks are NA; **the other one is still estimated**. The note says which sites are missing | `no_data` |
| Rows of different length, or a length that is not a multiple of 3 | nothing is estimated | `error` |

Lower case and RNA (`U`) are read as upper-case DNA. Codon columns are counted on the alignment, so a codon that a
gap splits (`A-G`) counts as a gap. When several statuses apply, the most serious one wins
(error > no_data > saturated > not_converged > warning > ok), and the note lists all of them. One exception: an
iteration that did not converge is `not_converged` whether or not its last round was saturated, because that
differs between platforms too (its note says what the last round was). The genetic code is set with `-c` /
`--codontable` and takes the NCBI table numbers that COAL's `--codontable` accepts (1-6, 9-16, 21-23).

## Output

**`coal kaks`** writes KaKs_Calculator 2.0's header (`KaKs.cpp`) verbatim:

```
Sequence  Method  Ka  Ks  Ka/Ks  P-Value(Fisher)  Length  S-Sites  N-Sites  Fold-Sites(0:2:4)  Substitutions
S-Substitutions  N-Substitutions  Fold-S-Substitutions(0:2:4)  Fold-N-Substitutions(0:2:4)  Divergence-Time
Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)  GC(1:2:3)  ML-Score  AICc  Akaike-Weight  Model
```

The header is the same in KaKs_Calculator 3.0's documentation; KaKs_Calculator 3.0 itself was not run (its licence
allows academic use only, and it has no public repository). Values follow KaKs_Calculator's rules:

- Ka, Ks, the site counts and Divergence-Time below 1e-6 print as `NA`. So Ka = 0 prints as `NA` there, as in
  KaKs_Calculator.
- `Length` is the number of nucleotides compared.
- `Fold-*` are given for the LWL family only.
- Substitution-Rate-Ratio is `kappa:kappa:1:1:1:1`. For MYN it is `kappa_TC:kappa_AG:1:1:1:1`, and for NG
  `1:1:1:1:1:1`.
- `GC` is computed over both rows as given, with gaps included in the length.
- ML-Score, AICc, Akaike-Weight and Model are `NA`; they belong to the GY, MS and MA methods.

Rows come per pair in KaKs_Calculator's method order (NG, LWL, MLWL, LPB, MLPB, YN, MYN). Pairs with status
`error` are left out, as KaKs_Calculator leaves out pairs it cannot read, and the summary on stderr counts them.
This format has no status column. A YN or MYN row whose iteration did not converge has `NA` in Ka, Ks, Ka/Ks,
P-Value(Fisher), S-Substitutions, N-Substitutions and Divergence-Time (KaKs_Calculator prints the values of its
last round there: difference 12); the summary on stderr says how many such rows there are, and `coal dnds` names
them.
`--digits N` sets the number of significant digits.

**`coal dnds`** (`--digits` default 8):

| Column | Meaning |
|---|---|
| `seq1`, `seq2` | the two IDs. For FASTA input, the record IDs. For AXT, the name line split at the first `&`; when there is no `&`, seq1 is the whole name line and seq2 is `NA` |
| `method` | NG, LWL, MLWL, LPB, MLPB, YN or MYN |
| `dN`, `dS`, `dNdS` | Ka, Ks, Ka/Ks. `NA` when not available; dNdS is `NA` when dS < 1e-6. Zeros print as 0 |
| `SE_dN`, `SE_dS` | standard errors (YN, MYN, NG); `NA` elsewhere |
| `N_sites`, `S_sites` | nonsynonymous and synonymous sites |
| `Nd`, `Sd` | nonsynonymous and synonymous differences |
| `codons_compared` | codon columns used after the policy |
| `status`, `note` | see Policy: `ok`, `warning`, `not_converged`, `saturated`, `no_data`, `error` |

**Input.** The format is detected from the first non-blank character (`>` means FASTA, anything else AXT), or set
with `--format`.

- **AXT** is KaKs_Calculator's format: a name line, the two aligned sequences, then a blank line. When a record
  has more than two sequence lines, they are joined and split in half, as KaKs_Calculator does. A record with
  **one** sequence line (the second sequence is missing) is an error for that record: `coal dnds` gives it status
  `error`, `coal kaks` leaves it out with a warning, and the other records are computed. (KaKs_Calculator would
  cut the one line in half and compare the halves.) Lines starting with `#` before a record are skipped.
- **Aligned FASTA**: `--pairs all` (the default) compares every pair of the alignment in AXTConvertor's order
  (1-2, 1-3, ... 2-3, ...). `--pairs consecutive` compares records 1 and 2, 3 and 4, and so on. In `coal kaks`, a
  FASTA pair is named `id1&id2`, as KaKs_Calculator's AXTConvertor names it.
- Standard input (`-`) is not read; give a file.
- **The output must not be the input.** `-o` naming the input file (by any spelling: `./x`, a symbolic or hard
  link) is refused with exit status 2 before anything is written; the output is opened while the input is still
  being read, so the input would be destroyed. The same holds for `coal translate` and `coal codon-dnds`.

**Exit status.** 0 when the run finished, even if some pairs have a status other than `ok`. 1 when the input
cannot be read or the output cannot be written. 2 for a command-line error, including a method that is not
implemented. `-j N` sets the number of threads (default: all hardware threads). Pairs are computed in parallel
and written in input order, and the output does not depend on `-j`.

`--yn-flavour paml` and `--yn-weighting 0|1` choose which reference YN00 follows where the two differ (below).
`--yn-weighting 0` gives equal pathway weights in one round, which is the default control file of PAML's yn00
(`weighting = 0`). KaKs_Calculator always weights pathways.

## C++ API

`#include <coal/kaks.hpp>` (installed with `coal::core`, docs/cpp-engine.md). The library writes nothing to
stdout or stderr and has no global or static mutable state, and it calls no C library function that has one:
the log-factorials of Fisher's test are computed with the engine's own log-gamma, because `lgamma()` writes the
global `signgam` (a data race between threads, which ThreadSanitizer reported for `coal kaks -j 4` before this was
changed; it reports none now). It never throws for bad data: problems are reported in `Result::status` and
`Result::note`. It can be called from many threads at once.

```cpp
#include <coal/kaks.hpp>
const coal::GeneticCode* code = coal::get_code(1);            // NCBI table; nullptr if not supported
coal::kaks::Result r = coal::kaks::compute(row1, row2, coal::kaks::Method::YN, *code);
// r.ka, r.ks, r.ka_ks, r.se_ka, r.se_ks (NaN = NA), r.S, r.N, r.Sd, r.Nd, r.codons_compared,
// r.status (ok / warning / not_converged / saturated / no_data / error), r.note

// several methods for one pair: apply the policy once
coal::kaks::Prepared p = coal::kaks::prepare(row1, row2, *code);
for (auto m : {coal::kaks::Method::YN, coal::kaks::Method::NG}) use(coal::kaks::compute(p, m, *code));

coal::kaks::Method m; std::string why;
if (!coal::kaks::parse_method("GY", m, why)) { /* why: "... does not implement natively ..." */ }
```

`Options{yn_flavour, yn_weighting}` is the fourth argument (the defaults follow KaKs_Calculator).
`Result` also carries the fields KaKs_Calculator prints: `p_fisher`, `divergence`, `rates[6]`, `fold_*`, `gc[4]`,
`substitutions` and `length`. It also has the policy counts (`gap_codons`, `ambiguous_codons`, `internal_stops`,
`terminal_stop_dropped`) and, for YN/MYN, `iterations`, `converged` and, when the iteration did not converge
(`ka` and `ks` are then NaN), `last_ka` and `last_ks`: the last round's values, for inspection only.

## Validation

The validation compares COAL with two references:

- **KaKs_Calculator 2.0.1** (bioconda `kakscalculator2`, the build orthologr embeds as gestimator) in two builds.
  The stock binary prints 6 significant digits. A local build of the same source with `stream.precision(15)` in
  `CONVERT` was used for the numeric comparison.
- **PAML 4.10.10 yn00**, built from source with the `rst` line printed at `%.12g`; the stock binary prints 4
  decimals. It was run once per pair, because yn00 carries t over from one data set to the next, with
  `weighting = 1` and `weighting = 0`. PAML refuses stop codons, so columns with a gap, an ambiguous base or a stop
  were removed first, which matches COAL's policy.

The scripts are in `benchmarks/kaks/`:

- `make_orthologr_pairs.py` builds the orthologr set.
- `simulate_pairs.py` makes known-omega pairs.
- `validate_kaks.py` runs the comparisons and classifies every difference.
- `known_omega.py` summarises the simulations.
- `bench_kaks.py` measures speed.

Results of 2026-11-11 (branch `feat/kaks-engine`), run again on 2026-11-12 after the changes of the branch
`fix/review-findings` (thread-safe log-gamma, differences 8 to 10): every count and every maximum difference in the
tables below is unchanged; only the Fisher p-values moved, in their 12th digit. **On these real genes the review
fixes changed nothing else:** the `.kaks` output for the 999 orthologr and 600 simulated pairs is byte-identical
before and after them (checked again by the independent review of 2026-11-13). The change of 2026-11-13 (an
iteration that did not converge gives no estimate, difference 12) changed one row of these 11,193: the MYN row of
the simulated pair `sim_w0.5_t2_k2_31`, whose iteration does not converge in 100 rounds. KaKs_Calculator prints
its last round there (Ka 0.58117, Ks 1.43957, as COAL did), COAL now prints NA. The orthologr file is
byte-identical, and so is every other row of the simulated one. `validate_kaks.py` was run again on both sets
after that change: every count and maximum in the tables below is as before, apart from that row.

**Scope of each statement.** The tables below are about real genes of normal length (999 pairs) and simulated
genes of 500 codons (600 pairs). They say nothing about sequences of a few codons. Those are covered by two sets
of random short pairs, and there the picture is different: see "Short and saturated pairs" below.

**What was compared, and how.** "KaKs_Calculator" is always the program run on the same AXT file as COAL, method by
method. "PAML" is `yn00` run on one pair at a time, on a PHYLIP file from which the columns with a gap, an
ambiguous base or a stop codon had been removed, read from its `rst` file at 12 digits, and compared with
`coal dnds -m YN --yn-flavour paml --yn-weighting 1` (and `0`). A difference is counted wherever both programs
print a number; where one prints a number and the other does not, the pair is in the "NA in one" counts of
`validate_kaks.py`'s report and belongs to one of the documented differences. So "PAML to 5e-12" means: the
largest absolute difference in dN or dS over these pairs, with numbers that PAML prints to 12 digits.

**Data sets**

- **orthologr 1k.** orthologr's example genes `inst/seqs/ortho_{thal,lyra}_cds_1000.fasta` (A. thaliana and
  A. lyrata) were paired by DIAMOND 2.2.8 reciprocal best hits, giving 999 pairs. The proteins were aligned
  globally (Biopython PairwiseAligner, BLOSUM62, gap open 11, extend 1, full length, no end truncation) and turned
  into codon alignments by coal-cpp.
- **Simulated.** 600 pairs of 500 codons under GY94 with kappa = 2 and unequal F3x4 frequencies: omega 0.1, 0.5,
  1 and 2, crossed with t = 0.3, 1 and 2 substitutions per codon, 50 replicates each.
- **Test fixtures** in `tests/kaks/`: 23 hand-made edge cases and 48 simulated pairs, with the reference outputs
  stored in `tests/kaks/golden/`.
- **Random short pairs, one changed position per codon** (`benchmarks/kaks/random_short_pairs.py --mode single`,
  seed 5): 4,000 pairs of 2 to 8 codons over full and reduced alphabets. They reach sites of one kind only, base
  frequencies of zero and a kappa that cannot be estimated (differences 8 to 10), but no codon differs at more
  than one position, so no pathway is ever weighted.
- **Random short pairs, codons drawn independently** (`--mode multi`, seed 20261112; the generator of the
  independent review of 2026-11-13): 3,000 pairs of 1 to 12 codons from ten codon pools. Codons differ at two and
  three positions, one kind of site is often saturated, and the YN / MYN iteration often does not converge
  (differences 11 to 13).

**COAL compared with KaKs_Calculator 2.0, 15 significant digits.** "Max |diff|" is over pairs where both tools
give a number.

| Method | Set | Ka max abs diff | Ks max abs diff | Ka/Ks max abs diff | Pearson r (Ka, Ks) | Pairs > 1e-4 |
|---|---|---:|---:|---:|---:|---:|
| NG, LWL, MLWL, LPB, MLPB | orthologr 999 | 0 | 0 | 0 | 1 | 0 |
| YN | orthologr 999 | 1.0e-15 | 1.0e-15 | 1.0e-14 | 1 | 0 |
| MYN | orthologr 999 | 1.0e-15 | 1.0e-14 | 2.0e-14 | 1 | 0 |
| NG, LWL, MLWL, LPB, MLPB | simulated 600 | 0 | 0 | 0 | 1 | 0 |
| YN | simulated 600 | 1.0e-15 | 1.4e-13 | 1.0e-14 | 1 | 0 |
| MYN | simulated 600 | 1.0e-15 | 2.9e-14 | 1.0e-14 | 1 | 0 |

Fisher's p-values agree to 7e-12 (2.4e-12 on the orthologr pairs, 6.7e-12 on the simulated ones) apart from the
two documented cases below (`p-obs`, `p-neg`). On the 4,000 random short pairs with one changed position per
codon, no value that both tools print differs by more than 1e-9 (Ka, Ks, Ka/Ks and the site counts, all seven
methods). **That is not true for short pairs in general:** on the 3,000 pairs with codons drawn independently,
YN and MYN print different numbers in 23 fields (Ka, Ks or Ka/Ks; NG and the LWL / LPB family in none), which
belong to differences 11 and 13 below, and in 210 rows (123 YN, 87 MYN) COAL prints no estimate at all because
the iteration did not converge (difference 12), where KaKs_Calculator prints its last round.

**COAL `--yn-flavour paml` compared with PAML yn00, 12 significant digits.** The result is the same for
weighting 1 and 0.

| Set | dN max abs diff | dS max abs diff | S, N sites max abs diff | SE max abs diff | Pairs > 1e-4 |
|---|---:|---:|---:|---:|---:|
| orthologr 999 | 4.8e-13 | 5.0e-13 | 5.0e-9 | 4.7e-13 | 0 |
| simulated 600 | 5.0e-13 | 5.0e-12 | 5.0e-9 | 4.5e-11 | 0 |

These differences are the 12-digit printing of yn00.

**Text of every row compared with the stock KaKs_Calculator, 6 digits.** Every differing field is classified;
the categories are explained in the next section.

| Set | Method | Rows identical | -0 | nan | p-obs | p-neg | saturated | not converged | Unexplained |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| orthologr | NG | 965 / 999 | 11 | 0 | 22 | 0 | 1 | 0 | 0 |
| orthologr | LWL / MLWL | 998 / 999 | 0 | 0 | 0 | 1 | 1 | 0 | 0 |
| orthologr | LPB / MLPB | 992 / 987 of 999 | 0 | 6 / 11 | 0 | 0 | 1 | 0 | 0 |
| orthologr | YN | 970 / 999 | 11 | 0 | 17 | 0 | 1 | 0 | 0 |
| orthologr | MYN | 975 / 999 | 11 | 0 | 13 | 0 | 0 | 0 | 0 |
| simulated | NG | 596 / 600 | 0 | 0 | 4 | 0 | 0 | 0 | 0 |
| simulated | LWL / MLWL / LPB / MLPB | 480 / 495 / 528 / 515 of 600 | 0 | 0 | 0 | 118 / 103 / 70 / 83 | 2 | 0 | 0 |
| simulated | YN | 599 / 600 | 0 | 0 | 0 | 0 | 1 | 0 | 0 |
| simulated | MYN | 598 / 600 | 0 | 0 | 0 | 0 | 1 | 1 | 0 |

### Differences from KaKs_Calculator

Differences 1 to 10 are deliberate and were each checked against the references. Differences 11 to 13 are about
short or saturated pairs on which the two programs print different numbers; for each, the text says which
program is right, or that neither number is an estimate.

1. **`saturated`.** KaKs_Calculator replaces an undefined correction with a number. YN and MYN (when TN93 fails)
   report JC69's capped distance: when P >= 0.75, P is set to 0.75(n-1)/n, or t to 99 when P + Q > 1. The LWL
   family sets the class's K to 0, which underestimates Ks. COAL gives NA with status `saturated`. Example: the
   orthologr pair AT1G04625.1, whose proportion of synonymous differences is about 0.8. There KaKs_Calculator
   reports YN Ks = 3.69 and LWL Ks = 1.24, and COAL reports NA for both. PAML yn00 reports the same capped
   values, so the same pair is the one "NA in one" row against PAML.
2. **`p-obs`.** KaKs_Calculator's Fisher test sums the tables that are "no more probable than the observed one".
   Its log-factorial approximation sometimes makes the observed table itself appear 1e-13 more probable than
   itself, so it is left out. COAL counts ties within a relative 1e-7, as R's `fisher.test` does. Every such case
   was checked: COAL's p minus KaKs_Calculator's p equals P(observed table). This affects 52 of 7 × 999 orthologr
   rows.
3. **`p-neg`.** LPB93 sets S = Sd/Ks, so once Ks > 1 the table [Sd, S-Sd; Nd, N-Nd] has a negative cell. The LWL
   family hits the same once Sd > S. KaKs_Calculator still computes a p-value, and its log-factorial returns 1 for
   the negative cell. COAL gives NA for the test only; Ka and Ks are unaffected. This is common in the
   high-divergence simulations.
4. **`-0`.** KaKs_Calculator prints Ka/Ks = `-0` when Ka is a negative zero. COAL prints `0`.
5. **`nan`.** Where KaKs_Calculator prints `nan` or `-nan`, COAL prints `NA`. This includes the `S-Sites` and
   `N-Sites` columns of LPB and MLPB: these methods have no site count of their own, KaKs_Calculator derives
   S = Sd/Ks and N = Nd/Ka for the output, and prints `-nan` when the divisor is 0 or NA (identical sequences,
   Ka = 0, a saturated Ks); COAL prints `NA` there. It also covers pairs with nothing to compare.
6. **RNA and lower case.** KaKs_Calculator's NG and LWL family compare codon letters, so `U` in one sequence
   against `T` in the other counts as a difference; its YN and MYN read U as T. KaKs_Calculator's GC column counts
   upper-case G and C only. COAL reads U as T and ignores case everywhere (edge case `lowercase_rna`).
7. **Pairs that KaKs_Calculator cannot read** (a length that is not a multiple of 3): both tools leave them out of
   the `.kaks` file. `coal dnds` lists them with status `error`.

8. **YN and MYN with sites of one kind only.** With no nonsynonymous sites (N = 0: `GCT` against `GCA`, where F3x4
   gives weight to synonymous neighbours only; 200 × `AAA` against `AAG`) Ks is still defined, and with no
   synonymous sites (S = 0: Met and Trp codons) Ka is. COAL estimates the one that has sites; the other is NA and
   the status is `no_data`. (Until the review of 2026-11-12 COAL gave NA for both.) When no codon differs at more
   than one position, KaKs_Calculator's YN and PAML's yn00 give the same number: Ks = 0.477386 for `GCT`/`GCA` and
   0.44084 for the Lys pair in all three programs. When a codon differs at two or three positions the references
   go wrong, because they add the NaN of the missing kind into t and the pathway weights P(t) become NaN: yn00
   prints `nan` with `weighting = 1`; KaKs_Calculator prints `NA` or `nan`, **or a number that is too small**
   (difference 11). COAL takes t from the kind that exists; omega has no effect there, since every change between
   codons of non-zero frequency is of that one kind. Example `ATGTGGATGTGG` against `ATGTGGTGGATG`: COAL's YN
   gives Ka = 0.513346, which is what yn00 prints with `weighting = 0` (0.513345885278); KaKs_Calculator prints
   `NA`.
9. **MYN: `nan` rows from the mean kappa.** KaKs_Calculator estimates kappa(TC) and kappa(AG) from two classes of
   sites (non-degenerate and four-fold) and gives a class weight 0 when its kappa cannot be estimated, but still
   adds `0 × NaN` to the weighted mean (and reads an uninitialised value when the class has no site at all). The
   mean kappa is then NaN and every field of the row is `nan`. COAL leaves a class of weight 0 out of the mean, as
   the weight says. Fixture `myn_nan_kappa`: COAL Ka = 0.136556, Ks = 0.971172.
10. **MYN: `nan` from TN93 with a missing base pair, or from undefined site frequencies.** The TN93 distance
    divides by the products of the frequencies of T and C and of A and G at the sites. When the sites have no T-C
    pair or no A-G pair that product is 0: KaKs_Calculator tests only for both being 0, divides 0 by 0 when there
    is no difference of that kind, and prints `nan` (with a difference of that kind it divides by 0 into `-inf`,
    takes that as "TN93 fails" and falls back to F84). With one transition kind impossible TN93 is exactly F84, so
    COAL computes F84 there. An identical short pair (`GCTAAAGGT`) and `GCTGCTGCTAAAAAA` against
    `GCAGCCGCGAAGAAG` (no nonsynonymous difference) have Ka = 0 in COAL and `nan` in KaKs_Calculator; before the
    review COAL reported "Ka saturated" for them. Likewise, the base frequencies at the synonymous
    (nonsynonymous) sites are the mean over the two sequences, and a sequence without any such site contributes
    0/0 in KaKs_Calculator and yn00; COAL averages over the sequences that have such sites (fixture
    `myn_sites_in_one_sequence`: Ks = 0, as KaKs_Calculator's own YN says).

Differences 9 and 10 do not occur in the 999 orthologr and 600 simulated pairs (genes of normal length have every
base at their sites), which is why the tables above have no `nan` for YN and MYN; in the 4,000 random pairs with
one changed position per codon KaKs_Calculator's MYN prints `nan` for the sites of about a quarter of the pairs.
A KaKs_Calculator built with two lines changed (the mean over the classes of non-zero weight; "TN93 fails" when
either product is 0), called "the corrected reference" below, gives COAL's numbers for every one of those 4,000
pairs where it prints a number (0 differences above 1e-9). On the 3,000 pairs with codons drawn independently it
does not: see the table under difference 13.

#### Short and saturated pairs: where the two programs print different numbers

These three differences need a codon that differs at two or three positions, so that pathways are weighted by
P(t). None occurs in the 999 orthologr and 600 simulated pairs. The examples are in `tests/test_kaks.py`
(`TestMultiPositionCodons`) with the reference values.

11. **KaKs_Calculator drops the differences of codons that differ at several positions once a NaN is in its
    iteration, and prints the rest as if it were the estimate. COAL is right here.** The pathway weights of such
    a codon come from P(t). When t is NaN (no site of one kind, difference 8; or 0/0 in TN93 in some round,
    difference 10; or the logarithm of a negative number when one kind is saturated), every weight is NaN, the
    test `sum of weights > 1e-20` is false, and the codon's differences are silently left out. YN and MYN then
    print a Ka from the remaining codons, with Ks `nan` or `NA` beside it. Two examples:
    - `CCTCTTTCT` against `TTTCCTTTT` (no synonymous site; 4 differences at 9 sites, one codon differs at two
      positions): KaKs_Calculator's YN and MYN print Ka = 0.263548, which is JC69 of 2 differences. COAL prints
      0.673456, JC69 of all 4, and so does yn00 with `weighting = 0` (0.673456194904). In all 6 such fields of
      the 3,000 pairs COAL's value is yn00's.
    - `TGCATGATATCTTTTAATAACCGAGAAGTTTCCGAC` against `TCCATGATATCTTTTAATTCTCGAGAATTTTCCGAC`: KaKs_Calculator's
      MYN prints Ka = 0.0781 and Ks = `nan`. COAL prints Ka = 0.1673 and Ks = 0.1231, as the corrected reference
      does; YN gives 0.1672 and 0.1231 in both programs.

    **COAL printed the wrong numbers too until 2026-11-12** (0.0781 for the pair above): it followed the
    reference into the 0/0. The fix for difference 10 removed the NaN and with it the lost codons, so some MYN
    values that COAL used to print changed (13 Ka and 3 Ks among the 3,000 pairs, all upwards to the value of
    the corrected reference or, where the iteration does not converge, to no value at all (difference 12); none
    in the real genes).
12. **An iteration that does not converge: KaKs_Calculator prints its last round, COAL prints NA.** YN00
    repeats "weights from the current estimates, estimates from the weighted counts" at most 10 times, MYN at
    most 100 times. On short or nearly saturated pairs the estimates can jump between rounds for ever, for
    instance when the proportion of differences crosses the 0.75 at which JC69 stops being defined (`AAA` against
    `ATG`: Ka alternates between about 0.52 and 3.17). KaKs_Calculator and PAML yn00 then print the state of the
    last round, without a sign. That state is not a fixed point of anything, and it is not reproducible: the
    jumps amplify a difference in the 16th digit into the third digit, or across the saturation limit. For
    `AAA` / `ATG` with YN, the same COAL source printed Ka = 3.16884648849903 on Linux (GCC 15), 3.17357163988094
    on Windows (MSVC), 3.1735716240756 on macOS (Apple Clang, arm64) and "saturated" with GCC 13 and an older C
    library (2026-11-13; "Platform reproducibility" below has the counts). A value that changes in the third
    digit with the compiler is not a result, so COAL gives none: Ka, Ks, Ka/Ks, the standard errors, the
    weighted difference counts and what is derived from them (Fisher's p, Divergence-Time) are NA, in `coal
    kaks` and in `coal dnds`; the status is `not_converged`; the note says `YN00 iteration did not converge in
    10 rounds: no estimate (last round: Ka ..., Ks ...; such values can differ between platforms)` (or `MYN ...
    100 rounds`), and the summary on stderr counts such results. The site counts and kappa do not come out of the
    iteration and are printed. COAL follows the published procedure and KaKs_Calculator's code round for round
    (traced on the examples); it does not try to repair the iteration. Among the 3,000 pairs: 123 YN and 87 MYN
    rows. Among the real genes: one simulated pair (MYN; KaKs_Calculator's last round is Ka 0.58117, Ks 1.43957,
    which COAL printed until 2026-11-13), none of the orthologr pairs. Not every such iteration wanders: that
    simulated pair's and the test fixture `sim_w0.1_t2_k2_3`'s (YN, dS above 4) creep towards a value more slowly
    than the round limit allows, and their last rounds agree between platforms and with KaKs_Calculator to many
    digits. COAL does not try to tell the two kinds apart, because no test of one round's values can: the rule
    is the round limit of the references. `--yn-weighting 0` needs no iteration, and neither do the other five
    methods.
13. **A saturated kind of site changes the path of the iteration.** When every site of one kind differs (Sd = S
    or Nd = N; MYN scales the counts down to exactly that), the proportion of differences is 1 up to rounding,
    and rounding decides whether a round takes "P + Q > 1" (t = 99) or JC69's cap. In both cases that kind is
    saturated: COAL prints NA (difference 1), KaKs_Calculator the cap. The other kind's estimate goes through t
    and converges to the same value by another path, so it agrees to the tolerance of the iteration (5e-8 for
    MYN) and not to 1e-9: 5 fields among the 3,000 pairs, the largest difference 4.4e-9.

What the comparison of the 3,000 pairs with codons drawn independently gives, per field (Ka, Ks, Ka/Ks) in which
both programs print a number and the numbers differ by more than 1e-9 (relative); KaKs_Calculator run one pair
per process:

| Class | COAL against KaKs_Calculator 2.0 | COAL against the corrected reference | Who is right |
|---|---:|---:|---|
| 11, no site of one kind (COAL status `no_data`) | 6 | 6 | COAL (= yn00, `weighting = 0`) |
| 11, a NaN in the reference's row (its Ka or Ks is `nan` / `NA`) | 12 | 0 | COAL (= the corrected reference) |
| 12, iteration not converged | 0 | 0 | nothing to compare: COAL prints NA. The reference prints its last round in 493 (corrected: 578) of these fields |
| 13, saturated kind, other path | 5 | 5 | the same value within the iteration's tolerance |
| unexplained | 0 | 0 | |

So "the corrected reference gives COAL's numbers" holds for the 4,000 pairs with one changed position per codon
and for the second row of the 3,000; it says nothing where the iteration does not converge (COAL has no number
there; the corrected reference prints more last rounds than the unpatched one, because it prints numbers where
the unpatched one prints `nan`), and it does not hold for pairs without a site of one kind, which the corrected
lines do not touch. NG, LWL, MLWL, LPB and MLPB have no difference between two numbers on either random set.

**KaKs_Calculator's MYN depends on the order of the pairs in a file.** It carries state from one pair to the
next: run once over the 3,000 pairs, 49 MYN rows differ (Ka, Ks or Ka/Ks) from the same pairs run one per process
(the other six methods: none; the corrected reference: none). The comparisons of this section therefore run
KaKs_Calculator once per pair (`random_short_pairs.py` does; `--one-file` reports the count). The tables for the
real genes above were made with one run per file; there the stock and the 15-digit outputs agree with COAL to
the stated tolerances, so no such row occurs in them. COAL has no state between pairs: every pair gives the same
row alone, in any order and with any number of threads (`tests/test_kaks.py`, `test_threads_same_output`).

`tests/test_kaks.py` holds the cases (`TestOneKindOfSite`, `TestMynNan`, `TestMultiPositionCodons`), and its
comparison with KaKs_Calculator's stored output does not accept a difference because of the status alone: a
number that KaKs_Calculator prints and COAL does not is a failure unless the status is `saturated`
(difference 1) or `not_converged` (difference 12; only in the columns that come out of the iteration).

The categories, tolerances and checks are those of `tests/test_kaks.py`, `benchmarks/kaks/validate_kaks.py` and
`benchmarks/kaks/random_short_pairs.py`. On the real-gene sets (999 + 600 pairs) and on the 4,000 random pairs
with one changed position per codon, all fields agree to 1e-9 apart from differences 1 to 10. On short pairs
whose codons differ at several positions they do not: differences 11 to 13.

### Platform reproducibility

What COAL prints for a pair should not depend on the compiler or the operating system. This is what is done
about it, what was measured, and what is not guaranteed.

**Guaranteed by the source.**

- No fused multiply-add and no fast-math: the engine is compiled without floating-point contraction on every
  compiler (a `#pragma` in `kaks.cpp` for Clang, GCC and MSVC, and `-ffp-contract=off` / `/fp:precise` in
  `CMakeLists.txt`). Clang contracts `a * b + c` by default, and arm64 always has the instruction, so without
  this a macOS build and a Linux build would differ in the last digits of nearly every YN / MYN value.
- No dependence on the number of threads, on the order of the pairs, or on the locale.
- An iteration that did not converge prints no value (difference 12).

With that, builds that use the same math library print the same bytes. Measured on 2026-11-13 on x86-64 Linux
with GCC 13 (`-O2`, `-O3`, `-O3 -march=x86-64-v3`), GCC 16, Clang 18 (`-O3`, `-O3 -march=x86-64-v3`) and two
builds without any build-system flag (the pragma alone, on FMA-capable targets): the output of all eight is
byte-identical on the 999 orthologr pairs (and on a second, independent alignment of the same genes), the 600
simulated pairs, the edge cases and both random sets (about 67,000 rows; `coal kaks` at 6 and 15 digits, `coal
dnds` at 8 digits, and for all but the first set also at 15 digits and with `--yn-flavour paml`).

**Not guaranteed: the C library.** `exp`, `log` and `hypot` are not correctly rounded in every C library, and
glibc (which changed them between 2.39 and 2.43), Apple's libm and Microsoft's UCRT return a different last bit
now and then. COAL uses the platform's functions. The effect was measured by running the same binary with these
three functions replaced by ones that compute in `long double` and round:

| Rows | Differ at the default digits (6 / 8) | Differ at 15 digits |
|---|---:|---:|
| 999 orthologr pairs, 7 methods (6,993 rows) | 0 | 528 (the last one to three digits) |
| 600 simulated pairs (4,200 rows) | 0 | 483 |
| 4,000 random short pairs, one change per codon (28,000 rows) | 0 | 7 |
| 3,000 random short pairs, codons drawn independently (21,000 rows), converged | 1 (MYN: SE and Sd in the 8th digit) | 744 |
| the same set, not converged (210 rows) | 2 (and the last round quoted in the note of 6 more) | 2 |

So on real genes two platforms agree in every printed digit at the default precision, and from about the 13th
digit on they need not (`--digits 15` shows that). This is an observation on these 11,193 rows and not a proof:
a value can sit on the rounding boundary of the sixth digit. MYN stops at a tolerance of 5e-8 on dS and dN, so
its eighth digit depends on the path (the one row above). The two rows of the last line are pairs of two
codons (`CGGGCC` / `GCAAAG` and `GCATGG` / `ATGGCT`, MYN) where the C library decides whether the iteration
settles within 100 rounds: one platform prints `not_converged`, the other `saturated`. No such row was found
among the real genes.

**Before 2026-11-13** the last round of an iteration that did not converge was printed. The same measurement
then: with FMA contraction 11 to 15 of the 210 such rows differed at the default digits and 197 to 203 at 15
digits (in 4 or 5 the convergence itself); with the other math library 8 to 10 and 156 to 161. And 2 to 4 of the
20,790 converged rows of that random set differed at the default digits with FMA: Fisher's p of a table with a
cell that is zero up to rounding (0 or NA; such a cell now counts as zero), the difference counts of a one-codon
pair that is saturated in both kinds, and eighth digits of MYN. Without FMA and with the same C library every
build was identical then too (GCC `-O0`, `-O2`, `-O3`, with AVX2; GCC 16; Clang `-O0`, `-O3`).

**Speed.** Giving up FMA costs nothing measurable: on the 999 orthologr pairs (YN + MYN, one thread, best of 5)
a GCC 13 `-march=x86-64-v3` build took 2.41 s with contraction and 2.44 s without, a Clang 18 build 3.54 s and
3.52 s. Default x86-64 builds never had the instruction, and their code is unchanged.

**Tests.** `tests/test_kaks.py` asserts for the pairs that do not converge only what holds on every platform
(status, NA fields, the form of the note), never a digit of a last round.

### Where KaKs_Calculator and PAML differ (YN00)

On the 999 orthologr pairs and 600 simulated pairs, KaKs_Calculator's YN and PAML yn00 with `weighting = 1` give
the same dN, dS and dN/dS (COAL matches both to about 1e-13). Both start the iteration at t = 0.4 and omega = 1
and use the same tolerance (5e-4, 10 rounds). They differ only in edge cases, which `--yn-flavour` selects:

- **kappa when none can be estimated.** KaKs_Calculator uses 2. PAML uses 4.6, the initial value of its main
  loop.
- **kappa above 99.** KaKs_Calculator resets it to 2. PAML caps it at 999.
- **dN/dS when dN = 0.** PAML prints its floor `omega = 1e-5`; COAL prints 0. This is the 1e-5 maximum difference
  in the dN/dS row of the PAML comparison.
- **Pathway weighting.** PAML's distributed `yn00.ctl` has `weighting = 0` (equal weights, one round).
  KaKs_Calculator always weights. COAL's default follows KaKs_Calculator; `--yn-weighting 0` matches PAML's
  default.

### Recovering a known omega (simulated data)

The table gives the median dN/dS over 50 pairs of 500 codons. Within each cell, the three values are for
t = 0.3, 1 and 2 substitutions per codon.

| Method | omega 0.1 | omega 0.5 | omega 1 | omega 2 | NA (of 600) |
|---|---|---|---|---|---:|
| YN | 0.11 / 0.08 / 0.09 | 0.53 / 0.49 / 0.48 | 1.06 / 1.03 / 1.04 | 2.02 / 2.01 / 1.95 | 1 |
| MYN | 0.10 / 0.08 / 0.10 | 0.53 / 0.49 / 0.46 | 1.07 / 1.04 / 1.06 | 2.14 / 2.09 / 2.05 | 1 |
| LPB | 0.12 / 0.12 / 0.13 | 0.53 / 0.51 / 0.54 | 1.05 / 1.01 / 0.93 | 1.92 / 1.72 / 1.61 | 2 |
| MLPB | 0.11 / 0.11 / 0.12 | 0.51 / 0.49 / 0.51 | 1.02 / 0.97 / 0.92 | 1.81 / 1.69 / 1.55 | 2 |
| NG | 0.09 / 0.08 / 0.07 | 0.46 / 0.45 / 0.45 | 0.91 / 0.89 / 0.86 | 1.71 / 1.64 / 1.49 | 27 |
| LWL | 0.10 / 0.10 / 0.11 | 0.44 / 0.43 / 0.43 | 0.85 / 0.86 / 0.78 | 1.58 / 1.44 / 1.37 | 2 |
| MLWL | 0.11 / 0.10 / 0.10 | 0.46 / 0.44 / 0.44 | 0.87 / 0.86 / 0.78 | 1.60 / 1.45 / 1.31 | 2 |

YN00 and MYN recover omega to within about 10 % at every divergence; at omega = 0.1 and t = 1 they run about 20 %
low. The methods that ignore the transition/transversion bias (NG, LWL) underestimate omega, increasingly so at
high divergence and high omega. This is the bias Yang & Nielsen (2000) describe. NG's 27 NA are Ks saturations
(Jukes-Cantor) at omega = 0.1, t = 2.

## Benchmark

The benchmark used the 999 orthologr pairs on the cloud workspace, a 2-core Xeon VM at 2.1 GHz that other jobs
were also using (load average about 5). Wall times there are noisy, so CPU time (user + sys) is the fair
single-thread measure. Thread scaling could not be measured on that machine, and the headline numbers are to be
taken on a 16-thread desktop. KaKs_Calculator is single-threaded. Each figure is the best of 3 runs
(`benchmarks/kaks/bench_kaks.py`).

| Method | KaKs_Calculator CPU s | coal -j 1 CPU s | Pairs per CPU-second, KaKs_Calculator | Pairs per CPU-second, COAL | CPU speed-up |
|---|---:|---:|---:|---:|---:|
| YN | 1.88 | 1.15 | 530 | 871 | 1.6x |
| MYN | 2.38 | 1.71 | 420 | 585 | 1.4x |
| NG | 0.53 | 0.07 | 1,877 | 13,787 | 7.3x |
| LWL | 0.76 | 0.09 | 1,321 | 11,249 | 8.5x |
| MLWL | 0.70 | 0.09 | 1,431 | 11,706 | 8.2x |
| LPB | 0.67 | 0.09 | 1,484 | 11,654 | 7.9x |
| MLPB | 0.81 | 0.09 | 1,236 | 11,435 | 9.2x |

YN00 and MYN spend about 90 % of their time on P(t) (a 61 × 61 symmetric eigendecomposition per iteration, as in
both reference tools), so their single-thread gain is modest. The gain comes from threads (`-j`): pairs are
independent and computed on the thread pool. In a pipeline, Ka/Ks takes a small share of the time (in the orthologr 0.4.3
profile of 2026-10-29, about 0.5 %); alignment dominates.

## Not done (yet)

- GY94 maximum likelihood, model selection and averaging (MS, MA), and the gamma variants.
- Standard input.
- A Python API: the commands are reachable from `coal`, and an in-process binding is planned with the R/Python
  bindings work.
- KaKs_Calculator 3.0 was not run (see Output).
- The textbook NG86 site count (see the note under "Methods"): `-m NG` is KaKs_Calculator's.
- A ThreadSanitizer job in the nightly checks: the race on `signgam` was found and its fix confirmed by hand
  (`coal-cpp kaks -j 4` in a `-fsanitize=thread` build: one report before, none after).
