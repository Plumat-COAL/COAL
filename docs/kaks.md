# Ka/Ks (dN/dS): `coal kaks` and `coal dnds`

COAL estimates Ka (dN) and Ks (dS) for pairs of aligned coding sequences in process, with a C++ engine
(`cpp/src/kaks.cpp`, `kaks_comeron.cpp`, `kaks_li.cpp`) derived from KaKs_Calculator 2.0, PAML's yn00,
libsequence and seqinr (see "Where the code comes from" below). No external program runs. It has two commands:

- **`coal kaks`** does the job of `KaKs_Calculator -i IN.axt -o OUT.kaks -m METHOD` for seven of KaKs_Calculator's
  methods (NG, LWL, MLWL, LPB, MLPB, YN, MYN). It takes KaKs_Calculator's `-i`, `-o`, `-m` and `-c` and writes its
  tab-separated format: the same 22 columns, the same NA rules and 6 significant digits. It is not a drop-in for
  every KaKs_Calculator command line: the other methods are refused, and without `-m` it computes YN, where
  KaKs_Calculator computes MA (model averaging). It also has two methods that KaKs_Calculator does not have,
  `Comeron` and `Li`, the two that the R package orthologr computes without KaKs_Calculator (see
  "Comeron and Li" below).
- **`coal dnds`** writes a tidy TSV from the same engine, with one row per pair and method, plus standard errors,
  a status and a note.

```
coal kaks -i pairs.axt -o pairs.kaks                     # YN00 (the default)
coal kaks -i pairs.axt -o pairs.kaks -m NG -m YN -m MYN  # several methods; rows come in KaKs_Calculator's order
coal dnds codon_aln.fasta -o dnds.tsv                    # every pair of an aligned FASTA (MSA), YN00
coal dnds pairs.fasta --pairs consecutive -m YN -m NG -j 8 -c 2
coal dnds codon_aln.fasta -m Comeron -m Li               # Comeron (1995) and Li (1993), COAL's policy
coal dnds codon_aln.fasta -m Comeron --comeron-flavour gestimator   # the numbers of orthologr's gestimator
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
| `Comeron` (alias `C95`) | Comeron (1995) as libsequence's `gestimator` computes it (the default method of the R package orthologr): non-degenerate, transitional and transversional 2-fold and 4-fold sites, Kimura's two-parameter correction per class, pathways weighted by Grantham's distance | J Mol Evol 41:1152-1159 | NA |
| `Li` | Li (1993) as `seqinr::kaks()` computes it (orthologr's `"Li"`): 0-, 2- and 4-fold sites, pathways weighted by an amino-acid similarity table | J Mol Evol 36:96-99 | yes (seqinr's variance formulas) |

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

`Comeron` and `Li` are not KaKs_Calculator methods. Their names are orthologr's, they are read in any case
(`-m comeron`), their rows come after KaKs_Calculator's seven, and they are implemented for the standard genetic
code only, as their references are (`-c` other than 1 with them is refused with exit status 2). `Li` is not
`LPB`: the two implement the same formulas of Li (1993) with different site and pathway counting, and give
different numbers (see "Comeron and Li").

**Where the code comes from.** The engine is not an independent implementation: `cpp/src/kaks.cpp` is derived
from the source code of KaKs_Calculator 2.0 (Copyright (C) 2009 Beijing Institute of Genomics of CAS; Da-Peng Wang,
Yu-Bin Zhang, Zhang Zhang, Kristian K Ullrich) and of PAML's `yn00.c` (Copyright 1998 Ziheng Yang), both under the
GNU GPL version 3. Its functions follow theirs in structure, constants and order of operations, which is why the
numbers are the same; they were rewritten without classes and global state, for any of COAL's genetic codes and
for threads, and changed where listed under "Differences from KaKs_Calculator". The header of `kaks.cpp` and the
`NOTICE` file give the details. This is the reason COAL is distributed under GPL-3.0-only from this release on.
Neither program is included as such, and `coal kaks` / `coal dnds` do not run them.

In the same way, `cpp/src/kaks_comeron.cpp` is derived from libsequence (Copyright (C) 2003-2009 Kevin Thornton;
GNU GPL version 3 or later) in the form the R package orthologr compiles into its `gestimator()` (modified by
Sarah Scharfenberg and Hajk-Georg Drost, 2014; their changes under GPL version 2 or later), and
`cpp/src/kaks_li.cpp` from `src/kaks.c` of the R package seqinr 4.2-44 (licence `GPL (>= 2)`; the file has no
notice of its own, and seqinr's help says the C program "was provided by Manolo Gouy"). Both are modified
versions of that code, with the same structure, constants and order of floating-point operations, which is why
their numbers are the references' to the last bit. The headers of the two files and `NOTICE` give the exact
licence statements and what was changed. All three licences allow distribution under GPL version 3.

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

For `Comeron` and `Li` the same 22 columns are written with what the method defines: Ka, Ks, Ka/Ks, Length,
Substitutions (nucleotide differences in the compared codons) and GC; for `Li` also `Fold-Sites(0:2:4)` (its
L0:L2:L4) and the two `Fold-*-Substitutions` columns (its transitional and transversional differences per
class). The other columns are `NA` (there are no S / N counts, so no Fisher test), and the rate-ratio column is
`1:1:1:1:1:1`, as for NG. The Method column says `Comeron` or `Li`.

Rows come per pair in KaKs_Calculator's method order (NG, LWL, MLWL, LPB, MLPB, YN, MYN), then Comeron, then
Li. Pairs with status
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
| `method` | NG, LWL, MLWL, LPB, MLPB, YN, MYN, Comeron or Li |
| `dN`, `dS`, `dNdS` | Ka, Ks, Ka/Ks. `NA` when not available; dNdS is `NA` when dS < 1e-6. Zeros print as 0 |
| `SE_dN`, `SE_dS` | standard errors (YN, MYN, NG, Li); `NA` elsewhere |
| `N_sites`, `S_sites` | nonsynonymous and synonymous sites (`NA` for Comeron and Li, which count sites by degeneracy class, not as S and N) |
| `Nd`, `Sd` | nonsynonymous and synonymous differences (`NA` for Comeron and Li) |
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

`--comeron-flavour gestimator` and `--li-flavour seqinr` make `Comeron` and `Li` read the pair as their
reference programs do ("Comeron and Li" below); the default for both is `coal`, the policy above.

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

`Options{yn_flavour, yn_weighting, comeron_flavour, li_flavour}` is the fourth argument (the defaults follow
KaKs_Calculator, and COAL's policy for Comeron and Li: `ComeronFlavour::coal`, `LiFlavour::coal`; the reference
flavours are `ComeronFlavour::gestimator` and `LiFlavour::seqinr`). `Method::Comeron` and `Method::Li` are the
two further methods; `parse_method()` reads `"Comeron"`, `"C95"` and `"Li"` in any case; with a genetic code
other than the standard one they return status `error`. Their degeneracy tables are built once, into immutable
function-local statics (thread-safe initialisation), so the statement about global state holds for them too.
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
  compiler (a `#pragma` in `kaks.cpp`, `kaks_comeron.cpp` and `kaks_li.cpp` for Clang, GCC and MSVC, and
  `-ffp-contract=off` / `/fp:precise` in `CMakeLists.txt`; coalr's copies of the three files carry the pragma,
  and so does its `coalr.cpp`, which computes the quartiles of the cluster summary). Clang contracts `a * b + c` by default, and arm64 always has the instruction, so without
  this a macOS build and a Linux build would differ in the last digits of nearly every YN / MYN value.
- No dependence on the number of threads, on the order of the pairs, or on the locale.
- An iteration that did not converge prints no value (difference 12).

With that, builds that use the same math library print the same bytes. Measured on 2026-11-13 on x86-64 Linux
with GCC 13 (`-O2`, `-O3`, `-O3 -march=x86-64-v3`), GCC 16, Clang 18 (`-O3`, `-O3 -march=x86-64-v3`) and two
builds without any build-system flag (the pragma alone, on FMA-capable targets): the output of all eight is
byte-identical on the 999 orthologr pairs (and on a second, independent alignment of the same genes), the 600
simulated pairs, the edge cases and both random sets (about 67,000 rows; `coal kaks` at 6 and 15 digits, `coal
dnds` at 8 digits, and for all but the first set also at 15 digits and with `--yn-flavour paml`).

**Comeron and Li** were measured apart (2026-11-13, x86-64 Linux, glibc 2.39): GCC 13.3 `-O3`, `-O3
-march=x86-64-v3` and `-O0`, Clang 18.1 `-O3` and `-O3 -march=x86-64-v3`, and two builds without any build-system
flag on an FMA target (GCC `-O2 -march=x86-64-v3 -ffp-contract=fast`, Clang `-O2 -march=x86-64-v3`: the pragma
alone, as R compiles coalr). On the 14,153 pairs of "Validation of Comeron and Li", both methods, both flavours,
`coal dnds` at 17 and at 8 digits and `coal kaks` at 15 digits (54 outputs, 169,828 rows per build), the seven
builds print the same bytes: 0 rows differ. So the statement of that section (every number is the same double as
gestimator's and `seqinr::kaks`'s) holds for builds that target FMA as well, on this C library. It did not before
the pragma: the same two FMA builds with the pragma removed differ from the default build in 3,045 of the 28,306
rows of the reference flavours at 17 digits (534 Comeron, 2,511 Li; in the last one or two digits) and in none at
8 digits. Both methods call `log` and `sqrt`, so the next paragraph applies to them too: with another C library
the last bit of a value can differ from the reference's, which is why the tests compare to a relative 1e-12 and
not bit for bit.

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

## Comeron and Li

`-m Comeron` and `-m Li` are the two dN/dS methods that the R package
[orthologr](https://github.com/drostlab/orthologr) computes without KaKs_Calculator: `"Comeron"`, its default,
with libsequence's `gestimator` (compiled into the package), and `"Li"` with `seqinr::kaks()`. COAL has both so
that an orthologr run can be reproduced, or continued, with COAL's engine. They are ports of those two
programs (see "Where the code comes from"), not new implementations of the papers.

**Comeron (1995).** Every codon position is put in one of four classes: non-degenerate (L0), 2-fold degenerate
where the synonymous change is the transition (L2S) or the transversion (L2V), and 4-fold degenerate (L4).
Transitions (P) and transversions (Q) are counted per class; for codons that differ at two or three positions
the two or six pathways are weighted by the inverse of their length in Grantham's (1974) amino-acid distance. With
Kimura's two-parameter correction per class, A = -1/2 ln(1 - 2P - Q) + 1/4 ln(1 - 2Q) and B = -1/2 ln(1 - 2Q):

    Ks = (L2S A2S + L4 A4) / (L2S + L4) + Bs,   Bs from Qs = (Q2V + Q4) / (L2V + L4)
    Ka = (L2V A2V + L0 A0) / (L2V + L0) + Ba,   Ba from Qa = (Q0 + Q2S) / (L0 + L2S)

There is no variance, so `SE_dN` and `SE_dS` are `NA`.

**Li (1993).** Positions are non-degenerate (L0), 2-fold (L2) or 4-fold (L4); transitional (A) and
transversional (B) distances per class, again Kimura's; Ks = (L2 A2 + L4 A4) / (L2 + L4) + B4 and
Ka = A0 + (L0 B0 + L2 B2) / (L0 + L2). seqinr's counting has rules of its own: the three isoleucine codons and
the arginine codons are 2-fold with corrections for which changes are synonymous, and the pathways between
codons that differ at several positions are weighted with a table of amino-acid similarities (a stop codon gets
the smallest value of the table). The standard errors are the square roots of seqinr's variances `vka` and `vks`.

**`Li` is not `LPB`.** Both are "Li (1993)", and KaKs_Calculator's LPB / MLPB use the same two formulas, but the
site classes and the pathway weights differ. On the 999 orthologr pairs Li and LPB agree to 1e-6 in 105 (dN)
and 102 (dS) pairs, Li and MLPB in 273 and 263; the largest differences are 0.034 (dN) and 0.069 (dS) against
LPB, 0.018 and 0.053 against MLPB, with correlations of 0.999 and 0.995 or better. So `Li` is a method of its
own, not an alias (`tests/test_kaks_comeron_li.py`, `test_li_is_not_lpb`).

### Two flavours

Each method can read a pair in two ways. The default (`coal`) applies the policy of this page, as for every
other method. The other (`--comeron-flavour gestimator`, `--li-flavour seqinr`) reads the pair as the reference
program does, to reproduce numbers that were made with orthologr. The counting and the formulas are the same
in both; on a pair without stop codons, without `U`, and with every correction defined, the two flavours give
the same doubles.

| | `Comeron`, default (`coal`) | `--comeron-flavour gestimator` (orthologr 0.4.3's gestimator) |
|---|---|---|
| Stop codon at the end of both rows | the column is dropped | compared like any codon: a stop codon has no sites (the column counts half the sites of the other codon), `TAA` against `TGA` or against a sense codon is counted as a change, and a pathway through a stop codon is kept, with Grantham distance 999 |
| Stop codon elsewhere | the column is skipped; status `warning` | compared, as above; no warning |
| `U` | read as `T` | the codon is not A/C/G/T: the column is skipped |
| Gap or ambiguous letter in a codon | column skipped | column skipped (the same) |
| A logarithm of the correction has an argument <= 0, in a class that has sites | the estimate is `NA`; status `saturated` | gestimator sets that logarithm to 0 and prints the finite number that results (too small). COAL gives that number, with status `warning` and the note `Ka: correction undefined, gestimator's value (...)` (or `Ks:`) |
| 1 - 2Qa <= 0 (Ba undefined) | Ka is `NA`; `saturated` | gestimator's rule: Ba = 0 when Kimura's distance over the whole rows is below 1 (a number comes out, flagged as above); otherwise it prints 999: `NA`, `saturated` |
| 1 - 2Qs <= 0 (Bs undefined) | Ks is `NA`; `saturated` | gestimator prints 999: `NA`, `saturated` |
| No synonymous (nonsynonymous) sites | `NA`; `no_data` | gestimator prints 999: `NA`, `no_data` |
| dN/dS | `NA` when dS < 1e-6, as for the other methods | gestimator's: 999 (`NA`) only when dS is 0, or dN or dS is 999 |
| Rows of different length, or not a multiple of 3 | status `error` | status `error` (gestimator 0.4.3 reads beyond the rows) |

gestimator's "cannot be calculated" value 999 is never a number in COAL; orthologr turns it into `NA` too. One
refinement of the default flavour: an undefined logarithm in a class that has **no** sites (weight 0 in the
formula) is not a saturation.

| | `Li`, default (`coal`) | `--li-flavour seqinr` (`seqinr::kaks`) |
|---|---|---|
| Stop codons (at the end, or elsewhere) | dropped; skipped with status `warning` | compared like any codon: a stop is an "amino acid" with the smallest similarity, and its positions count as sites (so even two identical stop codons change L0, L2 and with them Ka and Ks a little) |
| `U` | read as `T` | the codon is removed with its column |
| Gap or ambiguous letter | column skipped | column removed (the same for two sequences; for an alignment of several sequences `seqinr::kaks` removes a column when **any** sequence has a gap, COAL only when one of the two has) |
| Correction undefined | `NA`; `saturated` | seqinr returns 9.999999 for the value and its variance. COAL: `NA`, `saturated`, in both flavours |
| A class of sites is empty (Ks needs L2 and L4, Ka needs L0 and L2) | `NA`; `no_data` | seqinr returns 9.999999 here too ("saturated"); COAL: `NA`, `no_data`, in both flavours |
| Negative Ka or Ks (the formulas give one for some short pairs) | the number, status `warning` | the number, status `warning` (seqinr returns it; orthologr turns it into `NA`) |

`coal dnds` and `coal kaks` compare the rows they are given, so with a reference flavour a terminal stop codon
is compared when the rows contain it. `coal pipeline`, `coal codon-dnds` and coalr's `coal_dnds_pairs()` /
`coal_dnds_clusters()` make the codon alignment themselves, and it ends before a terminal stop codon. Their
option `--terminal-stops compare` (`terminal_stops = "compare"`) adds the two stop codons as the last compared
column for the reference flavours when both genes end in one. orthologr 0.4.3 needs it: it translates the stop
codon to `*`, aligns it, and PAL2NAL threads it, so its gestimator compares the two stop codons (the fork of
orthologr with COAL built in translates without the stop codon; for its numbers the option stays `drop`).

**A property of libsequence that both flavours keep.** For codons that differ at two positions, libsequence's
`TwoSubs` passes (codon 1, intermediate 1, codon 2, intermediate 2) to a function whose parameters are
(codon 1, intermediate 1, intermediate 2, codon 2). The four branches counted are therefore codon 1 to
intermediate 1, intermediate 1 to intermediate 2, codon 1 to codon 2, and codon 2 to intermediate 2, where
the two pathways are codon 1 to intermediate 1 to codon 2 and codon 1 to intermediate 2 to codon 2; and
`SingleSub`, given two codons that differ at two positions, classifies the change at the last of them. The
same nucleotide changes are counted, but the degeneracy class of some is taken from another codon than the
pathway's. The current libsequence (1.9.8) makes the same call, so this is how "Comeron as gestimator" has
always been computed, and COAL reproduces it in both flavours because the purpose of the method is to give
those numbers. It is one reason why Comeron and the other methods differ by a few percent (cross-checks below).
Use `YN` when the estimate itself matters more than comparability with orthologr.

### Validation of Comeron and Li

References (`benchmarks/kaks/ref/`): **gestimator** is orthologr 0.4.3's own copy of the libsequence headers,
compiled unmodified with a driver that calls them as `gestimator()` does (`maxHits = 3`) and prints 17
significant digits (`gestimator_ref.cpp`; the package function prints 4). **seqinr** is `seqinr::kaks()` 4.2-44
called per pair from R, 17 digits (`seqinr_ref.R`). `benchmarks/kaks/validate_comeron_li.py` runs the
comparisons; a value is "identical" when the 17-digit texts give the same double. Results of 2026-11-13,
Linux x86-64, gcc 13.3:

| Set | Pairs | Comeron `gestimator`: dN, dS identical to gestimator | 999 in gestimator = NA in COAL (dN, dS) | Li `seqinr`: Ka, Ks identical to seqinr | 9.999999 in seqinr = NA in COAL (Ka, Ks) | Mismatches |
|---|---:|---:|---:|---:|---:|---:|
| orthologr 999 (COAL's alignments, no stop codons) | 999 | 999, 998 | 0, 1 | 999, 998 | 0, 1 | 0 |
| orthologr 0.4.3's own 984 codon alignments (terminal stop codons in 968) | 984 | 984, 981 | 0, 3 | 984, 981 | 0, 3 | 0 |
| simulated 600 | 600 | 600, 600 | 0, 0 | 600, 599 | 0, 1 | 0 |
| `tests/kaks/comeron_li_cases.axt` (hand-made edge cases, stops, gaps, IUPAC, case, U, 1 to 3 codons, saturated) | 499 | 471, 408 | 27, 90 | 376, 254 | 122, 244 | 0 |
| `edge_cases.axt`, `sim_small.axt` | 71 | 69, 64 | 1, 6 | 65, 62 | 5, 8 | 0 |
| random short pairs, one changed position per codon | 4,000 | 4,000, 3,810 | 0, 190 | 3,512, 2,725 | 488, 1,275 | 0 |
| random short pairs, codons drawn independently | 3,000 | 2,866, 2,232 | 134, 768 | 1,773, 1,005 | 1,227, 1,995 | 0 |
| random short pairs from all 64 codons, with stops, gaps, N / IUPAC, lower case, U | 4,000 | 3,724, 3,179 | 276, 821 | 3,118, 2,055 | 882, 1,945 | 0 |

Every number that a reference returns is the same double in COAL (no value differs at all, so the 1e-12
tolerance of the tests is not used on this platform), dN/dS and the standard errors (against the square roots of
seqinr's variances) included; every 999 and 9.999999 is `NA`; the two pairs that both references refuse (a
length that is not a multiple of 3) have status `error`. gestimator's Kimura rule (Ba = 0) decided 248 of the
random and hand-made rows. The two files are compiled without floating-point contraction, so a build for a
target with FMA (`-march=x86-64-v3`, arm64) gives the same bytes; a C library with another `log` can differ in
the last bit ("Platform reproducibility" has the counts for both).

**What the default flavours change** (default against reference flavour, same input):

| Set | Method | Pairs the policy reads differently | dN changes | by > 1% | by > 20% | number to `NA` | dS changes | by > 1% | by > 20% | number to `NA` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| orthologr 999 (no stop codons) | Comeron | 0 | 1 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
| orthologr 999 | Li | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| orthologr 0.4.3's own 984 alignments | Comeron | 968 | 111 | 105 | 8 | 3 | 1 | 1 | 0 | 0 |
| orthologr 0.4.3's own 984 alignments | Li | 968 | 957 | 124 | 5 | 0 | 572 | 39 | 0 | 0 |
| simulated 600 | Comeron | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 2 |
| simulated 600 | Li | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

So on real genes as orthologr 0.4.3 aligns them, dropping the terminal stop codons changes Comeron's dN in the
108 pairs whose two stop codons differ (gestimator counts `TAA` against `TGA` as a substitution), by more than
1% in 105 and by more than 20% in 8 (0.8% of the pairs); 3 further values are `NA` where gestimator computes a
number from an undefined correction; dS changes in one pair. Li's values change a little in nearly all pairs, because seqinr
counts the positions of stop codons as sites. In the three tables, every change is explained by the policy
(a stop codon or `U` in a compared column) or by the saturation rule; the script fails otherwise ("unexplained":
0 in all sets, the random ones included).

**Cross-checks** (default flavours; COAL's MLPB, LPB and YN are KaKs_Calculator's numbers):

| Set | Comeron vs MLPB: r (dN, dS), median ratio (dN, dS) | Comeron vs YN | Comeron vs Li | Li vs MLPB | Li vs YN |
|---|---|---|---|---|---|
| orthologr 999 | 0.9997, 0.9920; 0.988, 1.027 | 0.9989, 0.9304; 0.992, 1.006 | 0.9998, 0.9931; 0.994, 1.016 | 0.9997, 0.9973; 1.000, 1.001 | 0.9991, 0.9226; 0.998, 0.988 |
| orthologr 0.4.3's 984 | 0.9997, 0.9924; 0.988, 1.027 | 0.9990, 0.9410; 0.992, 1.003 | 0.9998, 0.9930; 0.994, 1.016 | 0.9995, 0.9974; 1.000, 1.001 | 0.9991, 0.9402; 0.999, 0.988 |
| simulated 600 | 0.9999, 0.9974; 0.980, 1.064 | 0.9994, 0.9550; 0.966, 1.087 | 0.9999, 0.9968; 0.988, 1.025 | 0.9999, 0.9993; 0.993, 1.030 | 0.9995, 0.9478; 0.979, 1.046 |

(The audit of orthologr of 2026-11-11 measured Comeron against KaKs_Calculator's MLPB on its 984 pairs:
r = 0.9998 for dN and 0.9926 for dS.)

### With orthologr

`benchmarks/kaks/orthologr_e2e.R` runs `orthologr::dNdS(dnds_est.method = "Comeron")` with its default path
(reciprocal best hits, pairwise protein alignment, PAL2NAL, gestimator), keeps orthologr's protein and codon
alignment of every pair, and gives the same pairs to coalr. gestimator prints 4 significant digits, so "equal"
means: COAL's value printed with `%.4g` is orthologr's. orthologr 0.4.3 (CRAN / GitHub master), its examples of
20 and of 1,000 genes per species (*A. thaliana*, *A. lyrata*), one core, 2026-11-13:

| What coalr is given | 20 genes: pairs equal to orthologr (dN, dS, dN/dS) | 1,000 genes |
|---|---:|---:|
| orthologr's codon alignments: `coal_kaks(method = "Comeron", comeron_flavour = "gestimator")` | 20, 20, 20 of 20 | 984, 984, 984 of 984 |
| orthologr's protein alignments and the CDS: `coal_codon_align()`, then `coal_kaks()` as above (the codon alignments are PAL2NAL's in all pairs) | 20, 20, 20 of 20 | 984, 984, 984 of 984 |
| the CDS only: `coal_dnds_pairs(methods = "Comeron", comeron_flavour = "gestimator", terminal_stops = "compare", preset = "pwalign")`, pairs where COAL's protein alignment is orthologr's | 7, 7, 7 of 7 | 444, 444, 444 of 444 |
| the same call, pairs where the protein alignments differ | 4, 2, 3 of 13 | 204, 144, 135 of 540 |

So the engine gives orthologr's numbers for every pair when the alignment is the same, and the protein
alignment is the only source of difference. **It differs in more than half of the pairs, and not because of
ties:** orthologr 0.4.3 calls `pwalign::pairwiseAlignment(type = "global")` without a substitution matrix, which
scores identities and mismatches only; COAL aligns with BLOSUM62 (or BLOSUM45 / 80, PAM250). Among the 540 pairs,
409 have another dN or dS at 4 digits: the median difference is 2.7% (dN) and 3.3% (dS), above 1% in 267 (dN) and
331 (dS) pairs. Over all 981 pairs with values, r = 0.960 (dN) and 0.925 (dS). To reproduce a published orthologr
0.4.3 table exactly, give COAL orthologr's protein alignments (second row of the table).

The fork of orthologr with COAL built in (`integration/all`, 0.5.0) aligns with BLOSUM62 and translates without
the stop codon. Against its default path (995 pairs; `options(orthologr.legacy_comeron_saturation = TRUE)` for
the 0.4.3 treatment of saturation; `terminal_stops = "drop"`): from its protein alignments, 995 of 995 pairs
equal in dN, dS and dN/dS; from the CDS, COAL's protein alignment is the fork's in 505 pairs, all 505 equal;
in the other 490 pairs, 291 (dN) and 185 (dS) are equal, the median difference of the rest is 0.07% (dN) and
2.8% (dS), and r over all pairs is 0.9994 (dN) and 0.9943 (dS). (Why these 490 alignments differ, co-optimal
alignments or the scoring of end gaps, was not examined here.)

**Time**, 1,000 genes, one core, on a 2-core cloud VM that other jobs were using: orthologr 0.4.3's default path
1,176 s in this run (745 s in a quieter run of 2026-11-11, one file and several processes per pair);
`coal_dnds_pairs()` with Comeron 0.57 s for the same 984 pairs from the coding sequences (alignment, codon
threading and Ka/Ks); `coal_codon_align()` + `coal_kaks()` on orthologr's protein alignments 0.25 s. The
reciprocal-best-hit search is not in COAL's figure.

**What an orthologr build with `engine = "coal"` has to do to offer Comeron and Li** (the fork's
`dNdS_coal()` calls `coalr::coal_dnds_pairs()`; nothing in orthologr was changed for this work):

1. accept `dnds_est.method = "Comeron"` and `"Li"` for the COAL engine when
   `all(c("Comeron", "Li") %in% coalr::coal_info()$methods)` (its `coal_methods()` lists KaKs_Calculator's seven);
2. pass `comeron_flavour = "gestimator"` and `li_flavour = "seqinr"` to `coalr::coal_dnds_pairs()` (and to
   `coal_dnds_clusters()`, `coal_kaks()`) for orthologr's numbers, with `terminal_stops = "compare"` only to
   reproduce 0.4.3, whose alignments include the stop codons; or leave the defaults for COAL's policy and say so;
3. saturation: with `comeron_flavour = "gestimator"` a row whose note contains `Ka: correction undefined` (or
   `Ks:`) holds the finite value of orthologr 0.4.3. The fork returns `NA` for these unless
   `options(orthologr.legacy_comeron_saturation = TRUE)`: set dN (dS) and dN/dS to `NA` for such rows otherwise;
4. Li: coalr already returns `NA` where seqinr returns 9.999999, and flags a negative value with status
   `warning` (orthologr turns it into `NA`);
5. values are doubles, where gestimator's file has 4 significant digits: `signif(x, 4)` gives orthologr's
   printed values.

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

`Comeron` and `Li` are table look-ups per codon column: 0.08 CPU s each for the 999 pairs on the same machine
(2026-11-13; about 13,000 pairs per CPU-second, the same as NG and the LWL family; gestimator's code compiled
on its own takes 0.09 s). Reading the pair for the two reference flavours (the columns with stop codons, the
counts of Kimura's distance) is part of preparing every pair and adds about 0.01 s per 1,000 pairs to the
other methods.

YN00 and MYN spend about 90 % of their time on P(t) (a 61 × 61 symmetric eigendecomposition per iteration, as in
both reference tools), so their single-thread gain is modest. The gain comes from threads (`-j`): pairs are
independent and computed on the thread pool. In a pipeline, Ka/Ks takes a small share of the time (in the orthologr 0.4.3
profile of 2026-10-29, about 0.5 %); alignment dominates.

## Not done (yet)

- GY94 maximum likelihood, model selection and averaging (MS, MA), and the gamma variants.
- Comeron and Li for genetic codes other than the standard one (their references have no other; Comeron's
  tables could be built for any code, seqinr's site classes are written out for the standard code).
- Comeron with `maxHits` 1 or 2 (libsequence's option to leave out codons that differ at several positions;
  orthologr always uses 3).
- A corrected two-substitution count for Comeron (see "Comeron and Li": both flavours keep libsequence's).
- Standard input.
- A Python API: the commands are reachable from `coal`, and an in-process binding is planned with the R/Python
  bindings work.
- KaKs_Calculator 3.0 was not run (see Output).
- The textbook NG86 site count (see the note under "Methods"): `-m NG` is KaKs_Calculator's.
- A ThreadSanitizer job in the nightly checks: the race on `signgam` was found and its fix confirmed by hand
  (`coal-cpp kaks -j 4` in a `-fsanitize=thread` build: one report before, none after).
