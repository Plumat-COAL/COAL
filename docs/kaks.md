# Ka/Ks (dN/dS): `coal kaks` and `coal dnds`

COAL estimates Ka (dN) and Ks (dS) for pairs of aligned coding sequences in process, with a C++ engine
(`cpp/src/kaks.cpp`, `kaks_comeron.cpp`, `kaks_li.cpp`) derived from KaKs_Calculator 2.0, PAML's yn00,
libsequence and seqinr (see "Where the code comes from" below), and with a maximum-likelihood method of its own
(`kaks_gy.cpp`: `GY`, the codon model of Goldman & Yang 1994 as PAML's codeml fits it to a pair). No external
program runs. It has two commands:

- **`coal kaks`** does the job of `KaKs_Calculator -i IN.axt -o OUT.kaks -m METHOD` for eight of KaKs_Calculator's
  methods (NG, LWL, MLWL, LPB, MLPB, GY, YN, MYN). It takes KaKs_Calculator's `-i`, `-o`, `-m` and `-c` and writes its
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
coal dnds codon_aln.fasta -m GY --ml-columns             # maximum likelihood (GY94), with kappa, t and lnL
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
| `GY` (aliases `GY-HKY`, `GY94`) | Goldman & Yang (1994) by maximum likelihood: the branch length t, kappa and omega are estimated from the codon pairs with F3x4 codon frequencies, as PAML's codeml does for a pair (`runmode = -2`); dN and dS as Yang & Nielsen (2000) define them. Any genetic code. Slower than YN, and the method that stays accurate at high dS (see "GY: maximum likelihood") | Mol Biol Evol 11:725-736 | yes (observed information and the delta method, as codeml's `getSE`) |
| `Comeron` (alias `C95`) | Comeron (1995) as libsequence's `gestimator` computes it (the default method of the R package orthologr): non-degenerate, transitional and transversional 2-fold and 4-fold sites, Kimura's two-parameter correction per class, pathways weighted by Grantham's distance | J Mol Evol 41:1152-1159 | NA |
| `Li` | Li (1993) as `seqinr::kaks()` computes it (orthologr's `"Li"`): 0-, 2- and 4-fold sites, pathways weighted by an amino-acid similarity table | J Mol Evol 36:96-99 | yes (seqinr's variance formulas) |

[Tutorial 5](tutorials/05-choosing-a-method.md) compares YN, NG, LWL, MYN and Comeron on 3,142 real ortholog
pairs, by divergence.

**`NG` is KaKs_Calculator's NG, not the textbook Nei-Gojobori count.** KaKs_Calculator takes a codon's synonymous
sites from its first and third positions only, leaves the changes that lead to a stop codon out of the
nonsynonymous sites, and then rescales S + N to the compared length; COAL does the same, so that `-m NG` gives
KaKs_Calculator's numbers. The classical NG86 procedure (as in MEGA, or PAML's NG86 output) counts the sites of
all three positions without that rescaling, and its dS differs: by up to about 0.02 on the 999 orthologr pairs
(measured by the independent review of 2026-11-12). Use `YN` (the default) when the absolute values matter, and
do not mix NG values of COAL / KaKs_Calculator with NG86 values of other programs.

The other KaKs_Calculator methods are refused with exit status 2 and the message "COAL does not implement
natively". These are `MS` and `MA` (model selection and averaging), `NONE`, the gamma
variants `GNG GLWL GLPB GMLWL GMLPB GYN GMYN`, and `ALL`. KaKs_Calculator's default method is MA. COAL's default is
YN, so `coal kaks` without `-m` gives YN rows.

`Comeron` and `Li` are not KaKs_Calculator methods. Their names are orthologr's, they are read in any case
(`-m comeron`), their rows come after those of KaKs_Calculator's methods, and they are implemented for the
standard genetic code only, as their references are (`-c` other than 1 with them is refused with exit status 2).
`Li` is not `LPB`: the two implement the same formulas of Li (1993) with different site and pathway counting, and give
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
| Multiple-hit correction undefined (log of a number that is not above 1e-10, so that an argument which is 0 in exact arithmetic counts as 0 whatever the rounding: "Differences", 14; for YN00/MYN, F84, K80 and JC69 all fail, or P + Q > 1) | that estimate is NA; the note says which (Ka, Ks). A YN or MYN row in which neither Ka nor Ks has an estimate has NA also in Sd and Nd when the pathways were weighted ("Differences", 17) | `saturated` |
| GY: the estimate of dS (or dN) is above 10 substitutions per site; this includes every pair whose branch length ends at its upper bound, 50 substitutions per codon ("GY: maximum likelihood") | that estimate and Ka/Ks are NA; the other one is kept. The note gives the estimate that was left out | `saturated` |
| GY: the optimiser does not converge (not seen in validation; a pair that does not settle because it is saturated is `saturated`) | no estimate: Ka, Ks and Ka/Ks are NA; the note has the last point | `not_converged` |
| MYN: kappa_TC or kappa_AG estimated as infinite (the TN93 transition term of a class of sites exactly at its limit) | nothing is estimated: Ka, Ks and the site counts are NA, and so is `coal kaks`'s rate ratio (Substitution-Rate-Ratio); the note names the kappa ("Differences", 15) | `saturated` |
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
  `1:1:1:1:1:1`. It is `NA` when a kappa is infinite (MYN, "Differences", 15): COAL's output has no `nan` or `inf`.
- `GC` is computed over both rows as given, with gaps included in the length.
- ML-Score, AICc and Model are filled for `GY` rows: the maximised log-likelihood, KaKs_Calculator's AICc
  (-2 lnL + 2 K n / (n - K - 1) with K = 3 parameters and n codons; `NA` when n is 4 or less) and `HKY`, the name
  KaKs_Calculator gives the model of its GY. For a GY row, Substitution-Rate-Ratio holds the estimated kappa,
  Divergence-Time is t / 3 (as for every method: (S Ks + N Ka) / (S + N)), and S-Substitutions / N-Substitutions
  are S Ks and N Ka scaled to the observed number of nucleotide differences, as KaKs_Calculator's GY prints them.
  The Method column says `GY` (KaKs_Calculator writes `GY-HKY`). For the other methods these columns are `NA`, and
  Akaike-Weight (model averaging, MA) always is.

For `Comeron` and `Li` the same 22 columns are written with what the method defines: Ka, Ks, Ka/Ks, Length,
Substitutions (nucleotide differences in the compared codons) and GC; for `Li` also `Fold-Sites(0:2:4)` (its
L0:L2:L4) and the two `Fold-*-Substitutions` columns (its transitional and transversional differences per
class). The other columns are `NA` (there are no S / N counts, so no Fisher test), and the rate-ratio column is
`1:1:1:1:1:1`, as for NG. The Method column says `Comeron` or `Li`.

Rows come per pair in KaKs_Calculator's method order (NG, LWL, MLWL, LPB, MLPB, GY, YN, MYN), then Comeron, then
Li. Pairs with status
`error` are left out, as KaKs_Calculator leaves out pairs it cannot read, and the summary on stderr counts them.
This format has no status column. A YN or MYN row whose iteration did not converge has `NA` in Ka, Ks, Ka/Ks,
P-Value(Fisher), S-Substitutions, N-Substitutions and Divergence-Time (KaKs_Calculator prints the values of its
last round there: difference 12); the summary on stderr says how many such rows there are, and `coal dnds` names
them.
`--digits N` sets the number of significant digits.

**What `P-Value(Fisher)` tests, and what has a standard error.** The p-value is Fisher's exact test on the
2 x 2 table of substitutions against sites that did not change, synonymous in one row and nonsynonymous in the
other: `[[Sd, S - Sd], [Nd, N - Nd]]`. It asks whether the proportion of changed sites differs between
synonymous and nonsynonymous sites, that is whether dN/dS differs from 1 for this pair; a small p with dN/dS
below 1 is evidence of purifying selection on the pair, not of anything about single sites, and with thousands
of pairs the p-values need a correction for multiple testing. The counts are estimates and not whole numbers
(KaKs_Calculator's convention, kept; "Differences" below says how the tails are summed), so read it as
KaKs_Calculator's statistic, comparable with its output, more than as an exact test. It is in the `coal kaks`
format only; `coal dnds` does not print it. **Standard errors** are given by the methods that have a variance
formula: `NG`, `YN`, `MYN` and `Li`, and by `GY` from the curvature of its likelihood (`SE_dN`, `SE_dS` in
`coal dnds`; `NA` for LWL, MLWL, LPB, MLPB and Comeron).
They are the standard errors **of dN and of dS, each on its own**. No method here gives a standard error or a
confidence interval for dN/dS itself: the ratio of two estimates has no simple variance, and dividing the two
SEs is not one. For an interval on the ratio use a likelihood method (PAML's codeml, HyPhy) on the codon
alignment, or bootstrap the codon columns.

**`coal dnds`** (`--digits` default 8):

| Column | Meaning |
|---|---|
| `seq1`, `seq2` | the two IDs. For FASTA input, the record IDs. For AXT, the name line split at the first `&`; when there is no `&`, seq1 is the whole name line and seq2 is `NA` |
| `method` | NG, LWL, MLWL, LPB, MLPB, GY, YN, MYN, Comeron or Li |
| `dN`, `dS`, `dNdS` | Ka, Ks, Ka/Ks. `NA` when not available; dNdS is `NA` when dS < 1e-6. Zeros print as 0 |
| `SE_dN`, `SE_dS` | standard errors of dN and of dS (YN, MYN, NG, Li, GY); `NA` elsewhere. There is none for the ratio (below) |
| `N_sites`, `S_sites` | nonsynonymous and synonymous sites (`NA` for Comeron and Li, which count sites by degeneracy class, not as S and N) |
| `Nd`, `Sd` | nonsynonymous and synonymous differences (`NA` for Comeron and Li) |
| `codons_compared` | codon columns used after the policy |
| `status`, `note` | see Policy: `ok`, `warning`, `not_converged`, `saturated`, `no_data`, `error` |
| `kappa`, `t`, `lnL` | only with `--ml-columns`, after `note`: the transition / transversion rate ratio (GY: the estimate; YN: the kappa it used; MYN: kappa from TN93), and GY's branch length t (expected substitutions per codon) and maximised log-likelihood; `NA` where a method has none |

**Input.** The format is detected from the first non-blank character (`>` means FASTA, anything else AXT), or set
with `--format`.

- **AXT** is KaKs_Calculator's format: a name line, the two aligned sequences, then a blank line. When a record
  has more than two sequence lines (wrapped sequences), they are joined and split in half, as KaKs_Calculator
  does, provided the half falls at the end of a line and every line holds only letters and the gap / missing
  symbols `- . ? * ~ !`. A record with **one** sequence line (the second sequence is missing), or with more lines
  that break either rule (how a blank line lost between two records shows: the next record's name line and
  sequences are joined to this one), is an error for that record: `coal dnds` gives it status `error`, `coal kaks`
  leaves it out with a warning, and the other records are computed. (KaKs_Calculator would cut the line, or the
  joined records, in half and compare the halves.) Lines starting with `#` before a record are skipped.
- **Aligned FASTA**: `--pairs all` (the default) compares every pair of the alignment in AXTConvertor's order
  (1-2, 1-3, ... 2-3, ...). `--pairs consecutive` compares records 1 and 2, 3 and 4, and so on. In `coal kaks`, a
  FASTA pair is named `id1&id2`, as KaKs_Calculator's AXTConvertor names it.
- Standard input (`-`) is not read, nor a pipe, FIFO or process substitution (`<(...)`); give a file.
- **The output must not be the input.** `-o` naming the input file (by any spelling: `./x`, a symbolic or hard
  link) is refused with exit status 2 before anything is written; the output is opened while the input is still
  being read, so the input would be destroyed. The same holds for `coal translate` and `coal codon-dnds`.

**Exit status.** 0 when the run finished and at least one result is an estimate or a statement about the pair
(`ok`, `warning`, `saturated`, `not_converged`), even if other pairs have a status other than `ok`. 1 when the
input cannot be read or the output cannot be written; the message names the file and the cause (`write error on
'out.tsv': No space left on device`). 2 for a command-line error, including a method that is not implemented,
**and when no pair has an estimate because of the input**: every result is `error` or `no_data`, or the input holds
no pair at all (an empty file, or blank or `#` comment lines only: `coal dnds: error: the input holds no pair of
sequences (...)`, also with `-q`). The table is
still written with every row (`coal dnds`; `coal kaks` leaves `error` pairs out, as KaKs_Calculator does), and one
line on stderr, printed also with `-q`, names the usual cause: the input looks like protein sequences, or the two
sequences of the pairs differ in length (the input is not aligned), with the way to a codon alignment. When only
some pairs differ in length the exit status stays 0 and a warning says how many. (KaKs_Calculator returns 0
in these cases too, so this is a difference for a script that reads its exit status; the numbers and the file
format are KaKs_Calculator's as before.)

**Compressed input.** `coal kaks`, `coal dnds` and `coal translate` read a gzip-compressed input (recognised by
its first two bytes, whatever the file is called): the `coal` command decompresses it into a temporary folder,
runs the engine on the copy and removes the folder. `coal-cpp` itself has no decompressor (it uses the C++
standard library only): called directly on such a file it stops with "is gzip-compressed: decompress it first
(gunzip -k FILE)". `-j N` sets the number of threads (default: all hardware threads). A one-letter option takes its
value attached as well (`-j4`, `-j=4`, `-mYN`, `-ofile`; flags chain: `-qj4`), and a long one after `=`
(`--jobs=4`), in `coal kaks`, `coal dnds`, `coal translate`, `coal codon-dnds` and `coal-cpp-align`, as in the
main `coal` command. Pairs are computed in parallel
and written in input order, and the output does not depend on `-j`.

`--comeron-flavour gestimator` and `--li-flavour seqinr` make `Comeron` and `Li` read the pair as their
reference programs do ("Comeron and Li" below); the default for both is `coal`, the policy above.
`--comeron-twosubs comeron|libsequence` chooses Comeron's pathways for codons that differ at two positions, and
`--comeron-min-sites N` the smallest site class whose undefined correction makes dN or dS `NA` (a smaller class
cannot, unless the distance would then be 0; "Two-position codon differences" and "Small site classes" below).

`--yn-flavour paml` and `--yn-weighting 0|1` choose which reference YN00 follows where the two differ (below).
`--yn-weighting 0` gives equal pathway weights in one round, which is the default control file of PAML's yn00
(`weighting = 0`). KaKs_Calculator always weights pathways.

## What to cite

The numbers of `coal kaks`, `coal dnds` and `coal pipeline` come from published estimators, implemented in code
derived from KaKs_Calculator 2.0, PAML's yn00, libsequence and seqinr; the pipeline also runs DIAMOND and an
aligner. A paper should cite the methods and programs its numbers come from, with COAL (and the version
`coal --version` prints):

| What you ran | Cite |
|---|---|
| Anything | COAL |
| A codon alignment (`coal PEP NUC`, `pal2nal`, every Ka/Ks route of the pipeline) | Suyama et al. 2006 (PAL2NAL) |
| Ka/Ks, method NG | Nei & Gojobori 1986; Wang et al. 2010 (KaKs_Calculator 2.0) |
| Ka/Ks, method LWL | Li, Wu & Luo 1985; Wang et al. 2010 (KaKs_Calculator 2.0) |
| Ka/Ks, method MLWL | Li, Wu & Luo 1985; Tzeng, Pan & Li 2004; Wang et al. 2010 (KaKs_Calculator 2.0) |
| Ka/Ks, method LPB | Li 1993; Pamilo & Bianchi 1993; Wang et al. 2010 (KaKs_Calculator 2.0) |
| Ka/Ks, method MLPB | Li 1993; Pamilo & Bianchi 1993; Tzeng, Pan & Li 2004; Wang et al. 2010 (KaKs_Calculator 2.0) |
| Ka/Ks, method YN (the default) | Yang & Nielsen 2000; Wang et al. 2010 (KaKs_Calculator 2.0); Yang 2007 (PAML 4) |
| Ka/Ks, method MYN | Zhang, Li & Yu 2006; Yang & Nielsen 2000; Wang et al. 2010 (KaKs_Calculator 2.0); Yang 2007 (PAML 4) |
| Ka/Ks, method GY (maximum likelihood) | Goldman & Yang 1994; Yang & Nielsen 2000; Yang 2007 (PAML 4) |
| Ka/Ks, method Comeron | Comeron 1995; Thornton 2003 (libsequence) |
| Ka/Ks, method Li | Li 1993; Charif & Lobry 2007 (seqinr) |
| `coal pipeline pairs` with DIAMOND (orthologs) | Buchfink et al. 2021 (DIAMOND) |
| `coal pipeline clusters` with DIAMOND DeepClust | Buchfink et al. 2021 (DIAMOND); Buchfink et al. 2026 (DIAMOND DeepClust) |
| Protein alignments by Clustal Omega (`--aligner clustalo`, the cluster default) | Sievers et al. 2011 (Clustal Omega) |
| Protein alignments by MAFFT (`--aligner mafft`) | Katoh & Standley 2013 (MAFFT) |
| Through the R package orthologr | Drost et al. 2015 (orthologr) |

`coal cite` prints this table and the references; `coal cite YN` those of one method; `coal cite out/report.json`
those of one pipeline run, whose report lists the keys as `cite`.

- Suyama M, Torrents D, Bork P (2006). PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Research 34:W609-W612. <https://doi.org/10.1093/nar/gkl315>
- Nei M, Gojobori T (1986). Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. Molecular Biology and Evolution 3:418-426. <https://doi.org/10.1093/oxfordjournals.molbev.a040410>
- Li W-H, Wu C-I, Luo C-C (1985). A new method for estimating synonymous and nonsynonymous rates of nucleotide substitution considering the relative likelihood of nucleotide and codon changes. Molecular Biology and Evolution 2:150-174. <https://doi.org/10.1093/oxfordjournals.molbev.a040343>
- Li W-H (1993). Unbiased estimation of the rates of synonymous and nonsynonymous substitution. Journal of Molecular Evolution 36:96-99. <https://doi.org/10.1007/BF02407308>
- Pamilo P, Bianchi NO (1993). Evolution of the Zfx and Zfy genes: rates and interdependence between the genes. Molecular Biology and Evolution 10:271-281. <https://doi.org/10.1093/oxfordjournals.molbev.a040003>
- Tzeng Y-H, Pan R, Li W-H (2004). Comparison of three methods for estimating rates of synonymous and nonsynonymous nucleotide substitutions. Molecular Biology and Evolution 21:2290-2298. <https://doi.org/10.1093/molbev/msh242>
- Goldman N, Yang Z (1994). A codon-based model of nucleotide substitution for protein-coding DNA sequences. Molecular Biology and Evolution 11:725-736. <https://doi.org/10.1093/oxfordjournals.molbev.a040153>
- Yang Z, Nielsen R (2000). Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Molecular Biology and Evolution 17:32-43. <https://doi.org/10.1093/oxfordjournals.molbev.a026236>
- Zhang Z, Li J, Yu J (2006). Computing Ka and Ks with a consideration of unequal transitional substitutions. BMC Evolutionary Biology 6:44. <https://doi.org/10.1186/1471-2148-6-44>
- Comeron JM (1995). A method for estimating the numbers of synonymous and nonsynonymous substitutions per site. Journal of Molecular Evolution 41:1152-1159. <https://doi.org/10.1007/BF00173196>
- Wang D, Zhang Y, Zhang Z, Zhu J, Yu J (2010). KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies. Genomics, Proteomics & Bioinformatics 8:77-80. <https://doi.org/10.1016/S1672-0229(10)60008-3>
- Yang Z (2007). PAML 4: phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution 24:1586-1591. <https://doi.org/10.1093/molbev/msm088>
- Thornton K (2003). libsequence: a C++ class library for evolutionary genetic analysis. Bioinformatics 19:2325-2327. <https://doi.org/10.1093/bioinformatics/btg316>
- Charif D, Lobry JR (2007). SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis. In: Bastolla U, Porto M, Roman HE, Vendruscolo M (eds), Structural Approaches to Sequence Evolution: Molecules, Networks, Populations. Springer, pp. 207-232. <https://doi.org/10.1007/978-3-540-35306-5_10>
- Buchfink B, Reuter K, Drost H-G (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods 18:366-368. <https://doi.org/10.1038/s41592-021-01101-x>
- Buchfink B, Barbé É, Ashkenazy H, Reuter K, Kennedy JA, Drost H-G (2026). Clustering the protein universe of life using DIAMOND DeepClust. Nature Methods 23:724-727. <https://doi.org/10.1038/s41592-026-03030-z>
- Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Soding J, Thompson JD, Higgins DG (2011). Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Molecular Systems Biology 7:539. <https://doi.org/10.1038/msb.2011.75>
- Katoh K, Standley DM (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular Biology and Evolution 30:772-780. <https://doi.org/10.1093/molbev/mst010>
- Drost H-G, Gabel A, Grosse I, Quint M (2015). Evidence for active maintenance of phylotranscriptomic hourglass patterns in animal and plant embryogenesis. Molecular Biology and Evolution 32:1221-1231. (the R package orthologr) <https://doi.org/10.1093/molbev/msv012>

Each DOI is that of the cited work, checked against the publisher's or the institution's page of the paper. COAL
itself has no DOI: cite its version and repository (`CITATION.cff`). The DeepClust preprint (bioRxiv 2023,
<https://doi.org/10.1101/2023.01.24.525373>, key `buchfink2023` in reports written by earlier versions) is
superseded by the paper above; `coal cite` still resolves the old key and prints the paper after it.

## Filtering before interpretation

A row with status `ok` says that the method produced numbers, not that the numbers mean much. Before dN/dS values
are read, compared or averaged, four things are worth checking for every pair. They are the usual reasons for
impossible-looking values in a pairwise table.

| What | Why it matters | Where to look |
|---|---|---|
| **Few compared codons** | With a few dozen codons a single substitution moves dN/dS a long way, and the standard errors are of the size of the estimates. In orthologr's *Arabidopsis* example one pair is compared over 46 codons and has dN/dS 1.84 with status `ok`: that is noise, not positive selection | `codons_compared` |
| **Low coverage or identity** | When only a small part of the two proteins aligns, or they align at 30% identity, the pair may not be orthologous over its length and the alignment itself is uncertain; misaligned codons count as substitutions and inflate dN and dS | `codons_compared` against `prot_aln_length`; `prot_identity` (pipeline and coalr tables) |
| **Saturated dS** | From a dS of about 1 (about 0.5 with strong codon-usage bias) many synonymous sites have changed more than once and the correction for multiple hits dominates the estimate. In simulations with a known answer ("Accuracy against simulated truth" below) YN and MYN then overestimate dS (by 10 to 30% with even codon usage, 25 to 60% with skewed usage) and underestimate dN/dS (by about 20%, up to 50%), even at 1,000 codons; above a dS of about 2.5 the pairs that still have an estimate have dN/dS too high (+30% median at dS 3). `saturated` marks only the pairs where the correction formula has no value at all: a pair with dS 3 or 5 and status `ok` is not better (268 such pairs with dS above 2 stood next to 3,748 `saturated` ones in a yeast *S. cerevisiae* / *C. glabrata* run) | `dS` |
| **Very small dS** | dN/dS is a ratio: between very close sequences dS rests on one or two synonymous differences, and dividing by it gives values that are large, variable and often above 1 for no biological reason. In the simulations, at a true dS of 0.01 only 3 to 57% of the YN ratios of 300-codon genes are within ±50% of the truth. `dNdS` is `NA` only when dS is below 1e-6 | `dS`, `Sd`, `SE_dS` |

**Suggested starting thresholds.** These are conventions in common use, not properties of your data: state the
ones you used, and check that your conclusion does not depend on their exact values.

- at least 50 compared codons, and 100 for anything read gene by gene, dN/dS included;
- compared codons covering at least half of the protein alignment, and protein identity of at least 30 to 40%;
- dS at most 1 (`--max-ds 1`), or 0.5 where codon usage is strongly biased (GC3 far from 50%); 2 only as an outer
  bound beyond which nothing is usable; and for the ratio, dS at least about 0.05;
- status `ok` or `warning` only (`saturated`, `not_converged`, `no_data` and `error` have no estimate), and for
  `warning` a look at the note (an internal stop codon or a CDS that is not a multiple of 3 often means a
  pseudogene or an annotation error).

A whole-gene pairwise dN/dS averages over all sites and the whole time since divergence. A value above 1 is rare
and, from one pair, weak evidence; site and branch models (PAML codeml, HyPhy) test for positive selection, and
COAL's codon alignments are their input.

**The flags.** `--min-codons N` and `--max-ds X` (`coal dnds`, `coal pipeline`; `min_codons`, `max_ds` in coalr)
mark the rows concerned without removing anything. Where a protein alignment exists (`coal pipeline`,
`coal codon-dnds`; `min_identity` in coalr's `coal_dnds_pairs()` and `coal_dnds_clusters()`) there is a third,
`--min-identity P`: a pair whose protein identity (`prot_identity`, percent) is below P gets the note
`protein identity below P%`. `coal dnds` reads codon alignments and refuses it.

```
coal dnds codon.fasta --min-codons 100 --max-ds 1 -o dnds.tsv
coal pipeline pairs --query q_cds.fa --subject s_cds.fa -o out --min-codons 100 --max-ds 1
```

A row that has an estimate (status `ok` or `warning`) and fewer than N compared codons gets the note
`fewer than N codons compared`; one whose dS is above X gets `dS above X`; the notes are appended to what the
row already says (`; ` between them), and its status becomes `warning` if it was `ok`. **No row is dropped and no
number changes**: the table is the one without the options but for `status` and `note` of the flagged rows, and
the per-cluster summary, which counts `ok` and `warning` together, is the same. Rows without an estimate are left
as they are. Both options are off by default. stderr and, for the pipeline, `report.json`
(`counts.flagged`) say how many rows each flag marked. Filtering is then one line in whatever reads the table:

```r
d <- read.delim("out/dnds.tsv")
keep <- d$status == "ok" & d$dS >= 0.05 & d$codons_compared / d$prot_aln_length >= 0.5 & d$prot_identity >= 40
```

`coal kaks` has no such options: KaKs_Calculator's format has no status column.

Tutorial 2 applies these filters step by step
([docs/tutorials/02-pairs-dnds.md#5-filter-then-interpret](tutorials/02-pairs-dnds.md#5-filter-then-interpret)).

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
if (!coal::kaks::parse_method("MA", m, why)) { /* why: "... does not implement natively ..." */ }
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
`Method::GY` fills `lnl`, `t` and `omega` as well (NaN for the other methods), puts its estimate of kappa in
`kappa`, and counts the optimiser's steps in `iterations`.

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
methods; from 2026-11-14 to 15, 7 LPB / MLPB site counts of `saturated` rows, difference 14, which are `NA`
since difference 19). **That is not true for short pairs in general:** on the 3,000 pairs with codons drawn independently,
YN and MYN print different numbers in 23 fields (Ka, Ks or Ka/Ks), which belong to differences 11 and 13 below,
MLWL in 36 fields of 10 pairs whose kappa is exactly 2 or 0.5 (difference 18, since 2026-11-15; NG, LWL, LPB and
MLPB in none), and in 210 rows (123 YN, 87 MYN) COAL prints no estimate at all because
the iteration did not converge (difference 12), where KaKs_Calculator prints its last round.

**COAL `--yn-flavour paml` compared with PAML yn00, 12 significant digits.** The result is the same for
weighting 1 and 0.

| Set | dN max abs diff | dS max abs diff | S, N sites max abs diff | SE max abs diff | Pairs > 1e-4 |
|---|---:|---:|---:|---:|---:|
| orthologr 999 | 4.8e-13 | 5.0e-13 | 5.0e-9 | 4.7e-13 | 0 |
| simulated 600 | 5.0e-13 | 5.0e-12 | 5.0e-9 | 4.5e-11 | 0 |

These differences are the 12-digit printing of yn00.

**The K80 boundary 2P + Q = 1.** At the non-degenerate or 4-fold sites of a pair the proportions of transitions
(P) and transversions (Q) can be fractions whose sum 2P + Q is exactly 1 (P = 11/34, Q = 12/34). Whether
`1 - 2P - Q` then comes out as 0 (K80 undefined: kappa from the other kind of sites) or as a residue of 1e-17
(kappa about 60) depends on the order in which the floating-point proportions are formed. `--yn-flavour paml`
forms them as yn00's `GetKappa` does (Q = 1 − identical − P) and so agrees with yn00 at the boundary; the
default flavour keeps KaKs_Calculator's order (1 − P − identical) and its numbers. Such pairs are not rare in
distant species (319 of the 4,388 *S. cerevisiae* / *C. glabrata* pairs are at the boundary; for 10 of them the
two orders give different kappa), and their kappa from the 4-fold sites is a rounding artefact in either
program. Deep test of 2026-11-14, finding D01; the pair is `tests/kaks/yn_paml_boundary_pair.fasta`. The
tolerance of difference 14 below is not applied to the kappa estimate (of YN or MYN), so that both flavours keep
their reference's kappa.

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
program is right, or that neither number is an estimate. Differences 14 to 19 are boundaries, an infinite kappa,
genetic codes other than the standard one and the counts of a saturated LWL-family row (the section before 14).

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
   S = Sd/Ks and N = Nd/Ka for the output, and prints `-nan` when the divisor is 0 (identical sequences, Ka = 0,
   or a Ks of 0 because the classes it is made from are saturated and their K set to 0); COAL prints `NA` there.
   Where a saturated class leaves Ks (Ka) a number, KaKs_Calculator prints the quotient, and COAL `NA`
   (difference 19). It also covers pairs with nothing to compare.
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
| 14, LPB / MLPB S-Sites or N-Sites of a `saturated` row beside a class at an exact boundary (2026-11-14 to 15: 21 and 21; since then `NA`, difference 19) | 0 | 0 | COAL (the reference's quotient comes from a rounding residue) |
| 18, MLWL with kappa exactly 2 or 0.5: the other site formula (10 pairs; since 2026-11-15) | 36 | 36 | COAL (the formula of the exact kappa; the reference's choice follows the rounding of its column-order sums) |
| unexplained | 0 | 0 | |

So "the corrected reference gives COAL's numbers" holds for the 4,000 pairs with one changed position per codon
and for the second row of the 3,000; it says nothing where the iteration does not converge (COAL has no number
there; the corrected reference prints more last rounds than the unpatched one, because it prints numbers where
the unpatched one prints `nan`), and it does not hold for pairs without a site of one kind, which the corrected
lines do not touch. NG, LWL, LPB and MLPB have no difference between two numbers on either random set; MLWL has
the 36 fields of difference 18 on the second (10 pairs with kappa exactly 2 or 0.5) and none on the first.

**KaKs_Calculator's MYN depends on the order of the pairs in a file.** It carries state from one pair to the
next: run once over the 3,000 pairs, 49 MYN rows differ (Ka, Ks or Ka/Ks) from the same pairs run one per process
(the other six methods: none; the corrected reference: none). The comparisons of this section therefore run
KaKs_Calculator once per pair (`random_short_pairs.py` does; `--one-file` reports the count). The tables for the
real genes above were made with one run per file; there the stock and the 15-digit outputs agree with COAL to
the stated tolerances, so no such row occurs in them. COAL has no state between pairs: every pair gives the same
row alone, in any order and with any number of threads (`tests/test_kaks.py`, `test_threads_same_output`).

#### Boundaries, an infinite kappa and other genetic codes (2026-11-14)

Differences 14 to 17 (and 18, 2026-11-15; 19 is from the review of that day) came out of the property-based fuzz of the engine (`tests/test_kaks_properties.py`, findings
K01 to K06; the minimal pairs are in `tests/fixtures/kaks_fuzz/`) and of the simulation study above. With the
standard code they change only rows at an exact boundary of a correction, MYN rows with an infinite kappa and
the difference counts of YN / MYN rows without any estimate; every other row is byte-identical
(checked at 17 digits, all nine methods, on the *A. thaliana* / *A. lyrata*, human / mouse, *E. coli* K-12 /
Sakai, *Salmonella* / *Shigella* and *S. cerevisiae* / *C. glabrata* pairs of the deep test, the 600 simulated
pairs of "Validation" (no row changed), 4,000 pairs of the simulation study and the test fixtures).

14. **A correction whose argument is 0 in exact arithmetic is `saturated`.** Where the argument of a correction's
    logarithm is 0 in exact arithmetic (Kimura: 1 - 2P - Q = 0, for instance P = Q = 1/3; Jukes and Cantor: p =
    3/4), the floating-point value is 0, a residue of about 1e-16 or a small negative number, depending on the
    fractions and on the order of the sums. A residue gave a distance of 19 to 28 that was printed as an estimate:
    `GCTGCTGCT` / `GCCGCAGCT` had LWL dS 18.99 with status `ok` (the same pair with four codons, whose fractions
    are exact in binary, was `saturated`), and `TCA` / `TTT` NG dS 27.55. COAL now counts an argument at or below
    1e-10 as 0 (a distance from a larger one is below about 12 for K80 and 17 for JC69), in NG, the LWL family, the
    F84 / K80 / JC69 and TN93 distances of YN and MYN, and the default flavours of Comeron and Li; the result no
    longer depends on the order of the codon columns. KaKs_Calculator prints the residue's number (for YN / MYN it
    loses such pairs to NaN instead, difference 11), so these are rows where it prints a number and COAL `NA`
    (difference 1). For YN and MYN the fallback chain then works as at an exact 0: `TAT` / `TGG` (no synonymous
    site) has F84 and K80 undefined and JC69 defined, Ka = 1.648 (it was 18.99). Not changed: the kappa estimates
    (D01 above) and the reference flavours `--comeron-flavour gestimator` and `--li-flavour seqinr`, which keep
    their programs' numbers, residue included. In the same way, Comeron's "a class with fewer than
    `--comeron-min-sites` sites" is tested with a margin of 1e-9, because a class of exactly 3 sites can sum to
    2.9999999999999996 in one column order and 3 in another. What changed: of the 1,830 one-codon pairs of sense
    codons, 355 rows had a dN or dS above 10 (36 of them `ok`); one remains (MLWL, whose dS divides by 0.12 kappa-
    based sites: a real value of the method). In general MLWL's dS can be several times LWL's for the same pair
    without any residue, because its synonymous sites are kappa-based: with no transversion kappa is 1 and S is a
    fifth of the 2-fold sites (plus the 4-fold sites), so 513 Gln codon pairs CAG / CAA with 256 differences and
    no 4-fold site give dS = 5 x ln(513) / 2 = 15.60 (LWL 9.36), and as kappa approaches 0.5 S approaches the
    4-fold sites alone (400 CAG / CAA, 1,590 CAA / AAA and 10 CAA / CAA: kappa 0.5031, S 3.14, dS 162.70); both
    are what KaKs_Calculator prints (property fuzz seed 4243, FIXCI-B3, and the audit of 2026-11-15). Real genes: 22 rows of the 4,388 *S. cerevisiae* / *C. glabrata*
    pairs (LWL family; dS 12 to 23 with status `ok` before, or a Ka beside a class with the residue) and two
    Comeron rows of the *Salmonella* / *Shigella* pairs (a 3-site class: `NA`, `saturated`, as with 3 sites in any
    order); none in the other genomes. Against KaKs_Calculator on the random short pairs
    (`benchmarks/kaks/random_short_pairs.py`, 2026-11-14): 193 (one change per codon) and 231 (codons drawn
    independently) more fields where it prints a number and COAL `NA`, all with status `saturated`, and 7 and 21
    LPB / MLPB site counts of `saturated` rows that differ (S = Sd / Ks with the residue's class in Sd); no
    unexplained difference. (Since 2026-11-15 those LPB / MLPB site counts are `NA`, difference 19.)
15. **MYN: kappa_TC or kappa_AG infinite.** MYN estimates the two kappas with TN93 at the non-degenerate and the
    four-fold sites. When the A-G (T-C) term of one class is exactly at its limit, 1 - R P1 / (2 piA piG) -
    Q / (2 R) = 0 (P1 = 3/25, Q = 12/25 with A and G equally frequent), its logarithm is -inf and kappa_AG is
    infinite. The site counts were then NaN, and the row had status `no_data` with the note `no synonymous sites;
    no nonsynonymous sites` although the pair has its full number of sites; `coal kaks` on such pairs alone exited
    2 and suggested that the input might not be a codon alignment. The status is now `saturated`, with the note
    `MYN: kappa_AG estimated as infinite (TN93 transition correction undefined in a class of sites): no site
    counts, Ka and Ks NA`, and the Substitutions column keeps the number of differences. KaKs_Calculator prints
    `nan` for Ka, Ks, Ka/Ks and the site counts, and the rate ratio with `inf` for that kappa
    (`kappa_TC:inf:1:1:1:1`). COAL prints `NA` for all of them, the rate ratio included, as it prints an undefined
    Fold-Sites or GC field: the whole field `NA`, because its output has no `nan` or `inf` (since 2026-11-15; before,
    the rate ratio showed the `inf` as KaKs_Calculator does. The property fuzz's P3 found it, seed 8585, on a pair
    whose kappa_TC is infinite: `inf:0.994949:1:1:1:1`, now `tests/kaks/myn_kappa_tc_inf.axt`; among the deep
    test's genomes it is one *Salmonella* / *Shigella* row, p2067). The C++ API keeps the infinite value in
    `Result::rates`. In the simulation
    study 79 of 240,000 MYN rows (dS 0.5 to 3, all lengths); among real genes one *Salmonella* / *Shigella* and
    two *S. cerevisiae* / *C. glabrata* pairs. `tests/test_kaks_myn_kappa_inf.py` (and, through the C++ API,
    `tests/cpp/kaks_myn_kappa_inf.cpp`). A residue instead of the exact 0 gives a very large but finite kappa and
    an estimate, as in KaKs_Calculator: like D01, the kappa keeps the reference's exact test.
16. **Genetic codes other than the standard one: the special cases of the LWL family.** KaKs_Calculator
    hard-codes rules for codons of the standard code: LWL and LPB split a CGA <-> AGA (CGG <-> AGG) change at the
    first position into half a transition and half a transversion, and MLWL and MLPB (Tzeng et al. 2004) count
    ATA <-> ATG as nonsynonymous and ATA <-> ATC, ATT and CGA <-> AGA as synonymous; it applies them with every
    code. In the vertebrate mitochondrial code (2), where ATA is Met, MLWL and MLPB therefore counted half of a
    synonymous ATA -> ATG as nonsynonymous (ATAGCTAAACTGGGT / ATGGCCAAACTGGGT: dN 0.154 and 0.137 instead of 0),
    and in code 13, where AGA is Gly, LWL and LPB depended on which sequence was first (dN 0.109 against 0.051).
    COAL applies the LWL / LPB rule only where the code makes the two codons synonymous (an ordinary
    nonsynonymous transversion otherwise, and then the order does not matter), and the M variants' rules by what
    they are about: between two codons that are both 2-fold degenerate at the position, the change is synonymous
    exactly when the code says so. In the standard code those are exactly Tzeng et al.'s Ile and Arg cases, so its
    numbers are KaKs_Calculator's, unchanged; other codes have their own such cases (Asn AAA <-> AAC, AAT in codes
    9, 14 and 21; Cys TGA <-> TGC, TGT in code 10; Leu, Tyr and Ser codons in codes 12, 14, 16 and 22), which
    COAL now counts as synonymous too. ATA <-> ATG is nonsynonymous where the code has Ile and Met there. With
    those codes COAL's MLWL, MLPB, LWL and LPB differ from KaKs_Calculator's.
17. **No weighted difference counts without an estimate.** A YN or MYN row in which neither Ka nor Ks has an
    estimate (both saturated, or one saturated beside no sites of the other kind) used to print Sd and Nd, the
    pathway-weighted counts of an iteration whose t was a cap (99, or JC69's cap: difference 13). Like the last
    round of an iteration that did not converge, they depended on the path, and swapping the two sequences changed
    them (Nd 23.09 against 23.18 for a pair of code 5). They are now `NA` in such rows, as for `not_converged`
    (and so is Fisher's p, which needs them). Without codons that differ at several positions, and with
    `--yn-weighting 0`, the counts do not depend on t and are printed as before.

18. **MLWL at kappa exactly 2 or 0.5.** MLWL (Tzeng et al. 2004) picks its formula for the synonymous and
    nonsynonymous sites by comparing kappa = 2 ts / tv with 2 and with 0.5: S = (kappa - 1) L2 / (kappa + 1) + L4
    above 2, (kappa - 0.5) L2 / (kappa + 1.5) + L4 above 0.5, L2 / 3 + L4 otherwise (N likewise), and the formulas
    do not meet at the thresholds (at kappa = 2 the first gives L2 / 3, the second 3 L2 / 7). ts and tv are sums
    of halves and thirds of a difference per codon column, so kappa is often exactly 2 or 0.5, and exact kappa = 2
    is "not above 2": the middle formula. Summed in floating point in the order of the columns, kappa came out as 2
    or as 2.0000000000000004 (0.5 or 0.50000000000000011) according to that order, and KaKs_Calculator's choice of
    formula, like COAL's before, depended on it: permuting the codon columns of a pair could change N and S by a
    site or more, and dN and dS with them (property fuzz P6, seed 4242, FIXCI-B2: `SATACGTGCGTG---TCCAACCTA` /
    `AATACGTGACCTGGGAGACACATG`, code 13, N / S 13.67 / 4.33 or 15 / 3 by the order). COAL now sums ts and tv
    exactly as well (in integer units of 1/120 of a difference: a column adds j / 2 / d, d = 1 to 6 pathways) and
    makes the two comparisons on those sums, so every order of the columns takes the formula of the exact kappa;
    the numbers printed (kappa, the Fold-* counts) are the same doubles as before. Where the rounding goes the same
    way in every order, KaKs_Calculator's number is then not COAL's: `mlwl_kappa2_code1` of
    `tests/kaks/mlwl_kappa_boundary.axt` (code 1, ts = tv = 3.5, kappa = 2) has N / S 13.67 / 4.33 in
    KaKs_Calculator in every order of its columns, and 13.29 / 4.71 in COAL. In real genes: *S. cerevisiae* /
    *C. glabrata* pair p2754 of the deep test (code 1; L0 : L2 : L4 = 188 : 66.5 : 42.5, ts = tv = 49 exactly, tv
    summed as 48.99999999999999 in the file's order) had S 64.667, N 232.333, Ka 0.249535, Ks 1.81172, as
    KaKs_Calculator prints for that order; with the middle formula S = 3 x 66.5 / 7 + 42.5 = 71 and N = 188 +
    2 x 66.5 / 3.5 = 226, so Ka 0.256528 and Ks 1.65011, which is also what KaKs_Calculator prints for 38 of 40
    random orders of the same codon columns. It is the only MLWL row that changes in the deep test's genomes (none
    in *A. thaliana* / *A. lyrata*, human / mouse, *E. coli* K-12 / Sakai, *Salmonella* / *Shigella*), and none of
    the 999 orthologr and 600 simulated pairs of the tables above changes. Among random short pairs: none of the
    4,000 with one changed position per codon, 10 of the 3,000 with codons drawn independently (36 fields:
    `random_short_pairs.py` classes them as "MLWL kappa exactly 2 or 0.5"), and 10 to 22 of 20,000 per genetic
    code in a separate run; every such row has kappa exactly 2 or 0.5. No other method changes.
    `tests/test_engine_review_fixes.py` (`TestMlwlKappaBoundary`: 40 orders of each pair of
    `tests/kaks/mlwl_kappa_boundary.axt`).
19. **LWL family: no substitution counts from a saturated correction** (review of 2026-11-15). In a saturated
    class KaKs_Calculator, and COAL with it, sets A, B and K to 0, and Sd = L2 A2 + L4 K4, Nd = L0 K0 + L2 B2
    were summed with those zeros: numbers that look like estimates of substitutions and are not, and LPB / MLPB's
    S = Sd / Ks and N = Nd / Ka quotients of them (`TAT...` code 11, 11 codons: LPB S 11.5 and N 31, 42.5 sites;
    NP_459032.1 of the *Salmonella* / *Shigella* pairs: LPB N 537.72 beside Nd 245.54). COAL now prints `NA` for
    Sd (`S-Substitutions`) when Ks is saturated and for Nd (`N-Substitutions`) when Ka is, for LPB / MLPB's
    S-Sites and N-Sites under the same conditions, and so for Fisher's p, which needs them; KaKs_Calculator
    prints the numbers. Kept: LWL / MLWL's S and N (counts of the fold classes, with MLWL's kappa from the raw
    transitions and transversions: not estimates, as for a row that is not saturated), the Fold-S- and
    Fold-N-Substitutions (the raw differences per class), Sd or Nd on the side that is not saturated, and every
    row that is not `saturated`. It is the LWL family's counterpart of difference 17. All 9 methods at 17 digits,
    dnds and kaks formats: *Salmonella* / *Shigella* 530 of 25,758 rows change (LWL 129, LPB 129, MLWL 136, MLPB
    136: every saturated row of the four methods, and nothing else), *S. cerevisiae* / *C. glabrata* 8,694 of
    39,492 (2,079, 2,079, 2,268, 2,268; likewise). Against KaKs_Calculator these are fields where it prints a
    number and COAL `NA`; `tests/fuzz_kaks_vs_reference.py` classes them as difference 19 (405 fields of 500
    pairs, seed 77115; 1,242 of 1,500 pairs with codes 1, 4 and 11, seed 990017). `tests/test_kaks.py`
    (`TestSaturatedLwlCounts`), `tests/kaks/golden/edge_cases.coal.*`.

`tests/test_kaks.py` holds the cases (`TestOneKindOfSite`, `TestMynNan`, `TestMultiPositionCodons`), and its
comparison with KaKs_Calculator's stored output does not accept a difference because of the status alone: a
number that KaKs_Calculator prints and COAL does not is a failure unless the status is `saturated`
(difference 1) or `not_converged` (difference 12; only in the columns that come out of the iteration).

The categories, tolerances and checks are those of `tests/test_kaks.py`, `benchmarks/kaks/validate_kaks.py` and
`benchmarks/kaks/random_short_pairs.py`. On the real-gene sets (999 + 600 pairs) and on the 4,000 random pairs
with one changed position per codon, all fields agree to 1e-9 apart from differences 1 to 10 and 14. On short
pairs whose codons differ at several positions they do not: differences 11 to 14 and 18. On real genes MLWL can
differ at an exact kappa of 2 or 0.5 (difference 18: one *S. cerevisiae* / *C. glabrata* pair).

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
- Swapping the two sequences of a pair, or permuting its codon columns, gave the same statuses and the same
  numbers at the default digits in every pair of the fuzz campaigns (MLWL since 2026-11-15: difference 18); at `--digits 17` the last one or two digits of NG, YN, MYN and Comeron (default
  flavours) can change, because sums are formed in another order, and the YN / MYN iteration can carry that to
  about 3e-10 relative (70 of 1,593 rows of 300 random pairs at 17 digits, none at 8 or 6; fuzz finding K04).
- An iteration that did not converge prints no value (difference 12).
- On x86 with GCC or Clang, YN / MYN's transition-probability matrix is computed by an AVX2 copy of the code
  where the CPU has AVX2 (chosen at run time; the baseline copy elsewhere and on ARM and MSVC builds). It performs
  the same operations in the same order as the baseline copy (no FMA, no reassociation), so the numbers are the
  same bits: `tests/test_kaks_dispatch.py` compares the two copies byte for byte (the testing options
  `--kaks-baseline`, which forces the baseline copy, and `--kaks-isa`, which prints `avx2` or `baseline` and
  exits; not in the help, like `--yn-always-pmat`).

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

On these 600 pairs (one kappa and one codon usage, t up to 2 substitutions per codon) YN00 and MYN recover omega to within about
10 %; at omega = 0.1 and t = 1 they run about 20 % low, and NG and LWL are low. That does not generalise: the
larger study below finds YN00 and MYN unbiased only up to a dS of about 0.5, and the sign of the counting methods'
error set by the codon usage (too high with skewed usage). NG's 27 NA are Ks saturations (Jukes-Cantor) at
omega = 0.1, t = 2.

### Accuracy against simulated truth

The comparisons above say that COAL computes what KaKs_Calculator, yn00, gestimator and seqinr compute. How close
those estimators are to the truth was measured separately (2026-11-14, COAL 3232dd9; since then only the D9
refinement of Comeron's small classes and the rows of "Differences" 14 to 17 changed, in few of the 240,000
pairs, among them the 79 MYN rows of difference 15).

**Truth.** Pairs were simulated under GY94 / M0 with F3x4 codon frequencies (the rate matrix scaled to one
substitution per codon). The true dS and dN are those of the model in Yang & Nielsen's (2000) definition, which
YN00 and codeml estimate: the synonymous sites S are counted with omega set to 1 (S = 3 rhoS(omega = 1)),
dS = t rhoS / S and dN likewise, so that dN/dS is omega exactly. The counting methods (NG, the LWL family,
Comeron, Li) count sites without codon frequencies, so when kappa is not 1 or codon usage is skewed they estimate
a different quantity; their "error" here is their distance from the model's dS, dN and omega. The simulator was
checked against the model's algebra, against PAML's evolver (40 distribution tests) and with codeml, whose
pairwise maximum likelihood recovers the truth to within 8 % at every dS from 0.1 to 3.

**Grid.** dS 0.01, 0.03, 0.1, 0.2, 0.5, 1, 1.5, 2, 2.5 and 3; omega 0.05, 0.2, 0.5, 1 and 2; kappa 1, 2 and 4;
30, 100, 300 and 1,000 codons; even codon usage or a GC3-skewed one (third positions 80 % G/C); 200 pairs per
cell, 240,000 pairs.

**Bias.** Median relative error, over the 30 combinations of codon usage, kappa and omega (the median of the
cell medians; dS / dN/dS):

| Method | codons | dS 0.1 | dS 0.5 | dS 1 | dS 2 | dS 3 |
|---|---:|---|---|---|---|---|
| YN | 1,000 | -1% / +0% | +3% / -6% | +18% / -21% | +14% / -30% | -16% / -14% |
| MYN | 1,000 | -2% / +1% | +3% / -6% | +18% / -21% | +11% / -21% | -25% / -11% |
| NG | 1,000 | -9% / +11% | -12% / +5% | -17% / +4% | -24% / +5% | -55% / +22% |
| LPB | 1,000 | -11% / +23% | -6% / +18% | -8% / +16% | -11% / +20% | -30% / +32% |
| Comeron | 1,000 | -0% / -1% | +6% / -8% | +8% / -11% | +8% / -16% | -24% / +15% |
| Li | 1,000 | -7% / +8% | -0% / -1% | -0% / -1% | -2% / -3% | -20% / +8% |
| GY | 1,000 | -0% / +0% | +1% / +0% | +1% / +0% | +2% / +0% | +3% / -2% |
| YN | 300 | -2% / +0% | +4% / -7% | +16% / -20% | +7% / -21% | -19% / +0% |
| MYN | 300 | -3% / +1% | +4% / -7% | +18% / -22% | +4% / -15% | -21% / +8% |
| NG | 300 | -11% / +15% | -12% / +4% | -17% / +4% | -24% / +6% | -36% / +18% |
| LPB | 300 | -11% / +27% | -7% / +18% | -7% / +13% | -18% / +25% | -36% / +43% |
| Comeron | 300 | -0% / +0% | +5% / -10% | +8% / -12% | +1% / -12% | -30% / +14% |
| Li | 300 | -8% / +10% | -1% / -2% | +0% / -1% | -8% / -1% | -30% / +17% |
| GY | 300 | -0% / -1% | +1% / +0% | +3% / -0% | +4% / -1% | +4% / -3% |

The `GY` rows (maximum likelihood, added 2026-11-15) are from the same 30 combinations simulated again with 50
pairs per cell (`benchmarks/kaks/validate_gy.py grid`; YN on those pairs gives the numbers of its rows above to
within a few points): the estimator the simulation was checked with, now in COAL, without the bias of the YN00
approximation. Over the single combinations GY's cell medians stay within about 10 % up to dS 1 and within about
25 % at dS 2 to 3 for dS (cells of 50 pairs scatter that much), where YN's reach +95 %.

The medians hide the spread for the counting methods, whose error depends on the codon usage and kappa in both
directions. With skewed usage every counting method has dS 20 to 35 % too low and dN/dS 30 to 90 % too high
already at dS 0.03 to 0.1, rising to +100 to +190 % at dS 2 to 3 (kappa 4, omega 0.2, 300 codons: NG +31 % at
dS 0.1, +105 % at dS 2). With even usage and kappa 2 to 4, NG and LWL give dN/dS 10 to 30 % too low; with kappa 1,
NG is the most accurate method (within 3 % up to dS 1.5 at omega 0.2 and 1,000 codons) and LPB, MLPB and Li are within about 10 %. YN and MYN are
within about 5 % up to dS 0.5 in all 30 combinations; from dS about 1 they overestimate dS (+10 to +30 % with even
usage, +25 to +60 % skewed) and underestimate dN/dS (-20 to -40 % skewed), even at 1,000 codons. That is the bias
of the YN00 approximation, the same in yn00 and KaKs_Calculator; codeml's likelihood on the same pairs is within
8 %. Real genes usually lie between the two codon usages, so the counting methods' differences from YN on real
data (as in [tutorial 5](tutorials/05-choosing-a-method.md)) can be in either direction. YN's kappa is 3 to 10 %
high; MYN is as accurate as YN.

**Reliability limits.** YN and MYN have a median error within 15 % with at least 90 % of the pairs estimated up
to dS about 1 with even codon usage and about 0.5 with skewed usage, at every length from 100 codons. Status
`ok` does not mark the rest: at 1,000 codons, 29 % of the YN pairs at dS 1 and 35 % at dS 1.5 are more than 30 %
off, and none (dS 1) or 14 % (dS 1.5) of them are flagged. `saturated` becomes the majority only from dS 2.5;
`not_converged` is 3 to 7 % of the pairs from dS 1.5, at every length. Comeron's dN becomes NA when the synonymous
sites saturate (11 to 19 % of the pairs at dN 0.1 and dS 2, 41 % at dS 3), where YN still has a number.

**Gene length.** For 80 % of the YN dN/dS values to be within ±50 % of the truth, about 100 codons are needed at
dS 0.5 and about 300 at dS 0.1. At 30 codons it is 35 to 55 %, whatever the method; at dS 0.01, 3 to 57 % (300
codons) and 37 to 76 % (1,000 codons).

**Standard errors.** YN and MYN's SE of dS covers the truth (estimate ± 1.96 SE) in 88 to 95 % of the pairs for
dS 0.1 to 1 at 100 codons and more; beyond dS 1 it is too small (reported SE / actual SD 0.58 to 0.69 at dS 1 to
2, 1,000 codons), and also at dS 0.03 and below. The SE of dN covers less (0.58 to 0.75 at dS 1 to 2). NG's
standard errors match its scatter, but its bias brings the coverage down to 0.26 to 0.60 at 1,000 codons. Li's
standard errors are about 25 % too large.

**What follows for filtering** ("Filtering before interpretation" above): gene by gene, `--max-ds 1`, or 0.5 with
strong codon-usage bias, and 2 only as an outer bound; at least 100 codons and a dS of at least about 0.05 for a
ratio. YN stays the default (NG is better only with even codon usage and kappa near 1); Comeron's
`--comeron-min-sites 3` is supported (at 30 codons it cuts dN `NA` from 18 to 45 % to 1 to 2 % at dS 1 to 2 with
the same median error and RMSE; 10 is supported as well). The maximum-likelihood method `GY` does not have the
YN00 bias ("GY: maximum likelihood" below has its numbers on pairs simulated the same way).

**`--yn-weighting 0`** (equal pathway weights, one round: PAML yn00's default) is less biased than the default
weighting at omega below 1 and dS 1 to 1.5: kappa 1, 2 and 4, both codon usages, 1,000 codons, omega 0.2, dS 1:
dN/dS -8 to -32 % with the default, +5 to -22 % with `--yn-weighting 0`; and more pairs keep an estimate at dS 2
(70 % against 90 % over the 12 conditions). At omega 1 the two agree (the weights then hardly differ between pathways), and at omega
0.2 with even codon usage `--yn-weighting 0` overestimates dN/dS from dS 2 (+11 to +101 % at dS 2 to 3, against
-14 to +91 % with the default); 300 codons give the same picture. Measured on kappa 1, 2 and 4, both codon
usages, omega 0.2 and 1, dS 0.5 to 3, 300 and 1,000 codons, 200 pairs each. The default stays KaKs_Calculator's weighting: it is what makes `coal kaks -m YN` give KaKs_Calculator's
numbers (on the 4,000 reference pairs of the study, 0 of 3,120 differ; with `--yn-weighting 0`, 2,424 of 3,318
differ by more than 1e-5, the median dN/dS by 2 %).

**Reproducing it.** `benchmarks/kaks/sim_truth/` has the simulator (`codonsim.py`, `validate_sim.py`), the grid
(`gen.py`, fixed seeds; about 500 MB), the runs (`run_coal.sh`) and the tables and figures (`analyze.py`,
`perpair.py`, `make_tables.py`, `make_figs.py`); `run_coal.sh` lists the order. The comparisons with
KaKs_Calculator, yn00, codeml and orthologr of the study are not part of it.

## GY: maximum likelihood

`-m GY` estimates dN and dS by maximum likelihood under the codon model of Goldman & Yang (1994), the estimate
PAML's codeml gives for a pair of sequences (`runmode = -2`) and the model of KaKs_Calculator's `GY` (GY-HKY).
The approximate methods above count differences and correct them for multiple hits; GY fits the substitution
process itself, so it has no multiple-hit formula that fails and no YN00 approximation that drifts at high
divergence. It is the method to use when dS is above about 0.5 (see "Validation of GY" below), and it is slower:
about 15 to 25 times YN's time per pair.

```
coal dnds codon_aln.fasta -m GY --ml-columns -o dnds.tsv     # dN, dS, dN/dS, and kappa, t, lnL
coal kaks -i pairs.axt -m GY -o pairs.kaks                   # KaKs_Calculator's format (ML-Score, AICc, Model)
coal pipeline pairs --query a.fna --subject b.fna -m GY -o out
```

**Model.** Codon i changes to codon j, one nucleotide apart, at the rate pi_j, times kappa for a transition,
times omega for a nonsynonymous change; changes of two or three nucleotides at once and changes to or from a
stop codon of the genetic code have rate 0. The rate matrix is scaled to one expected substitution per codon,
so the branch length t between the two sequences is in substitutions per codon. The codon frequencies pi are
**F3x4**: the products of the base frequencies at the three codon positions, counted over both sequences in the
compared codons, with the stop codons removed (the frequencies YN uses; codeml's `CodonFreq = 2`, the value of
codeml's own control file). The log-likelihood is the sum over the compared codon columns of
log(pi_i P_ij(t)), and t, kappa and omega are the values that maximise it.

**What is reported.** dN and dS are those of Yang & Nielsen (2000), as codeml prints them: with rhoS and rhoN
the synonymous and nonsynonymous proportions of the rate at the estimates, and rhoS1, rhoN1 the same at
omega = 1, S = 3 L rhoS1 and N = 3 L rhoN1 are the sites (L codons), dS = t rhoS / (3 rhoS1) and
dN = t rhoN / (3 rhoN1). So dN/dS is omega. `coal dnds` prints dN, dS, dN/dS, the sites S and N, and as `Sd` and
`Nd` the expected differences S dS and N dN scaled to the observed number of nucleotide differences
(KaKs_Calculator's convention for its GY); `--ml-columns` adds kappa, t and lnL. `coal kaks` puts lnL in
ML-Score (see "Output").

**Standard errors.** `SE_dN` and `SE_dS` come from the observed information: minus the matrix of second
derivatives of the log-likelihood in (t, log kappa, log omega) at the maximum (taken by central differences of
the exact gradient) is inverted, and the variances of dS and dN follow by the delta method. This is what
codeml's `getSE = 1` reports, computed more accurately (below). They are given where the parameters are
estimated inside their bounds. A kappa or an omega below 1e-4 or above 99 is at its bound in all but name and is
held fixed; with omega there, the distance that the bound makes 0 has no standard error (`SE_dN` is `NA` for a
pair without nonsynonymous change, whose `SE_dS` is that of the synonymous changes alone). A saturated pair,
identical sequences and a pair with one kind of site have none. As for the other methods, they are standard
errors of dN and of dS, not of the ratio, and they describe the curvature at the maximum: for dS above about 2
the likelihood is skewed and a symmetric interval is a rough guide. In the simulations below the interval
dS ± 1.96 SE holds the true dS for 91 to 97 % of the pairs at every dS from 0.1 to 3 (95 % is the target), where
YN's interval holds it for 41 to 69 % from dS 2 with even codon usage.

**How it is computed.** The code is COAL's own (`cpp/src/kaks_gy.cpp`); only the symmetric eigen solver is
shared with YN. The rate matrix is symmetrised with pi and decomposed once per (kappa, omega); with the
decomposition the likelihood and its derivatives in t are cheap, so t is maximised out by a safeguarded Newton
iteration, and log kappa and log omega are found by a BFGS iteration with the exact gradient, inside the
bounds. Nothing is random: the same pair gives the same numbers in every run and with any number of threads,
and, because the likelihood depends only on the counts of unordered codon pairs, also when the two sequences
are swapped or the codon columns reordered (the same bytes at 17 digits; `tests/test_kaks_properties.py`).
The likelihood of a diverged pair can have two maxima, one at a moderate t and one on the "ridge" of saturation
(t at its bound, omega small); which is higher depends on the pair, and codeml reports whichever its random
start reaches. COAL starts a pair whose first maximum has t of 1 or more from three further points (one of them
on the ridge) and reports the highest maximum. On 4,388 *S. cerevisiae* / *C. glabrata* pairs the starts
disagreed for 22 pairs, all with t above 2. The ridge start takes its kappa and its omega t from where the first
start stopped. Where the first start had to be run a second time with lengthened steps ("Statuses", the row
`not_converged`), that is the second run's stop, which can lead the ridge start to a lower part of the ridge than
the first run's stop does (`GCTGATCCC` against `GATGCCGCG` under code 4: lnL -12.5851 from the one, -12.5748
from the other); such a pair gets a ridge start from each. Fixed starts are not a global search, and two kinds of
pair need more: on a pair of few codons the likelihood is flat and often has its maximum with kappa or omega at a
bound, where no start may lead; and a start that ends next to the plateau of saturation (t to infinity, where the
likelihood is the sum of log pi_a pi_b over the codon pairs and has no gradient in kappa or omega), on the ridge
or at a large t, cannot see a higher maximum elsewhere (`GTGGACAAC` against `TCCGGGAAG` under code 11: every start
ended on the ridge at lnL -17.3590 or below; the maximum is -16.5761 at t 7.96, kappa 1e-6). So after the starts
above (and the continuation of an unconverged one, "Statuses"), a pair of 40 codons or fewer, a pair whose result
did not converge, and a pair whose likelihood is within 0.05 per codon of the plateau is also scanned: the
likelihood, maximised over t, at 36 points of a grid in kappa and omega (1e-3, 0.01, 0.1, 1, 10, 99 each), and the
two highest points are started from (a pair of more than 40 codons only from a point higher than its result). What
the scan finds is taken only if the likelihood reported for it is higher than without it, after the
flat-likelihood rule ("Statuses"). Against an independent likelihood and a bounded global search, the pairs more
than 1e-3 below the global maximum, without and with the scan, in two sets of samples (this change's and its
audit's, each with its own search): of random pairs of 1 to 40 codons (nine codes), 67 and 4 of 600, 75 and 8 of
642; of 60 to 200 codons, 14 and 0 of 150, 11 and 2 of 150; of 300 to 800 codons, 5 and 1 of 60, 9 and 1 of 60.
Most of the remaining maxima have omega or kappa at a bound, on the ridge or near it, but not all (one pair of 39
codons has its maximum at kappa 2.0, omega 2.2, 0.055 above COAL's). No pair of more than 40 codons in the three
real data sets below is within 0.18 per codon of the plateau, so the scan does not touch them, and GY's time is
the same. COAL's likelihood is the highest of its starts, not a proven global maximum. In particular, t is
maximised as if the likelihood over t had one peak, and for a pair of very few codons it can have two, so the
reported t can be a local maximum (`GAGCTTTCT` against `AAGCTTTCT`, code 1: t 0.384 with lnL -12.765 is reported,
where t 10.12 has lnL -12.609, kappa and omega at 999 in both).

**Genetic codes.** `-c` sets the code for GY as for every method, and all 17 tables work (NCBI 1 to 6, 9 to 16,
21 to 23). The code enters the model in three places, each taken from the table: its stop codons are not states
of the chain (61 sense codons in the standard code, 60 in the vertebrate mitochondrial code, 63 in the ciliate
code) and are left out of the F3x4 frequencies before these are renormalised; a change is synonymous or
nonsynonymous by that code's translation; and a compared column with a stop codon of that code is skipped (TGA
is compared under codes 2 and 4, where it is Trp; AGA is skipped under code 2). Checked for every table
(`benchmarks/kaks/validate_gy.py codes`, 2026-11-15): on 40 random pairs per code the reported lnL is the
likelihood of the reported t, kappa and omega recomputed with numpy from the model's definition under that code
(to 1e-11 relative; S, dS and dN from the rates to the same; no parameter moved by 3 % has a higher likelihood),
and for the eleven tables codeml has (its `icode` 0 to 10: NCBI 1, 2, 3, 4, 5, 6, 9, 10, 12, 13, 15) COAL's lnL is
never below codeml's on 132 pairs, with omega, kappa and t the same to 3e-4 relative at most (codeml's printed
digits) in the 121 where both are at the same maximum inside codeml's bounds. Nothing had to be changed.

**Bounds.** t from 1e-6 to 50; kappa and omega from 1e-6 to 999 (codeml: t 1e-5 to 50, kappa 0.4 to 999, omega
0.001 to 99; KaKs_Calculator: t up to 3, kappa 0.01 to 30, omega 0.001 to 50). An omega at 1e-6 says that the
data support no nonsynonymous change (dN is then 1e-6 dS, where the counting methods print 0), and an omega at
999 that they support no synonymous change (dS is dN / 999; read dN/dS = 999 as "not estimable", as `NA` is for
the other methods when dS is 0). A parameter that the data want at a bound is reported at the bound itself: the
likelihood is linear in kappa or omega near 0, so its gradient in the logarithm fades before the bound is reached
and the iteration would stop about a factor 10 short of it (omega 1.1e-6 to 1.5e-5, up to 1e-5 log units below
the value at the bound, at a point that depends on the path). The same fading hides a maximum further inside from
an iteration that is near or on the bound. So when the iteration stops with kappa or omega in the zone of a bound
(below 1e-3 or above 99), both ends of the zone are tried: the bound, taken if the likelihood there is not lower,
and the inner edge (1e-3 or 99), taken if the likelihood there is higher; the iteration goes on from the better
one. An iteration that stops outside these zones is not affected, and later starts are seeded from where the
first one stopped by itself, never from a bound it was moved to.

**Which start wins.** Two log-likelihoods within 1e-7 of each other count as the same maximum (the optimiser's own
stopping rule leaves about 1e-11, starts from elsewhere up to 1e-8): the earlier, or the converged, start is kept.
A start is usable if it converged or ended on the ridge (t at its bound, converged or not: the pair is saturated
either way). A start that ran out of steps elsewhere never displaces a usable one by itself: if its likelihood is
higher than that of the best usable start, the iteration is continued from where it stopped (twice at most, and
not again once a continuation ends lower), and only if it still has not converged and is higher than the best
usable start by more than 1e-4 is the pair `not_converged`. (That last case has not occurred in 25,200 random
pairs; it is there so that a converged maximum known to be lower is not reported as the estimate.)

**Statuses.** The policy for gaps, ambiguity letters and stop codons is that of every method ("Policy").

| Situation | What GY gives | Status |
|---|---|---|
| Identical sequences | dN = dS = 0, t = 0, lnL = the sum of log pi_i; dN/dS, kappa and the site counts are `NA` (they have no estimate without a difference) | `ok` |
| The estimate of dS (or dN) is above 10 substitutions per site. Every pair whose t ends at its bound of 50 is such a pair (t / 3 is the mean of dS and dN over the sites) | that estimate and dN/dS are `NA`; the other estimate, the sites, kappa, t and lnL are kept. The note gives the estimate left out and says when t is at its bound | `saturated` |
| The likelihood is flat in t from the optimiser's stopping point to the bound (within 1e-6 log units): one codon, or every site different, so that the data put no upper limit on the distance | no estimate of either distance (`NA`); t is reported at its bound, 50. The note says so. (Without this rule the value would be wherever the optimiser's tolerance stopped it: 4.6 for `AAA` against `AAG`) | `saturated` |
| No synonymous site, or no nonsynonymous site (under F3x4: e.g. Met and Trp codons only) | omega cancels in the scaled model, so t and kappa are estimated and the one distance that exists is given, as for YN (if it is not saturated: the note then has both reasons) | `no_data` |
| The optimiser ends without a maximum, away from the bound of t | no estimate; the note has the last point. Before that: a start that uses up its 100 steps without converging, after steps of the gradient's size on a nearly flat likelihood, is run again with lengthened steps, and that run is taken if it converges with a likelihood not more than 1e-7 lower; a start that still has not converged is continued from its last point if it is the best one ("Which start wins"); `not_converged` is what remains: no usable start at all, or an unconverged point more than 1e-4 above the best usable start. Not seen on real and simulated genes; of 25,200 random pairs of 1 to 40 codons (all-random or heavily mutated, nine genetic codes), none (0.6.0: 103) | `not_converged` |

10 substitutions per site is far beyond what a pairwise comparison can resolve (the likelihood is flat there, and
an estimate of 15 or of 60 says the same thing: saturated); it is also the bound that the property tests put on
every method (P10). A dS of 3 to 10 with status `ok` is a legitimate maximum-likelihood estimate with a wide
error; `--max-ds` flags what you consider too high.

### Differences from codeml and KaKs_Calculator

1. **Deterministic and multi-start.** codeml draws its starting values at random (from the clock unless a seed
   is set) and stops at the maximum it reaches; COAL's starts are fixed, a diverged pair gets four (some
   five), and small or saturated pairs are scanned (above). On every pair of the comparison below COAL's
   likelihood is at least codeml's (lower by more than 1e-4 on none). That is what was measured, not a
   guarantee: both programs search locally, and on very short pairs a maximum can be missed (above).
2. **Wider bounds** (above). A pair with omega below 0.001 or kappa below 0.4 has those estimates in COAL and
   the bound in codeml, with a slightly higher likelihood in COAL.
3. **Saturation is a status.** codeml prints dS = 59.9 for a pair with t at 50; COAL prints `NA` and
   `saturated` for an estimate above 10. Together with point 1 this means that COAL can report `saturated`
   for a pair to which codeml gives a finite, moderate dS: where the likelihood has two maxima and the one on
   the ridge of saturation is the higher, COAL finds it and codeml, from most of its starts, stops at the lower,
   finite one (3 of the 160 yeast pairs below: codeml dS 7.7 to 9.4 with a likelihood 0.12 to 0.53 lower). The
   finite value is then not the maximum-likelihood estimate; the pair's lnL in `--ml-columns` can be compared
   with codeml's to see which case it is.
4. **Small probabilities are computed without cancellation.** The eigen decomposition gives each transition
   probability with an absolute error of about 1e-16, which is all of a probability that small: two codons two
   or three substitutions apart at a t below about 5e-4, as in thousands of nearly identical codons with one
   multiply changed codon. codeml then maximises rounding noise. One pair shows it: `t60000` of `e2.axt`
   (60,000 codons, one codon changed at all three positions). codeml reports lnL -246993.32 at its lower bound
   of t, but the likelihood recomputed independently at codeml's printed point is about -246999.7, and codeml's
   dN is 100 times below NG's. COAL replaces every pattern probability below 1e-8 by a series whose terms are
   all non-negative (uniformisation: exp(B tau) = e^(-x) sum x^k / k! R^k with R = I + B / lambda), likelihood,
   derivatives and gradient alike, and agrees with an independent matrix exponential to about 1e-12 relative
   there. For that pair COAL's maximum is lnL -246996.12 at t = 5.0e-5, with dN 2.4e-5 against NG's 2.0e-5.
   Probabilities above 1e-8 are computed as before.
5. **Identical sequences** have dN = dS = 0 and no kappa or omega, and lnL is the likelihood at t = 0, the sum
   of log pi over the codons (0.6.0 printed twice that); codeml prints the values its optimiser happened to
   stop at.
6. **Gaps and ambiguity.** COAL skips a codon column with a gap, an ambiguous base or a stop codon in either
   sequence and counts the F3x4 frequencies over the columns it compares; codeml's pairwise deletion
   (`cleandata = 0`) does the same for gaps and ambiguity, and refuses a file with a stop codon.
7. **Against KaKs_Calculator's GY**: the same model and definitions, but KaKs_Calculator bounds t at 3
   substitutions per codon (so its dS cannot exceed about 4), starts t at the proportion of differences and
   uses a random start for kappa and omega; its numbers were not reproduced (in the simulation study it took
   more than 20 minutes per file of 1,000-codon pairs and was run on 3 of 8 files). COAL writes `GY` in the
   Method column, KaKs_Calculator `GY-HKY` (accepted as a name by `-m`).
8. **Standard errors** differ from codeml's `getSE = 1` by 0.3 to 8 % in the median (the three data sets below;
   95 % of the pairs within 2 to 14 %, single pairs by more), and the difference is codeml's: it builds its
   Hessian from differences of the likelihood. COAL's were checked against a third computation (the likelihood
   rebuilt with numpy, the Hessian by five-point second differences, `tests/kaks/make_gy_se_reference.py`): the
   same to five digits on the 12 pairs of the test. Where omega is at codeml's bound (no nonsynonymous
   difference), codeml prints an SE of dS that is too small (0.0033 for a dS of 0.023 resting on four
   synonymous differences, where COAL gives 0.012).
9. **Not there:** F1x4 / F61 / Fequal frequencies, a fixed kappa or omega.

### Validation of GY

`benchmarks/kaks/validate_gy.py` compares `coal dnds -m GY` with codeml (PAML 4.10.10, `runmode = -2`,
`CodonFreq = 2`, one data set per pair) and with simulated truth (2026-11-15, a shared 2-core machine).

| Data | Pairs | COAL lnL below codeml's by more than 1e-4 | COAL higher by more than 1e-4 | Same maximum, codeml inside its bounds | Largest relative difference there: dS / dN/dS / kappa / t |
|---|---:|---:|---:|---:|---|
| Simulated (dS 0.1 to 3, omega 0.2 and 1, kappa 2 and 4, even and GC3-skewed usage, 300 codons) | 320 | 0 | 9 | 298 | 1.0e-3 / 5.3e-5 / 5.7e-6 / 4.6e-5 |
| *S. cerevisiae* / *C. glabrata* (median 417 codons) | 160 | 0 | 10 | 62 | 2.5e-4 / 4.7e-4 / 4.6e-5 / 2.5e-4 |
| *Salmonella* / *Shigella* (median 297 codons) | 159 | 0 | 1 | 149 | 7.6e-4 / 2.7e-4 / 5.3e-5 / 8.3e-5 |

codeml prints lnL with six decimals and the estimates with four or five, and the differences in the last
column are of that size (the medians are 2e-5 for dS, 5e-5 for dN/dS, 2e-6 for kappa); the median difference in
lnL is below 1e-6. COAL's likelihood is never the lower one. Where it is higher, codeml is at one of its bounds
(7 of the 9 simulated pairs: by 12 and 17 log units in the two where codeml ends with kappa at its bound of 0.4,
by less than 2 where its t is at 50; 66 of the yeast pairs have t = 50 in both programs, with the same
likelihood to 1e-4 in 60 of them) or stopped at a lower maximum (2 simulated pairs, by up to 0.13; 3 yeast
pairs, by up to 0.53; one *Salmonella* / *Shigella* pair, by 0.002). Which pairs these are changes from one
codeml run to the next, since its starting values are random. The golden test `tests/test_kaks_gy.py` holds 20 of these
pairs with codeml's output, three of them under the vertebrate mitochondrial code (`-c 2`, codeml `icode = 1`),
and 20 more pairs under the ten other tables codeml has; the tables it does not have (11, 14, 16, 21, 22, 23) are
tested against the likelihood recomputed from the definition.

**Against simulated truth** (the pairs of "Accuracy against simulated truth": GY94 with F3x4, truth in Yang &
Nielsen's definition; 300 codons, omega 0.2 and 1, kappa 2 and 4, 100 pairs per cell, 800 per column and codon
usage). Median relative error of dS / of dN/dS; in brackets the share of pairs with status `ok`, and the share
of those whose interval dS ± 1.96 SE holds the true dS (95 % if the standard errors are right):

| Method | usage | dS 0.1 | dS 0.5 | dS 1 | dS 1.5 | dS 2 | dS 2.5 | dS 3 |
|---|---|---|---|---|---|---|---|---|
| GY | even | -2% / +1% (100%; 92%) | +1% / +1% (100%; 96%) | +3% / +0% (100%; 97%) | +6% / -0% (100%; 96%) | +3% / -0% (94%; 94%) | -0% / +1% (83%; 94%) | +2% / +0% (75%; 91%) |
| YN | even | -3% / +2% (100%; 92%) | +0% / +2% (100%; 93%) | +2% / -2% (100%; 93%) | -3% / +6% (95%; 86%) | -17% / +23% (89%; 69%) | -28% / +41% (80%; 50%) | -37% / +62% (69%; 41%) |
| GY | skewed | +0% / +1% (100%; 93%) | +0% / +3% (100%; 94%) | +3% / -2% (100%; 94%) | +3% / +1% (99%; 94%) | +2% / +1% (100%; 94%) | +4% / -1% (93%; 95%) | +6% / -3% (88%; 94%) |
| YN | skewed | -1% / +0% (100%; 92%) | +2% / -10% (100%; 94%) | +27% / -31% (97%; 95%) | +22% / -35% (81%; 98%) | +4% / -27% (53%; 96%) | -12% / -28% (29%; 80%) | -22% / -27% (19%; 89%) |

GY's medians are within 6 % of the truth at every dS up to 3 with both codon usages, where YN is off by 20 to
40 % from dS 1 to 2; this is the "within 8 %" that the simulation study measured for codeml. What rises with dS
is the scatter, not the bias: at dS 3 and 300 codons, 18 % of the pairs are `saturated` (an estimate above
10), against YN's 56 % without an estimate. GY's standard errors are calibrated over the whole range.

**Speed** (CPU seconds per 1,000 pairs, `-j 1`, the same shared machine; codeml: wall time of one run with one
data set per pair):

| Data | GY | YN | codeml |
|---|---:|---:|---:|
| Simulated, 300 codons, dS 0.1 to 3 (1,600 pairs) | 32.5 | 1.6 | 82 |
| *Salmonella* / *Shigella*, median 297 codons (159 pairs) | 22.6 | 1.3 | 91 |
| *S. cerevisiae* / *C. glabrata*, median 417 codons, 97 of 160 saturated | 46.5 | 2.0 | 194 |

A pair costs 6 to 10 optimiser steps per start, each with one or two 61 x 61 symmetric eigendecompositions for
the likelihood and about as much work again for the gradient; a pair with t of 1 or more is started four times.
GY takes 15 to 25 times YN's time and is about 3 to 4 times faster than codeml on the same pairs (the standard
errors, six more decompositions, are included); with `-j` the pairs run in parallel as for every method.

**Reproducibility.** The optimiser stops at a tolerance (a gradient of 1e-5 per unit of log kappa or log omega),
so GY's numbers agree between platforms to about 1e-6 relative, not to the last bit as the counting methods do
(another C library's `exp` moves the last digits). The tests compare GY with tolerances; within one build the
output is the same bytes in every run. A pair at the edge of a status (an estimate of 9.99 or 10.01) can be
`ok` on one platform and `saturated` on another.

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
the same doubles when they use the same two-position pathways (`--comeron-twosubs`, whose default differs
between them: next section).

| | `Comeron`, default (`coal`) | `--comeron-flavour gestimator` (orthologr 0.4.3's gestimator) |
|---|---|---|
| Codons that differ at two positions | Comeron's two pathways (`--comeron-twosubs comeron`) | libsequence's branches (`--comeron-twosubs libsequence`) |
| Stop codon at the end of both rows | the column is dropped | compared like any codon: a stop codon has no sites (the column counts half the sites of the other codon), `TAA` against `TGA` or against a sense codon is counted as a change, and a pathway through a stop codon is kept, with Grantham distance 999 |
| Stop codon elsewhere | the column is skipped; status `warning` | compared, as above; no warning |
| `U` | read as `T` | the codon is not A/C/G/T: the column is skipped |
| Gap or ambiguous letter in a codon | column skipped | column skipped (the same) |
| A logarithm of the correction has an argument <= 0, in a class that has sites | the estimate is `NA`; status `saturated` (a class with fewer than 3 sites, `--comeron-min-sites`: the logarithm is 0, status `warning`; `NA`, `saturated` if the distance is then 0) | gestimator sets that logarithm to 0 and prints the finite number that results (too small). COAL gives that number, with status `warning` and the note `Ka: correction undefined, gestimator's value (...)` (or `Ks:`) |
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

### Two-position codon differences (`--comeron-twosubs`)

For codons that differ at two positions, libsequence's `TwoSubs` passes (codon 1, intermediate 1, codon 2,
intermediate 2) to a function whose parameters are (codon 1, intermediate 1, intermediate 2, codon 2). The four
branches it counts are therefore codon 1 to intermediate 1, intermediate 1 to intermediate 2, codon 1 to codon 2,
and codon 2 to intermediate 2, where Comeron's two pathways are codon 1 to intermediate 1 to codon 2 and codon 1
to intermediate 2 to codon 2; and `SingleSub`, given two codons that differ at two positions, classifies the
change at the last of them. The same nucleotide changes are counted, but the degeneracy class of some is taken
from another codon than the pathway's. The current libsequence (1.9.8) makes the same call, and its repository is
archived (read-only since 2023), so the call will not be changed there.

A consequence is that the estimate depends on which sequence is first: on 300 random 150-codon pairs with many
two-position differences, swapping the two rows changes Ka in 279 and Ks in 190 pairs (by up to 1.1%). So COAL
offers both:

| `--comeron-twosubs` | branches | default with | swap test (300 pairs) |
|---|---|---|---|
| `comeron` | the two pathways of Comeron (1995) | `--comeron-flavour coal` (the default) | 0 of 300 change |
| `libsequence` | libsequence's four branches, as above | `--comeron-flavour gestimator` | Ka 279, Ks 190 change |

`--comeron-flavour gestimator` alone still gives the numbers of orthologr 0.4.3's `gestimator()`;
`--comeron-flavour gestimator --comeron-twosubs comeron` gives the numbers of orthologr 0.5.0's native engine
(which corrected the same call). `coal pipeline` records the order a run used in `report.json`
(`parameters.comeron_twosubs`); coalr's argument is `comeron_twosubs` (`NULL`: by the flavour).

What it changes (default flavour, `coal pipeline pairs --preset pwalign -m Comeron`, 2026-11-14; new default
against `--comeron-twosubs libsequence --comeron-min-sites 0`, COAL's defaults before):

| Genome pair | Pairs | dS changed | dS median / p95 / max relative change | mean dS ratio new / old | dN/dS median / p95 change | status changes |
|---|---:|---:|---|---:|---|---:|
| *A. thaliana* / *A. lyrata* (orthologr's 1,000-gene example) | 995 | 474 of 994 | 0% / 2.9% / 12% | 1.006 | 0% / 5.4% | 0 |
| human / mouse | 946 | 809 of 932 | 1.0% / 4.8% / 50% | 1.015 | 2.1% / 8.2% | 6 |
| *S. cerevisiae* / *C. glabrata* | 4,164 | 807 of 847 | 5.9% / 39% / 166% | 1.108 | 9.9% / 36% | 314 (308 `ok` to `saturated`) |

Pairs of codons that differ at one or three positions are not affected. Against YN (a model-based
estimate, independent of both; measured with the gestimator flavour on the same alignments), the corrected dS is
slightly closer (median absolute difference 7.68 to 7.53% for
*A. thaliana*, 13.95 to 13.74% human / mouse, 25.97 to 25.47% yeast), and in a neutral simulation (kappa 2, 300
codons, 300 pairs) its median error of dN against the number of substitutions applied is 2.0% (libsequence's:
2.3%) and of dS -2.1% (-3.0%). (That truth is the count of substitutions per site as Comeron's own sites see them.
Against the model's dS in Yang & Nielsen's definition, "Accuracy against simulated truth" above, Comeron's dS at
kappa 2, omega 1 and 300 codons is +8 to +17% off with even codon usage and -17 to -21% with skewed usage; the
comparison of the two pathway orders is not affected.)

### Small site classes (`--comeron-min-sites`)

A class's correction needs 1 - 2P - Q > 0 and 1 - 2Q > 0 with P and Q per site of that class. The 2-fold class
whose synonymous change is a transversion (L2V) is tiny in short alignments (a median of 1.75 sites at 30 codons,
3 at 50, 6 at 100), and one change there makes its correction undefined. gestimator sets such a logarithm to 0;
the default flavour used to return `NA`. Now a class with fewer than `--comeron-min-sites` sites (default 3)
cannot make dN or dS `NA`: its undefined logarithm is set to 0, as gestimator does, and the row has status
`warning` and the note `Ka: a site class with fewer than 3 sites has an undefined correction (its logarithm set
to 0)` (or `Ks:`). A value that is not `NA` otherwise never changes; `--comeron-min-sites 0` restores the previous
rule. In the gestimator flavour the value is gestimator's either way, and the option only chooses the note (the
note above instead of `Ka: correction undefined, gestimator's value`), which tells orthologr that the row is not
saturated.

The rule applies only where the distance it gives is positive. In very short alignments, setting the small
class's logarithm to 0 can leave dS at exactly 0 although a synonymous change was seen (`ATGGCTAAA` /
`ATGGCCAAA`: one change, GCT to GCC; dS 0), and a 0 reads as "no synonymous change" (orthologr labels such a
pair `dS_zero`). Where a small class was the reason and the distance is not above 0, the row is what
`--comeron-min-sites 0` gives: `NA`, status `saturated`, in the default flavour; gestimator's 0 with the note
`Ks: correction undefined, gestimator's value (undefined logarithms set to 0)` (or `Ka:`) in the gestimator
flavour, which orthologr reads as `NA`. In 7,000 simulated neutral pairs of 5 to 100 codons (1,000 per length)
this concerns dS in 53 pairs at 5 codons, 4 at 10, 1 at 15 and none at 20 or more; dN never (none of its 168
small-class values is 0). No other value or status changes.

Simulated neutral pairs (kappa 2, 10-35% of codons differing, 600 pairs per length; the truth is the substitutions
applied): at 30 codons, dN is `NA` in 22 pairs (3.7%) with the previous rule and in none with the default; the
recovered values' mean absolute error is 8.8%, against 11.2% for all 30-codon pairs (50 codons: 2 pairs, 4.3%; 100
codons: 2 pairs; 300 codons: none). dS at 30 codons is `NA` in 58 pairs (9.7%) under every rule: that comes from
the pooled transversion term Bs, which has no class. Real data: *A. thaliana* 0 rows change, human / mouse 3,
yeast 14. Comeron (like every method) is unreliable below about 50-100 codons; no length cut-off is applied
(`--min-codons` flags short pairs).

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

**What an orthologr build with `engine = "coal"` has to do to offer Comeron and Li.** The fork's
`integration/all` (0.5.0) does all five (`R/coal.R`; its `dNdS_coal()` calls `coalr::coal_dnds_pairs()`), so
`dNdS(engine = "coal")` and `kaks_tool = "coal"` accept `"Comeron"` and `"Li"` there; the list is for other builds:

1. accept `dnds_est.method = "Comeron"` and `"Li"` for the COAL engine when
   `all(c("Comeron", "Li") %in% coalr::coal_info()$methods)` (as the fork's `coal_methods()` does);
2. pass `comeron_flavour = "gestimator"` and `li_flavour = "seqinr"` to `coalr::coal_dnds_pairs()` (and to
   `coal_dnds_clusters()`, `coal_kaks()`) for orthologr's numbers, with `terminal_stops = "compare"` only to
   reproduce 0.4.3, whose alignments include the stop codons; or leave the defaults for COAL's policy and say so;
3. saturation: with `comeron_flavour = "gestimator"` a row whose note contains `Ka: correction undefined` (or
   `Ks:`) holds the finite value of orthologr 0.4.3. The fork returns `NA` for these unless
   `options(orthologr.legacy_comeron_saturation = TRUE)`: set dN (dS) and dN/dS to `NA` for such rows otherwise;
4. Li: coalr already returns `NA` where seqinr returns 9.999999, and flags a negative value with status
   `warning` (orthologr turns it into `NA`);
5. values are doubles, where gestimator's file has 4 significant digits: `signif(x, 4)` gives orthologr's
   printed values;
6. a development coalr keeps its version (0.6.0.9000) when behaviour changes, so test
   `coalr::coal_info()$features` instead: `"comeron_twosubs"` and `"comeron_min_sites"` (the arguments) and
   `"comeron_min_sites_positive"` (the min-sites rule only where it gives a distance above 0, COAL #76; a coalr
   without it can give dS = 0 where the native engine gives `NA`). A coalr without `features` predates all three.

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

The YN and MYN times of this table predate the faster eigendecomposition of 2026-11-15 (four columns at a time,
neighbour-only sums, the AVX2 copy above): about 10 to 14% less CPU time (YN, MYN) on the AVX2 path on 2,862
*Salmonella* / *Shigella* pairs and 2 to 6% less on the baseline path (`-j 1`, a shared 2-core machine;
`benchmarks/perf/bench_pipeline.py`); 5 to 9% for a `coal pipeline pairs` run across three data sets; less on
other data. The numbers are the same bits.

`Comeron` and `Li` are table look-ups per codon column: 0.08 CPU s each for the 999 pairs on the same machine
(2026-11-13; about 13,000 pairs per CPU-second, the same as NG and the LWL family; gestimator's code compiled
on its own takes 0.09 s). Reading the pair for the two reference flavours (the columns with stop codons, the
counts of Kimura's distance) is part of preparing every pair and adds about 0.01 s per 1,000 pairs to the
other methods.

YN00 and MYN spend about 90 % of their time on P(t) (a 61 × 61 symmetric eigendecomposition per iteration, as in
both reference tools), so their single-thread gain is modest. P(t) is only read to weight the pathways between
two codons that differ at two or three positions, so **for a pair of sequences in which no codon differs at more
than one position it is not computed** (from 0.6). The result is the same bits: nothing that is computed depends
on the matrix then, and this was checked at 17 significant digits, before against after, on the 999 orthologr
pairs, the 600 simulated pairs, 4,000 + 3,000 random short pairs, 1,500 random pairs of 30 to 400 codons at
divergences from 0 to 0.3 per site and the three test sets of `tests/kaks`, for YN and MYN in both flavours and
both weightings and in both output formats (80 output files per build, byte-identical), with GCC 13 `-O3`,
Clang 18 `-O3` and, for four of the sets, GCC `-O2 -march=native` (2026-11-13); `tests/test_kaks_yn_shortcut.py`
repeats the comparison on every platform CI builds, against `--yn-always-pmat` (a testing option that computes
the matrix regardless). How much it saves depends on how close the sequences are: 141 of the 999 orthologr pairs
(*A. thaliana* / *A. lyrata*) have no such codon and none of the 600 simulated ones; on the 999 pairs YN took
about 10 to 20% less CPU time and MYN no measurable amount less (the machine's run-to-run noise was about 10%;
the pairs concerned are the short and nearly identical ones, which were cheap already), and on the 1,500 random
pairs, 913 of which have no such codon, YN took 46% and MYN 34% less (best of 9 runs each, `-j 1`). Closely
related genomes (strains, populations, recent duplicates) are where it pays. The gain comes from threads (`-j`): pairs are
independent and computed on the thread pool. In a pipeline, Ka/Ks takes a small share of the time (in the orthologr 0.4.3
profile of 2026-10-29, about 0.5 %); alignment dominates.

## Not done (yet)

- Model selection and averaging (MS, MA) and the gamma variants. For GY: codon frequencies other than F3x4
  (F1x4, F61) and a fixed kappa or omega.
- Comeron and Li for genetic codes other than the standard one (their references have no other; Comeron's
  tables could be built for any code, seqinr's site classes are written out for the standard code).
- Comeron with `maxHits` 1 or 2 (libsequence's option to leave out codons that differ at several positions;
  orthologr always uses 3).
- Standard input.
- A Python API: the commands are reachable from `coal`, and an in-process binding is planned with the R/Python
  bindings work.
- KaKs_Calculator 3.0 was not run (see Output).
- The textbook NG86 site count (see the note under "Methods"): `-m NG` is KaKs_Calculator's.
- A ThreadSanitizer job in the nightly checks: the race on `signgam` was found and its fix confirmed by hand
  (`coal-cpp kaks -j 4` in a `-fsanitize=thread` build: one report before, none after).
