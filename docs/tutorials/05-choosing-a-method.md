# Tutorial 5: choosing a Ka/Ks method, and knowing when none will do

**Goal:** see, on the pairs of [tutorial 2](02-pairs-dnds.md), how much the choice between YN, NG, LWL, MYN and
Comeron changes dN and dS, where the methods agree and where they part, and what saturation and short genes do to
any of them. Everything here is grounded in [kaks.md](../kaks.md), which has the validation behind each claim.

| | |
|---|---|
| Data | the 3,142 *E. coli* K-12 / *S. enterica* LT2 ortholog pairs of tutorial 2 (`ecsa_filtered/`); the *glpK* alignment of tutorial 1 (`t01/`); run from the scratch folder `coal-tutorials/` |
| Needs | the platform wheel (`coal-cpp`), Python 3 for the comparison script |
| Runtime | 16 s for five methods on 3,142 pairs with one thread (YN and MYN take most of it) |
| Tested with | COAL 0.6.0.dev0 (commit 009bdd5), Linux x86_64 |

## The short answer

| Situation | Use | Why |
|---|---|---|
| Anything, unless you have a reason | **YN** (the default) | models transition/transversion bias and codon frequencies; in COAL's simulations unbiased up to dS about 0.5, then dN/dS 10 to 40% too low by dS 1 to 2; standard errors; the numbers of KaKs_Calculator and PAML yn00 |
| Transition rates that differ between purines and pyrimidines, or a second opinion on YN | **MYN** | YN with TN93; the same accuracy in the simulations; flagged `saturated` slightly less often and `not_converged` slightly more often |
| Reproducing or continuing an orthologr analysis | **Comeron** (orthologr's default), with `--comeron-flavour gestimator` for orthologr 0.4.3's exact numbers | same numbers as orthologr's gestimator; no standard errors |
| Comparing with older work that used them | **NG**, **LWL**, **LPB** (and MLWL, MLPB) | counting methods: fast, simple; they count sites without codon frequencies, so codon-usage bias biases their dN/dS, in either direction (with skewed usage, as in section 1, too high) |
| Positive selection at sites or on branches | none of these | pairwise whole-gene dN/dS averages over sites and time; use codon models (codeml, HyPhy) on the codon alignments |

Whatever you choose, filter first ([tutorial 2](02-pairs-dnds.md#5-filter-then-interpret)), use one method for
everything you compare, and name it (`coal cite` gives the references).

## 1. Five methods on the same pairs

Take the codon alignments that tutorial 2 kept (`--keep-alignments` in its step 5) and compute five methods with
`coal dnds`; records 1+2, 3+4, ... are the pairs:

```bash
coal dnds ecsa_filtered/codon_alignments.fasta --pairs consecutive \
    -m YN -m NG -m LWL -m MYN -m Comeron -j 1 -o methods.tsv
```

```text
coal dnds: 3142 pair(s), 5 method(s): 14623 ok, 20 warning, 23 not_converged, 1044 saturated, 0 no_data, 0 error
coal dnds: warning: 23 result(s): the YN/MYN iteration did not converge, so there is no estimate: Ka and Ks are NA (status not_converged; the note has the last round)
```

(The same as `coal pipeline pairs ... --pairs ecoli_vs_salmonella/orthologs.tsv -m YN -m NG -m LWL -m MYN -m
Comeron`: the numbers are identical, because both routes use the same engine on the same alignments.)

```bash
python3 $T/scripts/compare_methods.py methods.tsv
```

```text
statuses per method:
  method              ok       warning     saturated not_converged       no_data         error
  NG                2924             3           215             0             0             0
  LWL               2998             3           141             0             0             0
  YN                2973             3           155            11             0             0
  MYN               2853             3           274            12             0             0
  Comeron           2875             8           259             0             0             0

2694 pairs have a dN/dS from every method; median dS / median dN/dS, by YN's dS:
  YN dS      pairs                NG               LWL                YN               MYN           Comeron
  0-0.5        222    0.331 / 0.0418    0.343 / 0.0410    0.294 / 0.0452    0.289 / 0.0445    0.301 / 0.0478
  0.5-1        638    0.809 / 0.0390    0.793 / 0.0408    0.795 / 0.0371    0.872 / 0.0332    0.710 / 0.0411
  1-1.5        985    1.216 / 0.0444    1.123 / 0.0507    1.280 / 0.0408    1.436 / 0.0355    1.043 / 0.0489
  1.5-2        538    1.525 / 0.0486    1.321 / 0.0599    1.686 / 0.0410    1.778 / 0.0380    1.303 / 0.0535
  >= 2         311    1.504 / 0.0480    1.303 / 0.0607    2.381 / 0.0267    1.722 / 0.0382    1.232 / 0.0549

YN, pairs with status ok and 0.5 <= dS <= 1.5: median relative standard error, by compared codons:
  0-100     codons:   107 pairs   SE_dS/dS 0.31   SE_dN/dN 0.32
  100-200   codons:   342 pairs   SE_dS/dS 0.23   SE_dN/dN 0.25
  200-400   codons:   711 pairs   SE_dS/dS 0.17   SE_dN/dN 0.20
  >= 400    codons:   487 pairs   SE_dS/dS 0.13   SE_dN/dN 0.16
```

How to read it:

- **At low divergence the methods agree.** Below a YN dS of 0.5, the five medians of dS lie between 0.29 and
  0.34 and those of dN/dS between 0.041 and 0.048. For close species any method gives the same picture.
- **As dS grows, they part, and always in the same direction.** At YN dS 1.5 to 2, Comeron's median dS is 1.30
  and LWL's 1.32 against YN's 1.69 and MYN's 1.78; their dN/dS is accordingly 30 to 45% higher (0.054 and 0.060
  against 0.041). The counting methods (NG, LWL, Comeron) count sites without the codon
  frequencies, so codon-usage bias shifts their dN/dS; with skewed usage, as here, upwards. In COAL's simulations
  the gap is there already at dS 0.03, so it is not only a weaker multiple-hit correction. YN and MYN model codon
  usage and the transition/transversion bias, but they are not exact at this distance either: they recover the
  true omega within about 10% only up to dS about 0.5, and at dS 1 to 2 their dN/dS is 10 to 40% too low. At YN
  dS 1.5 to 2 the true value most likely lies between the two groups, and no method here is accurate
  ([kaks.md, "accuracy vs simulated truth"](../kaks.md)).
- **Above YN dS 2, nobody is right.** YN's median dS jumps to 2.38 while the others stay at 1.2 to 1.7: the
  counting methods cannot exceed their saturation ceiling, YN can, and neither number is reliable. These are the
  pairs that tutorial 2 filters out with `dS <= 2`. That is an outer bound only: for conclusions about single
  genes, cut at dS 1 (`--max-ds 1`; about 0.5 with strong codon-usage bias).
- **The statuses differ too.** For the same 3,142 pairs, from 141 (LWL) to 274 (MYN) are `saturated`. A method
  that saturates less is not more accurate there: it gives a number where the correction is about to break down.
  Only YN and MYN iterate, so only they have `not_converged`.
- **Short genes are noisy whatever the method.** For pairs of the same divergence, the relative standard error of
  YN's dS falls from 31% under 100 codons to 13% over 400. A gene of 80 codons with dN/dS 0.06 ± a third is
  compatible with a wide range of selective regimes; this is why 50 codons is the bare minimum (about 100 for ratios gene by gene) and why single
  short genes with dN/dS > 1 (tutorial 2's *yqgC*) are not evidence of anything.

Choose the method before looking at the result, and do not switch methods to make a gene cross 1.

## 2. Saturation: when no pairwise method will do

Tutorial 1's *glpK* of three mycoplasmas, with YN and with NG / LWL:

```bash
cd t01
coal dnds glpK.codon.fasta -c 4 -j 1 | cut -f1-6,13-15
coal dnds glpK.codon.fasta -c 4 -j 1 -m NG -m LWL | cut -f1-6,13-15
cd ..
```

```text
seq1	seq2	method	dN	dS	dNdS	codons_compared	status	note
Mgen_glpK	Mpne_glpK	YN	0.17381333	NA	NA	508	saturated	Ks saturated (multiple-hit correction undefined)
Mgen_glpK	Mgal_glpK	YN	0.27543137	NA	NA	504	saturated	Ks saturated (multiple-hit correction undefined)
Mpne_glpK	Mgal_glpK	YN	0.25845413	NA	NA	504	saturated	Ks saturated (multiple-hit correction undefined)
```

```text
seq1	seq2	method	dN	dS	dNdS	codons_compared	status	note
Mgen_glpK	Mpne_glpK	NG	0.19610998	1.8710231	0.1048143	508	ok	
Mgen_glpK	Mpne_glpK	LWL	0.2043005	1.4556418	0.14035081	508	ok	
Mgen_glpK	Mgal_glpK	NG	0.33428012	1.4539696	0.22990861	504	ok	
Mgen_glpK	Mgal_glpK	LWL	0.34348963	NA	NA	504	saturated	Ks saturated (multiple-hit correction undefined)
Mpne_glpK	Mgal_glpK	NG	0.29610777	1.8100969	0.16358669	504	ok	
Mpne_glpK	Mgal_glpK	LWL	0.29800271	1.5917973	0.18721147	504	ok	
```

YN says the synonymous sites are saturated. NG prints dS of 1.45 to 1.87 with status `ok`, and LWL saturates
on one pair and not on the others. The NG numbers are not wrong in the arithmetic, but **`ok` is not a verdict on
reliability**: these genes are AT-rich and far apart, the synonymous sites have turned over more than once, and
NG's Jukes-Cantor correction has not run out of room yet. One *M. gallisepticum* pair gets a smaller NG dS
(1.45) than the much closer *M. genitalium* / *M. pneumoniae* pair (1.87; compare their dN, 0.33 and 0.20): a sign
that the numbers no longer measure distance. Treat dS above about 1 (about 0.5 with strong codon-usage bias) as unreliable and above 2 as saturated,
whatever the status says ([kaks.md, "Filtering before interpretation"](../kaks.md#filtering-before-interpretation)),
and for such distances:

- report dN and protein identity, which are still measurable;
- compare closer species if dS is what you need ([tutorial 3](03-clusters.md#3-split-the-pairs-by-species-before-reading-them)
  shows the effect of the outgroup);
- for dN/dS across deep divergence, fit codon models on a tree of several sequences (PAML codeml, HyPhy), with
  COAL's codon alignments as input (`-output paml`, [tutorial 1](01-codon-alignment.md#3-thread-the-dna-onto-the-alignment)).

## 3. Very close sequences

The opposite end has its own trap. Between the two *C. trachomatis* strains of tutorial 3 (median dS 0.0098), dS
rests on a few synonymous differences per gene and dN/dS is high (median 0.29) for reasons that have nothing to do
with the method: slightly deleterious variation that selection has not yet removed. The methods agree at low
divergence (section 1), so changing the method does not change this. `dNdS` is `NA` only when dS is below 1e-6; a floor on dS (about 0.05 for ratios gene by
gene; tutorial 2's 0.01 removes only the nearly identical pairs), or comparing like with like (strains with
strains, species with species), is the remedy.

## 4. Reproducing other programs

| To match | Options | Notes |
|---|---|---|
| KaKs_Calculator 2.0 (`-m YN`, `NG`, `LWL`, `MLWL`, `LPB`, `MLPB`, `MYN`) | the same `-m`; `coal kaks` writes its 22-column format | to 1e-9 on the validation sets, apart from documented edge cases (kaks.md, "Differences from KaKs_Calculator"). KaKs_Calculator's default, MA, and GY are not implemented |
| PAML yn00 | `-m YN --yn-flavour paml`; `--yn-weighting 0` for yn00's distributed `weighting = 0` | KaKs_Calculator's YN and yn00 (`weighting = 1`) give the same numbers on the validation sets; they differ in edge cases of kappa, and in the pathway weighting of yn00's default control file |
| orthologr's `"Comeron"` (gestimator) | `-m Comeron --comeron-flavour gestimator` (plus `--terminal-stops compare` in the pipeline for orthologr 0.4.3) | standard genetic code only |
| orthologr's `"Li"` (`seqinr::kaks`) | `-m Li --li-flavour seqinr` | `Li` is not `LPB`: same paper, different counting |
| textbook NG86 (MEGA, PAML's NG86 output) | not available | `-m NG` is KaKs_Calculator's NG; its dS differs by up to about 0.02 from the textbook count: do not mix the two |

## 5. Cite what you used

```bash
coal cite YN NG LWL MYN Comeron
```

```text
[coal] Newton C, Murati K (2026). COAL: codon alignments from protein alignments, version 0.6.0.dev0. https://github.com/plumat-coal/COAL
[suyama2006] Suyama M, Torrents D, Bork P (2006). PAL2NAL: robust conversion of protein sequence alignments into the corresponding codon alignments. Nucleic Acids Research 34:W609-W612. https://doi.org/10.1093/nar/gkl315
[yang2000] Yang Z, Nielsen R (2000). Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Molecular Biology and Evolution 17:32-43. https://doi.org/10.1093/oxfordjournals.molbev.a026236
[wang2010] Wang D, Zhang Y, Zhang Z, Zhu J, Yu J (2010). KaKs_Calculator 2.0: a toolkit incorporating gamma-series methods and sliding window strategies. Genomics, Proteomics & Bioinformatics 8:77-80. https://doi.org/10.1016/S1672-0229(10)60008-3
[yang2007] Yang Z (2007). PAML 4: phylogenetic analysis by maximum likelihood. Molecular Biology and Evolution 24:1586-1591. https://doi.org/10.1093/molbev/msm088
[nei1986] Nei M, Gojobori T (1986). Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. Molecular Biology and Evolution 3:418-426. https://doi.org/10.1093/oxfordjournals.molbev.a040410
[li1985] Li W-H, Wu C-I, Luo C-C (1985). A new method for estimating synonymous and nonsynonymous rates of nucleotide substitution considering the relative likelihood of nucleotide and codon changes. Molecular Biology and Evolution 2:150-174. https://doi.org/10.1093/oxfordjournals.molbev.a040343
[zhang2006] Zhang Z, Li J, Yu J (2006). Computing Ka and Ks with a consideration of unequal transitional substitutions. BMC Evolutionary Biology 6:44. https://doi.org/10.1186/1471-2148-6-44
[comeron1995] Comeron JM (1995). A method for estimating the numbers of synonymous and nonsynonymous substitutions per site. Journal of Molecular Evolution 41:1152-1159. https://doi.org/10.1007/BF00173196
[thornton2003] Thornton K (2003). libsequence: a C++ class library for evolutionary genetic analysis. Bioinformatics 19:2325-2327. https://doi.org/10.1093/bioinformatics/btg316
```

`coal cite` without arguments prints the whole table of what to cite for what; the README has it too.

## Clean up

```bash
rm -f methods.tsv
```
