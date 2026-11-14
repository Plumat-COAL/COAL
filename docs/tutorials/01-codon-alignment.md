# Tutorial 1: a codon alignment from a protein alignment and its coding DNA

**Goal:** take the coding sequences (CDS) of one gene from three bacteria, align their proteins, and thread the
DNA onto that alignment to get a codon alignment that PAL2NAL would also produce. Then break it on purpose in
the three ways that real data usually breaks it, and fix it.

| | |
|---|---|
| Data | glycerol kinase (*glpK*) of *Mycoplasmoides genitalium* G37, *M. pneumoniae* M129 and *M. gallisepticum* R(low) (NCBI RefSeq) |
| Needs | `coal` (any install for the codon alignment; the platform wheel for `coal translate` and `coal dnds`), Clustal Omega or MAFFT |
| Runtime | under 5 seconds for everything after the download (about 0.7 MB) |
| Tested with | COAL 0.6.0.dev0 (commit 009bdd5), Clustal Omega 1.2.4, PAL2NAL v14, Linux x86_64 |

These bacteria were chosen on purpose: mycoplasmas read **TGA as tryptophan** (NCBI genetic code 4), not as a
stop codon. That is the most common reason a codon alignment "does not match" (vertebrate mitochondria, code 2,
are the other classic), and it makes the errors below real ones rather than made-up ones.

## 0. Set up

From a checkout of COAL (or any folder: the download script needs `curl` or `wget` and `sha256sum` or `shasum`):

```
bash docs/tutorials/data/get_mycoplasma.sh          # -> coal-tutorial-data/, checksums verified
mkdir t01 && cd t01
```

```
ok                 coal-tutorial-data/mgenitalium_G37.cds.fna.gz
ok                 coal-tutorial-data/mpneumoniae_M129.cds.fna.gz
ok                 coal-tutorial-data/mgallisepticum_Rlow.cds.fna.gz
```

The files are NCBI's `cds_from_genomic.fna.gz` of assemblies GCF_000027325.1, GCF_000027345.1 and
GCF_000092585.1. Their SHA-256 checksums are in the script; a file that does not match is refused.

## 1. Pull out one gene

Each file holds every CDS of its genome. Take the *glpK* record of each and give it a short ID (the protein and
the DNA must end up with the **same IDs**; step 5 shows why):

```
D=../coal-tutorial-data
for sp in mgenitalium_G37:Mgen mpneumoniae_M129:Mpne mgallisepticum_Rlow:Mgal; do
  f=${sp%%:*}; s=${sp##*:}
  zcat $D/$f.cds.fna.gz | awk -v s=$s '/^>/{keep = /\[gene=glpK\]/; if (keep) print ">" s "_glpK"; next} keep'
done > glpK.cds.fasta
grep '>' glpK.cds.fasta
```

```
>Mgen_glpK
>Mpne_glpK
>Mgal_glpK
```

(On macOS use `gzcat` or `gunzip -c` instead of `zcat`.)

## 2. Translate and align the proteins

COAL needs a **protein alignment**; it does not align anything in this step. Translate with the right genetic
code and align with any aligner:

```
coal translate glpK.cds.fasta -c 4 -o glpK.faa
clustalo -i glpK.faa -o glpK.aln.fasta --threads 1
```

```
coal translate: 3 gene(s): 3 ok, 0 warning, 0 failed
```

`coal translate` drops the terminal stop codon and writes an internal stop as `*` (none here, because `-c 4`
reads TGA as W). `mafft glpK.faa > glpK.aln.fasta` works as well. If your proteins come from elsewhere (UniProt,
the annotation's `.faa`), that is fine too, as long as each one is the translation of its CDS.

## 3. Thread the DNA onto the alignment

```
coal glpK.aln.fasta glpK.cds.fasta -codontable 4 -o glpK.codon.fasta
head -4 glpK.codon.fasta
```

No message, exit status 0:

```
>Mgen_glpK
ATGGATCTAAAAAAACAATACATTATTGCCTTAGATGAAGGTACTAGTTCTTGTCGATCA
ATTGTTTTTGATCACAATCTTAACCAAATAGCAATAGCACAAAACGAATTTAACACTTTT
TTTCCTAATAGTGGTTGAGTTGAACAAGATCCACTAGAAATTTGATCAGCCCAACTAGCT
```

Each sequence is now 1,536 columns long (512 codon columns = the 512 columns of the protein alignment); a gap
in the protein alignment became `---`. The terminal stop codon of the CDS is not in the output, because it is
not in the protein alignment.

Other output formats, as PAL2NAL writes them:

```
coal glpK.aln.fasta glpK.cds.fasta -codontable 4 -output paml -o glpK.paml       # codeml's input
coal glpK.aln.fasta glpK.cds.fasta -codontable 4 -nogap -o glpK.nogap.fasta      # drop columns with a gap
```

```
$ head -c 120 glpK.paml
   3   1536
Mgen_glpK
ATGGATCTAAAAAAACAATACATTATTGCCTTAGATGAAGGTACTAGTTCTTGTCGATCA
```

`-nogap` leaves 1,512 columns: the 8 codon columns with a gap in any of the three sequences (4 in each) are gone.

### Same output as PAL2NAL

`pal2nal` is COAL's drop-in for `pal2nal.pl`: the same command line, the same defaults (CLUSTAL output). With
PAL2NAL v14 installed, the two outputs are byte-identical:

```
pal2nal    glpK.aln.fasta glpK.cds.fasta -codontable 4 > coal_pal2nal.clustal
pal2nal.pl glpK.aln.fasta glpK.cds.fasta -codontable 4 > real_pal2nal.clustal
cmp coal_pal2nal.clustal real_pal2nal.clustal && echo IDENTICAL
```

```
IDENTICAL
```

The few places where COAL deliberately differs from PAL2NAL are listed in [differences.md](../differences.md).

## 4. Error 1: the wrong genetic code (and the hint)

Leave out `-codontable 4`, as one easily does:

```
coal glpK.aln.fasta glpK.cds.fasta -o wrong_code.fasta
echo "exit $?"
```

```
coal: notice: wrong genetic code? Under genetic code 4 (Mold, protozoan and coelenterate mitochondrial; Mycoplasma), all 3 sequences match their protein; under code 1, none do. Re-run with -codontable 4.
coal: warning: Mgen_glpK pepAlnPos 46: W does not correspond to TGA
coal: warning: Mgen_glpK pepAlnPos 55: W does not correspond to TGA
coal: warning: Mgen_glpK pepAlnPos 91: W does not correspond to TGA
...
coal: warning: Mpne_glpK pepAlnPos 322: W does not correspond to TGA
coal: warning: ... and 14 more warnings (see them all with --max-warnings 0 or --report FILE.json)
coal: summary: 3 of 3 sequences have codons that do not match their amino acid (kept in the alignment). Under genetic code 4 (Mold, protozoan and coelenterate mitochondrial; Mycoplasma), all 3 sequences match their protein; under code 1, none do. Re-run with -codontable 4.
exit 0
```

How to read it:

- **The notice comes first and names the fix**: COAL re-threaded the genes under every code it supports and
  found one under which all of them fit. Every warning below it (`W does not correspond to TGA`) is the same
  problem, one codon at a time.
- **The exit status is 0 and the file is written**, as PAL2NAL does: a few mismatching codons are often real
  (sequencing errors, selenocysteine). The alignment is still the right one here, column for column; what would
  go wrong is everything downstream that uses code 1 (Ka/Ks under code 1 reads the 10 to 13 TGA codons of each gene as internal stop codons and skips those columns).
- `--report run.json` has the same hint in machine-readable form (`hints[0].suggested_codontable: 4`) and every
  mismatch per sequence.
- `--strict` does **not** stop this run: it refuses genes where more than half of the checkable codons
  disagree (the wrong DNA), and here 13 of 508 do. Check the exit status *and* the notices, or read the report.

The fix is the one the notice gives: `-codontable 4`.

## 5. Error 2: IDs that do not match (pairing by order)

A protein alignment from one source and a CDS file from another rarely share IDs. Make such a file, with NCBI's
own IDs and in another order:

```
for f in mgallisepticum_Rlow mgenitalium_G37 mpneumoniae_M129; do
  zcat $D/$f.cds.fna.gz | awk '/^>/{keep = /\[gene=glpK\]/; if (keep) print $1; next} keep'
done > glpK.cds.ncbi.fasta
coal glpK.aln.fasta glpK.cds.ncbi.fasta -codontable 4 -o by_order.fasta
echo "exit $?"
```

```
coal: error: Mgal_glpK: protein and DNA are inconsistent (no match found) (paired by file order with DNA 'lcl|NC_000912.1_cds_WP_010874407.1_50', because no protein ID matches a DNA ID)
coal: summary: 1 of 3 sequences has codons that do not match their amino acid (kept in the alignment).
coal: summary: 1 sequence of 3 was written to the output but is probably wrong (first: Mgen_glpK: 479 of the 508 codons COAL could check do not match this protein, so this is almost certainly the wrong DNA for it (or the wrong genetic code). ...) Remove it before analysis, or use --strict to refuse such genes.
coal: summary: Mgal_glpK: The DNA does not spell the start of this protein in any reading frame.
coal: summary: The run stopped at Mgal_glpK, so no output was written. Add --keep-going to write the other sequences and leave this one out.
exit 1
```

How to read it:

- **"paired by file order ... because no protein ID matches a DNA ID"** is the cause. Like PAL2NAL, COAL pairs
  by ID when every protein ID is in the DNA file and by position otherwise, so the first protein got the
  *M. gallisepticum* DNA.
- One gene could not be placed at all (`no match found`), so the run failed with exit 1 and wrote nothing.
- One gene *could* be placed, codon by position, with 479 of 508 codons wrong. PAL2NAL would write that row
  without a word; COAL flags it, and with `--strict` refuses it:

  ```
  coal glpK.aln.fasta glpK.cds.ncbi.fasta -codontable 4 -o by_order.fasta --strict
  ```
  ```
  coal: error: Mgen_glpK: protein and DNA are inconsistent (479 of 508 checkable codons do not match the protein) (paired by file order with DNA 'lcl|NC_004829.2_cds_WP_011113244.1_17', because no protein ID matches a DNA ID)
  ...
  exit 1
  ```

The fix: give the DNA the protein alignment's IDs (step 1 does that), or put the records in the same order.
COAL also pairs IDs that differ only in a version suffix (`.1`), upper/lower case or `sp|...|` fields, and says
so; anything else needs renaming.

## 6. Error 3: mismatching codons in a gene that is otherwise right

Single mismatching codons, as opposed to a whole gene, usually mean that the protein is not exactly the
translation of this CDS: another assembly or strain, a UniProt sequence corrected for a frameshift, a
selenocysteine (U, from TGA), or an `X` in the protein. COAL then falls back to PAL2NAL's placement (10-residue
chunks) and prints one warning per codon that disagrees: `Mgen_glpK pepAlnPos 46: W does not correspond to TGA`
names the sequence, the **alignment column** (1-based) and what disagrees. What to do:

| Warnings | Likely cause | Fix |
|---|---|---|
| the same residue / codon pair over and over (W / TGA, M / ATA, stops / AGA) | wrong genetic code | the code the notice names (`-codontable N`) |
| most of one gene's codons, with a "probably wrong" notice | wrong DNA for that protein (IDs, order, isoform) | fix the pairing; `--strict` to refuse such genes in a pipeline |
| a handful of scattered codons | protein and DNA from different sources or versions | translate the CDS yourself (step 2) and realign, or `-nomismatch` to drop those columns |
| "does not spell the start of this protein in any reading frame" | the DNA is not this gene's, or is truncated | check the record; `--keep-going` writes the others and leaves it out (exit 2) |

Extra sequence around the gene (a 5' UTR, the stop codon, a poly-A tail) is **not** an error: COAL finds where
the protein starts in the DNA, in any reading frame.

## 7. A first look at dN/dS

With the platform wheel (`coal-cpp`), the codon alignment goes straight into Ka/Ks for every pair of its rows:

```
coal dnds glpK.codon.fasta -c 4 -j 1 | cut -f1-6,13-15
```

```
coal dnds: 3 pair(s), 1 method(s): 0 ok, 0 warning, 0 not_converged, 3 saturated, 0 no_data, 0 error
seq1	seq2	method	dN	dS	dNdS	codons_compared	status	note
Mgen_glpK	Mpne_glpK	YN	0.17381333	NA	NA	508	saturated	Ks saturated (multiple-hit correction undefined)
Mgen_glpK	Mgal_glpK	YN	0.27543137	NA	NA	504	saturated	Ks saturated (multiple-hit correction undefined)
Mpne_glpK	Mgal_glpK	YN	0.25845413	NA	NA	504	saturated	Ks saturated (multiple-hit correction undefined)
```

How to read it: dN is estimated (0.17 to 0.28 substitutions per nonsynonymous site), but **dS is not**: these
species have diverged so far at synonymous sites that the correction for multiple substitutions has no value
(`saturated`). That is a result, not a bug, and it means that dN/dS cannot be measured between these genes by a
pairwise method. Without `-c 4` the same command also skips 13 to 15 "internal stop" columns per pair: one more
reason to get the code right first.

[Tutorial 5](05-choosing-a-method.md) shows why other methods (`-m NG`) do print a dS here, and why that number
should not be trusted either. [Tutorial 2](02-pairs-dnds.md) does all of this for whole genomes in one command.

## Clean up

```
cd .. && rm -rf t01
```

## What to cite

`coal cite` prints the references for what you ran. For this tutorial: COAL, and PAL2NAL (Suyama et al. 2006),
whose method the codon alignment follows; Clustal Omega (Sievers et al. 2011) for the protein alignment; and for
step 7, Yang & Nielsen (2000) and KaKs_Calculator 2.0 (Wang et al. 2010).
