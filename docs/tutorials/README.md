# COAL tutorials

Worked examples on real, small, public data: from coding sequences to an interpreted dN/dS table. Every command
was run as written, and the outputs shown are the real ones (COAL 0.6.0.dev0, commit 009bdd5, Linux x86_64).

| # | Tutorial | Data | Runtime (1 thread) |
|---|---|---|---|
| 1 | [A codon alignment from a protein alignment and its CDS](01-codon-alignment.md): PAL2NAL-compatible output, the wrong genetic code and the hint, IDs that do not match, mismatching codons | *glpK* of three *Mycoplasmoides* (genetic code 4) | seconds |
| 2 | [Two genomes in, a dN/dS table out](02-pairs-dnds.md): `coal pipeline pairs`, statuses, saturation, filtering, `report.json`, reproducing a run | *E. coli* K-12 and *S. enterica* Typhimurium LT2 | 16 s |
| 3 | [Gene families across genomes](03-clusters.md): `coal pipeline clusters`, the per-cluster summary, species pairs, tighter clusters, large clusters, picking the aligner | four *Chlamydia* genomes | 85 s per run |
| 4 | [The same from R](04-r-coalr.md): `coalr` in process and through Python, orthologr's `dNdS(engine = "coal")` | as 2 and 3 | 7 to 71 s |
| 5 | [Choosing a method](05-choosing-a-method.md): YN, NG, LWL, MYN, Comeron; saturation, short genes, very close sequences | as 1 and 2 | 16 s |

Start with 1 if you have a protein alignment already, with 2 if you have genomes. 3 and 4 build on 2; 5 uses the
outputs of 1 and 2.

## Setup

1. **COAL with its C++ engine.** Tutorial 1's codon alignment works with any install; everything else needs
   `coal-cpp`: the platform wheel of a release, or a source build with `COAL_CPP` set ([install.md](../install.md)).
   `coal --cpp-path` prints where it is.
2. **External programs:** DIAMOND 2.1.0 or later (tutorials 2-4), Clustal Omega (3, 4), MAFFT (one step of 3).
   They can be on `PATH`, or installed into COAL's cache, checksum-verified:

   ```
   coal pipeline install diamond clustalo mafft
   ```
3. **The scripts of this folder.** The tutorials call the download scripts and two small Python scripts
   (standard library only) through `$T`:

   ```
   git clone https://github.com/plumat-coal/COAL
   export T=$PWD/COAL/docs/tutorials
   mkdir coal-tutorials && cd coal-tutorials      # a scratch folder for data and results
   ```
4. **R** (tutorial 4 only): R 4.1 or later, a C++17 compiler, and the packages named there.

On macOS, use `gzcat` (or `gunzip -c`) where the tutorials say `zcat`, and `shasum -a 256` for `sha256sum`.

## Data

All data are NCBI RefSeq genome annotations (`*_cds_from_genomic.fna.gz`), downloaded by the scripts in
[`data/`](data/) from `https://ftp.ncbi.nlm.nih.gov/genomes/all/` and checked against fixed SHA-256 checksums.
Nothing is committed to the repository. A RefSeq assembly version never changes once released, so the checksums
stay valid; a file that does not match is refused.

| Script | Organism | Assembly | CDS | Size (gz) | SHA-256 |
|---|---|---|---:|---:|---|
| `get_mycoplasma.sh` | *Mycoplasmoides genitalium* G37 | GCF_000027325.1 (ASM2732v1) | 524 | 187 kB | `e7e499e12911e044fdb642a9056e64cf2a173e9837c867c492b99fa2abb7a24d` |
| | *Mycoplasmoides pneumoniae* M129 | GCF_000027345.1 (ASM2734v1) | 724 | 256 kB | `1ed0e736c2a3d163ef072aa5ffd46364718ef7e1427a231431ecd72b569ef6da` |
| | *Mycoplasmoides gallisepticum* R(low) | GCF_000092585.1 (ASM9258v1) | 775 | 295 kB | `6e1e66665c4cf6d5169fb9a309400ec58cce6c4aab591843226f43abffdb04af` |
| `get_enterobacteria.sh` | *Escherichia coli* K-12 MG1655 | GCF_000005845.2 (ASM584v2) | 4,318 | 1.5 MB | `8b7b48b7b905313082141db4b77807b7aafebca7a5c40bbd42832baacd1cfd9b` |
| | *Salmonella enterica* Typhimurium LT2 | GCF_000006945.2 (ASM694v2) | 4,554 | 1.5 MB | `111a4660298f1d1edd99e5e7171538dd805f36606aab1ac9aa0bc0ac986ee0e5` |
| `get_chlamydia.sh` | *Chlamydia trachomatis* D/UW-3/CX | GCF_000008725.1 (ASM872v1) | 887 | 331 kB | `5676b7dbcfc374ce46719ff4898e50fd36c94eae294a14ce1b6c267c59baef98` |
| | *Chlamydia trachomatis* L2/434/Bu | GCF_000068585.1 (ASM6858v1) | 901 | 332 kB | `0a1b6862a3042da26364efb92f0f5af002e006b84bc60d31a013ed6ef4007a46` |
| | *Chlamydia muridarum* Nigg | GCF_000006685.1 (ASM668v1) | 898 | 347 kB | `9d66dd5e7b0d931d24db4386565effb7d906116b005ae4f82c75e0e97b6a6894` |
| | *Chlamydia pneumoniae* CWL029 | GCF_000008745.1 (ASM874v1) | 1,052 | 387 kB | `e08181a6794d0050ff3b2e123333a5978e864d06a61fea9930f1579eeea6a981` |

```
bash $T/data/get_mycoplasma.sh [DEST]        # default DEST: coal-tutorial-data
```

## Scripts

| Script | What it does |
|---|---|
| [`scripts/dnds_summary.py`](scripts/dnds_summary.py) | filters a `dnds.tsv` with the thresholds of [kaks.md](../kaks.md#filtering-before-interpretation), step by step; quartiles; highest and lowest dN/dS with gene names; `--by-species` for cluster runs; `--only-ids`; `--write-kept` |
| [`scripts/compare_methods.py`](scripts/compare_methods.py) | compares the methods of a `dnds.tsv` made with several `-m` |

Both read the tables only; neither changes a number.

## The numbers will differ when...

- **another DIAMOND version** finds slightly different hits or clusters (tutorials 2-4 used 2.2.8);
- **another Clustal Omega or MAFFT version** aligns differently (tutorial 3);
- **another COAL version** changes a documented policy (the CHANGELOG says when).

With the same versions and inputs, the tables are byte-identical between runs and whatever `-j` is (tutorial 2
shows the check).
