#!/usr/bin/env bash
# Data of tutorial 03 (clusters): the coding sequences of four Chlamydia genomes: two C. trachomatis strains
# (very close), C. muridarum (mouse; a sister species) and C. pneumoniae (distant). About 1.4 MB in all.
#   bash docs/tutorials/data/get_chlamydia.sh [DEST]       (default DEST: coal-tutorial-data)
DEST=${1:-coal-tutorial-data}
. "$(dirname "$0")/fetch_ncbi_cds.sh"
mkdir -p "$DEST"
fetch_cds GCF_000008725.1 ASM872v1  ctrachomatis_D.cds.fna.gz       5676b7dbcfc374ce46719ff4898e50fd36c94eae294a14ce1b6c267c59baef98
fetch_cds GCF_000068585.1 ASM6858v1 ctrachomatis_L2.cds.fna.gz      0a1b6862a3042da26364efb92f0f5af002e006b84bc60d31a013ed6ef4007a46
fetch_cds GCF_000006685.1 ASM668v1  cmuridarum_Nigg.cds.fna.gz      9d66dd5e7b0d931d24db4386565effb7d906116b005ae4f82c75e0e97b6a6894
fetch_cds GCF_000008745.1 ASM874v1  cpneumoniae_CWL029.cds.fna.gz   e08181a6794d0050ff3b2e123333a5978e864d06a61fea9930f1579eeea6a981
