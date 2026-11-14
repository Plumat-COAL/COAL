#!/usr/bin/env bash
# Data of tutorial 01 (codon alignment): the coding sequences of three Mycoplasmoides genomes, which use
# NCBI genetic code 4 (TGA = Trp). About 0.7 MB in all.
#   bash docs/tutorials/data/get_mycoplasma.sh [DEST]      (default DEST: coal-tutorial-data)
DEST=${1:-coal-tutorial-data}
. "$(dirname "$0")/fetch_ncbi_cds.sh"
mkdir -p "$DEST"
fetch_cds GCF_000027325.1 ASM2732v1 mgenitalium_G37.cds.fna.gz     e7e499e12911e044fdb642a9056e64cf2a173e9837c867c492b99fa2abb7a24d
fetch_cds GCF_000027345.1 ASM2734v1 mpneumoniae_M129.cds.fna.gz    1ed0e736c2a3d163ef072aa5ffd46364718ef7e1427a231431ecd72b569ef6da
fetch_cds GCF_000092585.1 ASM9258v1 mgallisepticum_Rlow.cds.fna.gz 6e1e66665c4cf6d5169fb9a309400ec58cce6c4aab591843226f43abffdb04af
