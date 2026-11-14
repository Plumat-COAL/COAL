#!/usr/bin/env bash
# Data of tutorial 02 (pairs): the coding sequences of Escherichia coli K-12 MG1655 and Salmonella enterica
# serovar Typhimurium LT2. About 3 MB in all.
#   bash docs/tutorials/data/get_enterobacteria.sh [DEST]  (default DEST: coal-tutorial-data)
DEST=${1:-coal-tutorial-data}
. "$(dirname "$0")/fetch_ncbi_cds.sh"
mkdir -p "$DEST"
fetch_cds GCF_000005845.2 ASM584v2 ecoli_K12_MG1655.cds.fna.gz     8b7b48b7b905313082141db4b77807b7aafebca7a5c40bbd42832baacd1cfd9b
fetch_cds GCF_000006945.2 ASM694v2 styphimurium_LT2.cds.fna.gz     111a4660298f1d1edd99e5e7171538dd805f36606aab1ac9aa0bc0ac986ee0e5
