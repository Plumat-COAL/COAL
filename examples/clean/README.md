# clean: nuclear genes, standard genetic code

Eight related protein-coding genes (synthetic, 140-144 codons each, aligned in 150 columns) with a
ready-made protein alignment and their unaligned coding DNA, each ending in a stop codon. Everything fits the
standard genetic code (NCBI code 1), the IDs match exactly, and there are no frameshifts,
so COAL should report every sequence as OK. Use it as a first test:
`python -m coal examples/clean/protein.aln.fasta examples/clean/dna.fasta`.
