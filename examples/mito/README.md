# mito: vertebrate mitochondrial genes (genetic code 2)

Twenty-four synthetic mitochondrial-style genes (COX1, COX2, COX3 and CYTB from six
vertebrates, about 110 codons each). Mitochondria read DNA slightly differently from the
nucleus: TGA means tryptophan (W) instead of "stop", and ATA means methionine (M). Every gene
here uses TGA for at least one W, so with the default standard code (1) COAL finds codons
that don't match their amino acid. Re-running with the vertebrate mitochondrial code fixes
it: `python -m coal examples/mito/protein.aln.fasta examples/mito/dna.fasta -codontable 2`.
