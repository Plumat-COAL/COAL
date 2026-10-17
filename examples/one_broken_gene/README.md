# one_broken_gene: 19 good genes and one frameshift

Twenty synthetic orthologous genes (200 codons). Nineteen are fine. In `OG0001_sp14` one DNA base
has been deleted halfway along the gene (a frameshift, as caused by a sequencing or
assembly error), so from that point on its codons no longer spell the protein. A good tool
should flag that one gene clearly, explain what's wrong in plain words, and (if asked)
still produce the alignment for the other nineteen.
