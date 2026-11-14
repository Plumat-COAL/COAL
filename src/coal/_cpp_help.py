"""The --help texts of the commands coal-cpp runs (kaks, dnds, translate, codon-dnds), for an install
without coal-cpp. Generated from cpp/src/kaks_cli.cpp and cpp/src/dnds_cli.cpp by
`python tests/test_cli_commands.py --write`; do not edit by hand (tests/test_cli_commands.py compares
this copy with the sources and with the built binary)."""

HELP = {
    'kaks': (
        'usage: coal kaks -i IN -o OUT [-m METHOD ...] [-c TABLE] [-j N] [--pairs all|consecutive]\n'
        '                 [--digits N] [--yn-flavour kaks|paml] [--yn-weighting 0|1]\n'
        '                 [--comeron-flavour coal|gestimator] [--comeron-twosubs comeron|libsequence]\n'
        '                 [--comeron-min-sites N] [--li-flavour coal|seqinr]\n'
        '\n'
        "Ka/Ks for pairs of aligned coding sequences, written in KaKs_Calculator 2.0's tab-separated format.\n"
        "Takes KaKs_Calculator's -i -o -m -c for the seven KaKs_Calculator methods below (the default is YN;\n"
        "KaKs_Calculator's is MA), and has Comeron (1995) and Li (1993) as the R package orthologr computes them.\n"
        'The engine is derived from KaKs_Calculator 2.0, PAML yn00, libsequence and seqinr (GPL); see docs/kaks.md.\n'
        '\n'
        '  -i IN          AXT (name line, two aligned sequences, blank line) or aligned FASTA\n'
        "  -o OUT         output file ('-' or none: standard output)\n"
        "  -m METHOD      NG, LWL, MLWL, LPB, MLPB, YN (default), MYN; repeat for several. KaKs_Calculator's\n"
        '                 GY, MS, MA, NONE and gamma (G*) methods are not implemented natively and are refused.\n'
        "                 Comeron (alias C95): Comeron 1995 as libsequence's gestimator; Li: Li 1993 as\n"
        "                 seqinr::kaks (both: standard genetic code only; rows come after KaKs_Calculator's)\n"
        "  -c TABLE       NCBI genetic code (default 1), the tables COAL's --codontable accepts\n"
        '  -j, --jobs N   threads (default: all hardware threads)\n'
        '  --pairs P      FASTA input: all (every pair, default) or consecutive (records 1+2, 3+4, ...)\n'
        '  --digits N     significant digits of numbers (default 6, as KaKs_Calculator)\n'
        '  --yn-flavour F YN00 iteration as KaKs_Calculator (kaks, default) or PAML yn00 (paml)\n'
        '  --yn-weighting W  1: weight pathways (default); 0: equal weights (PAML yn00 weighting = 0)\n'
        "  --comeron-flavour F  Comeron: coal (default: COAL's stop-codon policy, NA when the correction is\n"
        "                 undefined) or gestimator (the numbers of orthologr 0.4.3's gestimator: stop codons\n"
        "                 compared, undefined logarithms set to 0, libsequence's two-position pathways)\n"
        '  --comeron-twosubs T  Comeron, codons that differ at two positions: comeron (the two pathways of\n'
        '                 Comeron 1995; default with --comeron-flavour coal) or libsequence (the branches\n'
        '                 libsequence counts, which depend on which sequence is first; default with gestimator)\n'
        '  --comeron-min-sites N  Comeron, coal flavour: a site class with fewer than N sites cannot make dN or\n'
        '                 dS NA (its undefined logarithm is 0, status warning; default 3; 0: any class can)\n'
        "  --li-flavour F Li: coal (default: COAL's stop-codon policy) or seqinr (stop codons compared, as\n"
        "                 seqinr::kaks); seqinr's sentinel 9.999999 is NA in both\n"
        '  -q, --quiet    no summary on stderr\n'
    ),
    'dnds': (
        'usage: coal dnds IN [-o OUT] [-m METHOD ...] [-c TABLE] [-j N] [--pairs all|consecutive]\n'
        '                 [--format auto|axt|fasta] [--digits N] [--yn-flavour kaks|paml] [--yn-weighting 0|1]\n'
        '                 [--comeron-flavour coal|gestimator] [--comeron-twosubs comeron|libsequence]\n'
        '                 [--comeron-min-sites N] [--li-flavour coal|seqinr] [--min-codons N] [--max-ds X]\n'
        '\n'
        'dN/dS for pairs of aligned coding sequences as a tidy TSV, one row per pair and method:\n'
        'seq1 seq2 method dN dS dNdS SE_dN SE_dS N_sites S_sites Nd Sd codons_compared status note\n'
        'status: ok, warning (internal stop codons skipped), not_converged (YN/MYN iteration: no estimate,\n'
        'dN and dS NA), saturated (dN or dS NA), no_data, error.\n'
        '\n'
        '  IN             aligned FASTA or AXT file (or -i IN)\n'
        '  -o OUT         output file (default: standard output)\n'
        '  -m METHOD      YN (default), NG, LWL, MLWL, LPB, MLPB, MYN, Comeron (alias C95), Li; repeat for\n'
        '                 several (Comeron and Li: standard genetic code only; SE for YN, MYN, NG and Li)\n'
        '  -c, --codontable TABLE  NCBI genetic code (default 1)\n'
        '  -j, --jobs N   threads (default: all hardware threads)\n'
        '  --pairs P      FASTA input: all (every pair, default) or consecutive (records 1+2, 3+4, ...)\n'
        "  --format F     input format (default auto: '>' first = FASTA, else AXT)\n"
        '  --digits N     significant digits of numbers (default 8)\n'
        '  --yn-flavour F YN00 iteration as KaKs_Calculator (kaks, default) or PAML yn00 (paml)\n'
        '  --yn-weighting W  1: weight pathways (default); 0: equal weights (PAML yn00 weighting = 0)\n'
        "  --comeron-flavour F  Comeron: coal (default: COAL's stop-codon policy, NA when the correction is\n"
        "                 undefined) or gestimator (the numbers of orthologr 0.4.3's gestimator: stop codons\n"
        "                 compared, undefined logarithms set to 0, libsequence's two-position pathways)\n"
        '  --comeron-twosubs T  Comeron, codons that differ at two positions: comeron (the two pathways of\n'
        '                 Comeron 1995; default with --comeron-flavour coal) or libsequence (the branches\n'
        '                 libsequence counts, which depend on which sequence is first; default with gestimator)\n'
        '  --comeron-min-sites N  Comeron, coal flavour: a site class with fewer than N sites cannot make dN or\n'
        '                 dS NA (its undefined logarithm is 0, status warning; default 3; 0: any class can)\n'
        "  --li-flavour F Li: coal (default: COAL's stop-codon policy) or seqinr (stop codons compared, as\n"
        "                 seqinr::kaks); seqinr's sentinel 9.999999 is NA in both\n"
        '  --min-codons N flag estimates from fewer than N compared codons: status warning and the note\n'
        '                 "fewer than N codons compared" (default 0: off). No row is dropped, no number changes\n'
        '  --max-ds X     flag estimates with dS above X (saturation): status warning and the note "dS above X"\n'
        '                 (default: off). See docs/kaks.md, "Filtering before interpretation"\n'
        '  -q, --quiet    no summary on stderr\n'
        '\n'
        'Exit status: 0 done; 1 the input cannot be read or the output cannot be written; 2 command-line error,\n'
        'or no pair has an estimate because of the input (every result error or no_data: a protein file,\n'
        'sequences that are not aligned).\n'
    ),
    'translate': (
        'usage: coal translate CDS.fasta [CDS.fasta ...] [-o PROTEINS.fasta] [--genes GENES.tsv] [-c TABLE] [-q]\n'
        '\n'
        'Translate coding sequences from their first base (a step of `coal pipeline`, docs/pipeline.md).\n'
        "A stop codon at the end is left out; an internal stop is written as '*'; a codon with an ambiguous base is X\n"
        'unless every reading gives the same residue; bases after the last complete codon are ignored. Genes that\n'
        'cannot be translated (empty, a character that is no nucleotide code, a repeated ID, an ID with a\n'
        'non-ASCII character) are left out of the protein file and listed as failed in GENES.tsv.\n'
        '\n'
        '  -o FILE        protein FASTA, one line per sequence (default: standard output)\n'
        '  --genes FILE   per-gene table: id nt_length aa_length terminal_stop internal_stops ambiguous_codons\n'
        '                 trailing_nt status code note   (status: ok, warning, failed)\n'
        '  -c, --codontable TABLE  NCBI genetic code (default 1)\n'
        '  -q, --quiet    no summary on stderr\n'
    ),
    'codon-dnds': (
        'usage: coal codon-dnds [-i GROUPS] -o PAIRS.tsv [-m METHOD ...] [-c TABLE] [-j N] [--group-column NAME]\n'
        '                 [--codon-out FILE] [--codon-dir DIR] [--kaks-out FILE] [--failures FILE] [--digits N]\n'
        '                 [--matrix NAME] [--gap-open X] [--gap-extend X] [--end-gaps free|penalized] [--max-cells N]\n'
        '                 [--tie-break leftmost|pwalign] [--flanks gap|align] [--yn-flavour kaks|paml]\n'
        '                 [--yn-weighting 0|1] [--comeron-flavour coal|gestimator] [--li-flavour coal|seqinr]\n'
        '                 [--comeron-twosubs comeron|libsequence] [--comeron-min-sites N]\n'
        '                 [--terminal-stops drop|compare] [--min-codons N] [--max-ds X] [--min-identity P]\n'
        '                 [--new-files] [-q]\n'
        '\n'
        'Protein alignment -> codon alignment -> Ka/Ks for groups of coding sequences (the engine of\n'
        '`coal pipeline`, docs/pipeline.md). GROUPS (a file, or standard input when -i is - or absent) is\n'
        'tab-separated, one group after the other:\n'
        '  G <group id> <all|rep|list> [nocodon] [swap]  start of a group; which pairs: all, member 1 against the\n'
        "                                  others, P lines; nocodon: do not write this group's codon (or protein)\n"
        '                                  alignment; swap: two members aligned here are aligned as member 2 against\n'
        '                                  member 1 (the rows come out in member order)\n'
        '  M <id> <CDS> [<protein row>]    a member; rows (when given, for every member) are the protein alignment\n'
        '  F <id> <reason>                 a member already known to be unusable (its pairs get status error)\n'
        '  H <qstart> <qend> <sstart> <send> <qseq_gapped> <sseq_gapped>   2 members: use this DIAMOND hit\n'
        '  P <i> <j>                       a pair of member numbers (from 0), with `list`\n'
        '  X <code> <reason>               the group is already known to have failed (its pairs get status error)\n'
        '  E                               end of the group\n'
        'Two members without rows are aligned here (global, BLOSUM62, needle defaults).\n'
        '\n'
        '  -o FILE          the pair table: [NAME] seq1 seq2 method dN dS dNdS SE_dN SE_dS N_sites S_sites Nd Sd\n'
        '                   codons_compared prot_aln_length prot_identity status note\n'
        '  --group-column NAME  write the group id as a first column called NAME\n'
        '  -m METHOD        YN (default), NG, LWL, MLWL, LPB, MLPB, MYN, Comeron (alias C95), Li; repeat for\n'
        '                   several (Comeron and Li: genetic code 1 only)\n'
        '  --tie-break T    which of several protein alignments with the best score: leftmost (gaps as far left\n'
        '                   as possible; the default) or pwalign (the one pwalign::pairwiseAlignment returns)\n'
        '  --comeron-flavour F, --li-flavour F  as in `coal dnds` (coal, the default; gestimator; seqinr)\n'
        '  --comeron-twosubs T, --comeron-min-sites N  as in `coal dnds`\n'
        '  --terminal-stops S  with the flavours gestimator / seqinr: drop (default) ends a pair before the\n'
        '                   terminal stop codons; compare takes the two stop codons as the last compared column\n'
        '                   when both genes end in one (as orthologr 0.4.3 does)\n'
        '  -c, --codontable TABLE  NCBI genetic code (default 1)\n'
        '  -j, --jobs N     threads (default: all hardware threads)\n'
        '  --codon-out FILE the codon alignments, one after the other (pairs: records 1+2, 3+4, ...)\n'
        '  --codon-dir DIR  one codon alignment file per group in DIR\n'
        '  --protein-out FILE, --protein-dir DIR  the protein alignments the codon alignments were made from, in\n'
        '                   the same layout: the same groups, records and order (a codon row is 3 x its protein row)\n'
        "  --kaks-out FILE  the same results in KaKs_Calculator's format (pairs named seq1&seq2)\n"
        '  --failures FILE  genes and groups that could not be used: group seq code message\n'
        '  --new-files      create every output file anew: a file or symbolic link of that name is removed first,\n'
        '                   never written through (`coal pipeline` writes NAME.tmp this way and renames it)\n'
        '  --digits N       significant digits (default 8)\n'
        '  --min-codons N, --max-ds X  as in `coal dnds`: flag estimates from fewer than N compared codons or with\n'
        '                   dS above X (status warning and a note; no row dropped, no number changed; default off)\n'
        '  --min-identity P the same kind of flag for a pair whose protein identity (the column prot_identity,\n'
        "                   percent of the alignment's columns) is below P: note `protein identity below P%`\n"
        '  -q, --quiet      no summary on stderr\n'
    ),
}
