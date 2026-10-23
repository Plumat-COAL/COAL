# Stress test 2026-10-23

COAL v0.4.2 (commit `c6f5728`, branch `main` of 2026-10-23), both engines, against the simulated truth, PAL2NAL v14
and each other, with a bias towards very long genes and very large files. Everything here was run on one machine in
one day; the tools that generated the inputs and ran the cases are in `tests/tools/` (`stress_gen.py`,
`stress_run.py`, `stress_resume.py`, `stress_mcp.py`, `stress_windows.py`, `stress_winjob.py`, `stress_table.py`),
so every row of the table can be rerun from its command. No generated input is committed (the largest were 1.6 GB).

## Machine

| | |
|---|---|
| CPU | AMD Ryzen 7 5800X3D, 8 cores / 16 threads |
| RAM | 32 GB (about 11 GB free during the runs; Windows and the browser held the rest) |
| Disk | one NVMe shared by Windows and WSL |
| Linux | Ubuntu 26.04.1 under **WSL 1** (kernel 4.4.0-26100-Microsoft) on Windows 11 Home 10.0.26200 |
| Python | 3.14.4 (WSL), 3.13.5 (Windows) |
| Compiler | g++ 15.2.0 (`build/coal-cpp`, Release); x86_64-w64-mingw32-g++-posix 15 for a static `coal-cpp.exe` |
| PAL2NAL | v14 from the official site (`scripts/fetch_pal2nal.py`), perl 5.40.1 |
| App | mcp 2.2.0 (WSL), gradio 5.49.1 (Windows only, see finding 8) |

WSL 1 is not a Linux kernel, and that cost a few cases (listed at the end): Python 3.14's default `forkserver`
start method cannot start there, so the Python engine's worker pools were run with `fork` through a one-line `.pth`
file in the virtualenv (the baseline suite shows 84 errors from this alone; all pass under `fork`); `ulimit -v`,
`ulimit -f`, `/dev/full` and tmpfs size limits are not enforced, so the full-disk cases could not be produced and the
memory-limit case was done on Windows with a Job Object instead.

**Measurements.** Wall time is the whole process, start-up included. Peak memory is the proportional set size
(PSS) summed over the process tree, sampled every 50 ms with psutil, so pages shared by forked workers count once
(the summed RSS, which counts them once per worker, is in the JSON logs as `peak_rss_mb` and is 3-4x larger for
17-process runs). "Correct" means byte-identical to the simulated truth (plain `-output fasta` runs; the truth is
built from the same codons), byte-identical to PAL2NAL v14's stdout where PAL2NAL could run, and identical between
the two engines in stdout (or the `-o` file), stderr, exit status and `--report` JSON without its engine fields.

## Every case

| case | input shape | input MB | engine | time s | peak MB | exit | correct? |
|---|---|---:|---|---:|---:|---:|---|
| long_fasta_10000 | 3 x 10000 | 0.1 | python | 0.09 | 13 | 0 | truth =, PAL2NAL =, engines = |
| long_fasta_10000 | 3 x 10000 | 0.1 | cpp | 0.08 | 14 | 0 | truth =, PAL2NAL =, engines = |
| long_fasta_10000 | 3 x 10000 | 0.1 | PAL2NAL v14 | 0.88 | 11 | 0 | output |
| long_gaps_10000 | 3 x 10000, gaps | 0.1 | python | 0.09 | 14 | 0 | truth =, PAL2NAL =, engines = |
| long_gaps_10000 | 3 x 10000, gaps | 0.1 | cpp | 0.10 | 13 | 0 | truth =, PAL2NAL =, engines = |
| long_gaps_10000 | 3 x 10000, gaps | 0.1 | PAL2NAL v14 | 0.80 | 11 | 0 | output |
| long_clustal_in_10000 | 3 x 10000, gaps, CLUSTAL in | 0.1 | python | 0.10 | 12 | 0 | truth =, PAL2NAL =, engines = |
| long_clustal_in_10000 | 3 x 10000, gaps, CLUSTAL in | 0.1 | cpp | 0.09 | 13 | 0 | truth =, PAL2NAL =, engines = |
| long_clustal_in_10000 | 3 x 10000, gaps, CLUSTAL in | 0.1 | PAL2NAL v14 | 0.83 | 11 | 0 | output |
| long_clustal_out_10000 | 3 x 10000, gaps, CLUSTAL out | 0.1 | python | 0.10 | 13 | 0 | PAL2NAL =, engines = |
| long_clustal_out_10000 | 3 x 10000, gaps, CLUSTAL out | 0.1 | cpp | 0.10 | 12 | 0 | PAL2NAL =, engines = |
| long_clustal_out_10000 | 3 x 10000, gaps, CLUSTAL out | 0.1 | PAL2NAL v14 | 0.81 | 11 | 0 | output |
| long_codon_out_10000 | 3 x 10000, gaps, codon out | 0.1 | python | 0.10 | 14 | 0 | PAL2NAL =, engines = |
| long_codon_out_10000 | 3 x 10000, gaps, codon out | 0.1 | cpp | 0.10 | 14 | 0 | PAL2NAL =, engines = |
| long_codon_out_10000 | 3 x 10000, gaps, codon out | 0.1 | PAL2NAL v14 | 0.83 | 10 | 0 | output |
| long_report_10000 | 3 x 10000, gaps, --report | 0.1 | python | 0.10 | 14 | 0 | truth =, engines = |
| long_report_10000 | 3 x 10000, gaps, --report | 0.1 | cpp | 0.10 | 13 | 0 | truth =, engines = |
| long_fasta_50000 | 3 x 50000 | 0.6 | python | 0.14 | 19 | 0 | truth =, PAL2NAL failed, engines = |
| long_fasta_50000 | 3 x 50000 | 0.6 | cpp | 0.10 | 14 | 0 | truth =, PAL2NAL failed, engines = |
| long_fasta_50000 | 3 x 50000 | 0.6 | PAL2NAL v14 | 9.51 | 16 | 255 | no output (fails above about 22,000 residues: documented) |
| long_gaps_50000 | 3 x 50000, gaps | 0.6 | python | 0.14 | 17 | 0 | truth =, PAL2NAL failed, engines = |
| long_gaps_50000 | 3 x 50000, gaps | 0.6 | cpp | 0.10 | 14 | 0 | truth =, PAL2NAL failed, engines = |
| long_gaps_50000 | 3 x 50000, gaps | 0.6 | PAL2NAL v14 | 8.54 | 15 | 255 | no output (fails above about 22,000 residues: documented) |
| long_clustal_in_50000 | 3 x 50000, gaps, CLUSTAL in | 0.6 | python | 0.15 | 17 | 0 | truth =, PAL2NAL failed, engines = |
| long_clustal_in_50000 | 3 x 50000, gaps, CLUSTAL in | 0.6 | cpp | 0.10 | 14 | 0 | truth =, PAL2NAL failed, engines = |
| long_clustal_in_50000 | 3 x 50000, gaps, CLUSTAL in | 0.6 | PAL2NAL v14 | 8.71 | 16 | 255 | no output (fails above about 22,000 residues: documented) |
| long_clustal_out_50000 | 3 x 50000, gaps, CLUSTAL out | 0.6 | python | 0.15 | 18 | 0 | PAL2NAL failed, engines = |
| long_clustal_out_50000 | 3 x 50000, gaps, CLUSTAL out | 0.6 | cpp | 0.14 | 15 | 0 | PAL2NAL failed, engines = |
| long_clustal_out_50000 | 3 x 50000, gaps, CLUSTAL out | 0.6 | PAL2NAL v14 | 8.56 | 15 | 255 | no output (fails above about 22,000 residues: documented) |
| long_codon_out_50000 | 3 x 50000, gaps, codon out | 0.6 | python | 0.19 | 19 | 0 | PAL2NAL failed, engines = |
| long_codon_out_50000 | 3 x 50000, gaps, codon out | 0.6 | cpp | 0.14 | 15 | 0 | PAL2NAL failed, engines = |
| long_codon_out_50000 | 3 x 50000, gaps, codon out | 0.6 | PAL2NAL v14 | 8.58 | 15 | 255 | no output (fails above about 22,000 residues: documented) |
| long_report_50000 | 3 x 50000, gaps, --report | 0.6 | python | 0.15 | 19 | 0 | truth =, engines = |
| long_report_50000 | 3 x 50000, gaps, --report | 0.6 | cpp | 0.10 | 13 | 0 | truth =, engines = |
| long_fasta_200000 | 3 x 200000 | 2.4 | python | 0.34 | 36 | 0 | truth =, engines = |
| long_fasta_200000 | 3 x 200000 | 2.4 | cpp | 0.14 | 18 | 0 | truth =, engines = |
| long_gaps_200000 | 3 x 200000, gaps | 2.4 | python | 0.32 | 32 | 0 | truth =, engines = |
| long_gaps_200000 | 3 x 200000, gaps | 2.4 | cpp | 0.14 | 18 | 0 | truth =, engines = |
| long_clustal_in_200000 | 3 x 200000, gaps, CLUSTAL in | 2.4 | python | 0.32 | 32 | 0 | truth =, engines = |
| long_clustal_in_200000 | 3 x 200000, gaps, CLUSTAL in | 2.4 | cpp | 0.14 | 14 | 0 | truth =, engines = |
| long_clustal_out_200000 | 3 x 200000, gaps, CLUSTAL out | 2.4 | python | 0.37 | 33 | 0 | engines = |
| long_clustal_out_200000 | 3 x 200000, gaps, CLUSTAL out | 2.4 | cpp | 0.26 | 40 | 0 | engines = |
| long_codon_out_200000 | 3 x 200000, gaps, codon out | 2.4 | python | 0.46 | 31 | 0 | engines = |
| long_codon_out_200000 | 3 x 200000, gaps, codon out | 2.4 | cpp | 0.27 | 40 | 0 | engines = |
| long_report_200000 | 3 x 200000, gaps, --report | 2.4 | python | 0.34 | 32 | 0 | truth =, engines = |
| long_report_200000 | 3 x 200000, gaps, --report | 2.4 | cpp | 0.15 | 17 | 0 | truth =, engines = |
| long_fasta_1000000 | 3 x 1000000 | 12.2 | python | 1.21 | 232 | 0 | truth =, engines = |
| long_fasta_1000000 | 3 x 1000000 | 12.2 | cpp | 0.34 | 106 | 0 | truth =, engines = |
| long_gaps_1000000 | 3 x 1000000, gaps | 11.7 | python | 1.17 | 226 | 0 | truth =, engines = |
| long_gaps_1000000 | 3 x 1000000, gaps | 11.7 | cpp | 0.32 | 102 | 0 | truth =, engines = |
| long_clustal_in_1000000 | 3 x 1000000, gaps, CLUSTAL in | 12.1 | python | 1.21 | 196 | 0 | truth =, engines = |
| long_clustal_in_1000000 | 3 x 1000000, gaps, CLUSTAL in | 12.1 | cpp | 0.36 | 108 | 0 | truth =, engines = |
| long_clustal_out_1000000 | 3 x 1000000, gaps, CLUSTAL out | 11.7 | python | 1.51 | 214 | 0 | engines = |
| long_clustal_out_1000000 | 3 x 1000000, gaps, CLUSTAL out | 11.7 | cpp | 1.00 | 102 | 0 | engines = |
| long_codon_out_1000000 | 3 x 1000000, gaps, codon out | 11.7 | python | 1.76 | 198 | 0 | engines = |
| long_codon_out_1000000 | 3 x 1000000, gaps, codon out | 11.7 | cpp | 1.19 | 102 | 0 | engines = |
| long_report_1000000 | 3 x 1000000, gaps, --report | 11.7 | python | 1.24 | 232 | 0 | truth =, engines = |
| long_report_1000000 | 3 x 1000000, gaps, --report | 11.7 | cpp | 0.41 | 103 | 0 | truth =, engines = |
| long1M_frameshift | 1 x 1,000,000, frameshift at 999,900, UTRs | 4.1 | python | 2.15 | 220 | 0 | truth =, engines = |
| long1M_frameshift | 1 x 1,000,000, frameshift at 999,900, UTRs | 4.1 | cpp | 0.21 | 61 | 0 | truth =, engines = |
| long1M_utr | 1 x 1,000,000, 5' UTR 100 nt | 4.1 | python | 1.58 | 222 | 0 | truth =, engines = |
| long1M_utr | 1 x 1,000,000, 5' UTR 100 nt | 4.1 | cpp | 0.20 | 23 | 0 | truth =, engines = |
| long1M_late_mismatch | 1 x 1,000,000, last residue wrong | 4.1 | python | 30.20 | 694 | 0 | truth =, engines = |
| long1M_late_mismatch | 1 x 1,000,000, last residue wrong | 4.1 | cpp | 5.73 | 175 | 0 | truth =, engines = |
| long1M_iupac | 1 x 1,000,000, 1% NNN-type codons | 4.1 | python | 1.82 | 222 | 0 | truth =, engines = |
| long1M_iupac | 1 x 1,000,000, 1% NNN-type codons | 4.1 | cpp | 0.22 | 60 | 0 | truth =, engines = |
| repeat_polyA_late_10000 | 1 x 10000 poly-A, last residue wrong, DNA repeat 2x | 0.1 | python | 0.10 | 13 | 0 | truth DIFF, engines =; first run: the generator's truth was wrong here (COAL's leftmost exact fit is right); superseded by repeat_polyA_prefix |
| repeat_polyA_late_10000 | 1 x 10000 poly-A, last residue wrong, DNA repeat 2x | 0.1 | cpp | 0.42 | 13 | 0 | truth DIFF, engines =; first run: the generator's truth was wrong here (COAL's leftmost exact fit is right); superseded by repeat_polyA_prefix |
| repeat_polyA_late_50000 | 1 x 50000 poly-A, last residue wrong, DNA repeat 2x | 0.4 | python | 0.26 | 24 | 0 | truth DIFF, engines =; first run: the generator's truth was wrong here (COAL's leftmost exact fit is right); superseded by repeat_polyA_prefix |
| repeat_polyA_late_50000 | 1 x 50000 poly-A, last residue wrong, DNA repeat 2x | 0.4 | cpp | 9.00 | 14 | 0 | truth DIFF, engines =; first run: the generator's truth was wrong here (COAL's leftmost exact fit is right); superseded by repeat_polyA_prefix |
| repeat_polyA_late_200000 | 1 x 200000 poly-A, last residue wrong, DNA repeat 2x | 1.4 | python | 0.78 | 70 | 0 | truth DIFF, engines =; first run: the generator's truth was wrong here (COAL's leftmost exact fit is right); superseded by repeat_polyA_prefix |
| repeat_polyA_late_200000 | 1 x 200000 poly-A, last residue wrong, DNA repeat 2x | 1.4 | cpp | 137.89 | 18 | 0 | truth DIFF, engines =; first run: the generator's truth was wrong here (COAL's leftmost exact fit is right); superseded by repeat_polyA_prefix |
| repeat_GPP_late_50000 | 1 x 50,000 (GPP)n, last residue wrong, DNA repeat 2x | 0.4 | python | 0.23 | 30 | 0 | truth DIFF, engines =; first run: the generator's truth was wrong here (COAL's leftmost exact fit is right); superseded by repeat_GPP_prefix |
| repeat_GPP_late_50000 | 1 x 50,000 (GPP)n, last residue wrong, DNA repeat 2x | 0.4 | cpp | 2.20 | 13 | 0 | truth DIFF, engines =; first run: the generator's truth was wrong here (COAL's leftmost exact fit is right); superseded by repeat_GPP_prefix |
| repeat_polyA_prefix_10000 | 1 x 10000 poly-A gene after 10000 codons of poly-A | 0.1 | python | 0.11 | 17 | 0 | truth =, PAL2NAL =, engines =; FINDING 10: quadratic scan in coal-cpp |
| repeat_polyA_prefix_10000 | 1 x 10000 poly-A gene after 10000 codons of poly-A | 0.1 | cpp | 0.42 | 13 | 0 | truth =, PAL2NAL =, engines =; FINDING 10: quadratic scan in coal-cpp |
| repeat_polyA_prefix_10000 | 1 x 10000 poly-A gene after 10000 codons of poly-A | 0.1 | PAL2NAL v14 | 1.43 | 6 | 0 | output |
| repeat_polyA_prefix_50000 | 1 x 50000 poly-A gene after 50000 codons of poly-A | 0.4 | python | 0.24 | 28 | 0 | truth =, engines =; FINDING 10: quadratic scan in coal-cpp |
| repeat_polyA_prefix_50000 | 1 x 50000 poly-A gene after 50000 codons of poly-A | 0.4 | cpp | 8.51 | 14 | 0 | truth =, engines =; FINDING 10: quadratic scan in coal-cpp |
| repeat_polyA_prefix_200000 | 1 x 200000 poly-A gene after 200000 codons of poly-A | 1.4 | python | 0.73 | 71 | 0 | truth =, engines =; FINDING 10: quadratic scan in coal-cpp |
| repeat_polyA_prefix_200000 | 1 x 200000 poly-A gene after 200000 codons of poly-A | 1.4 | cpp | 138.03 | 18 | 0 | truth =, engines =; FINDING 10: quadratic scan in coal-cpp |
| repeat_GPP_prefix_50000 | 1 x 50,000 (GPP)n gene after 50,000 codons of (GPP)n | 0.4 | python | 0.22 | 27 | 0 | truth =, engines =; FINDING 10: quadratic scan in coal-cpp |
| repeat_GPP_prefix_50000 | 1 x 50,000 (GPP)n gene after 50,000 codons of (GPP)n | 0.4 | cpp | 3.34 | 13 | 0 | truth =, engines =; FINDING 10: quadratic scan in coal-cpp |
| repeat_polyA_late_50000 | 1 x 50000 poly-A, last residue wrong | 0.2 | python | 3.12 | 28 | 0 | truth =, engines = |
| repeat_polyA_late_50000 | 1 x 50000 poly-A, last residue wrong | 0.2 | cpp | 0.27 | 13 | 0 | truth =, engines = |
| repeat_polyA_late_200000 | 1 x 200000 poly-A, last residue wrong | 0.8 | python | 14.92 | 66 | 0 | truth =, engines = |
| repeat_polyA_late_200000 | 1 x 200000 poly-A, last residue wrong | 0.8 | cpp | 0.99 | 31 | 0 | truth =, engines = |
| late_mismatch_50000 | 1 x 50000, last residue wrong | 0.2 | python | 1.61 | 48 | 0 | truth =, engines = |
| late_mismatch_50000 | 1 x 50000, last residue wrong | 0.2 | cpp | 0.35 | 15 | 0 | truth =, engines = |
| late_mismatch_200000 | 1 x 200000, last residue wrong | 0.8 | python | 6.72 | 157 | 0 | truth =, engines = |
| late_mismatch_200000 | 1 x 200000, last residue wrong | 0.8 | cpp | 1.31 | 41 | 0 | truth =, engines = |
| late_mismatch_200000_q | 1 x 200,000, last residue wrong, -q | 0.8 | python | 0.45 | 54 | 0 | truth =, engines = |
| late_mismatch_200000_q | 1 x 200,000, last residue wrong, -q | 0.8 | cpp | 0.17 | 23 | 0 | truth =, engines = |
| long1M_late_mismatch_q | 1 x 1,000,000, last residue wrong, -q | 4.1 | python | 1.67 | 266 | 0 | truth =, engines = |
| long1M_late_mismatch_q | 1 x 1,000,000, last residue wrong, -q | 4.1 | cpp | 0.51 | 91 | 0 | truth =, engines = |
| long1M_late_mismatch_nomismatch | 1 x 1,000,000, last residue wrong, -nomismatch | 4.1 | python | 34.47 | 705 | 0 | engines = |
| long1M_late_mismatch_nomismatch | 1 x 1,000,000, last residue wrong, -nomismatch | 4.1 | cpp | 5.59 | 192 | 0 | engines = |
| repeat_polyA_prefix_1000000 | 1 x 1,000,000 poly-A gene after 1,000,000 codons of poly-A | 7.1 | python | 3.99 | 295 | 0 | truth =, engines DIFF; FINDING 10: quadratic scan, killed at 30 min |
| repeat_polyA_prefix_1000000 | 1 x 1,000,000 poly-A gene after 1,000,000 codons of poly-A | 7.1 | cpp | 1800.05 | 67 | timeout | did not finish in 30 min (killed); FINDING 10: quadratic scan, killed at 30 min |
| fmt_crlf | 200 x 300, CRLF | 0.2 | python | 0.10 | 12 | 0 | truth =, PAL2NAL =, engines = |
| fmt_crlf | 200 x 300, CRLF | 0.2 | cpp | 0.10 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_crlf | 200 x 300, CRLF | 0.2 | PAL2NAL v14 | 0.26 | 7 | 0 | output |
| fmt_cr_only | 200 x 300, CR line ends | 0.2 | python | 0.11 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_cr_only | 200 x 300, CR line ends | 0.2 | cpp | 0.10 | 12 | 0 | truth =, PAL2NAL =, engines = |
| fmt_cr_only | 200 x 300, CR line ends | 0.2 | PAL2NAL v14 | 0.24 | 6 | 0 | output |
| fmt_bom | 200 x 300, UTF-8 BOM | 0.2 | python | 0.10 | 13 | 0 | truth =, PAL2NAL failed, engines = |
| fmt_bom | 200 x 300, UTF-8 BOM | 0.2 | cpp | 0.09 | 13 | 0 | truth =, PAL2NAL failed, engines = |
| fmt_bom | 200 x 300, UTF-8 BOM | 0.2 | PAL2NAL v14 | 0.04 | 0 | 0 | no output |
| fmt_blank_lines | 200 x 300, blank lines | 0.2 | python | 0.10 | 12 | 0 | truth =, PAL2NAL =, engines = |
| fmt_blank_lines | 200 x 300, blank lines | 0.2 | cpp | 0.09 | 14 | 0 | truth =, PAL2NAL =, engines = |
| fmt_blank_lines | 200 x 300, blank lines | 0.2 | PAL2NAL v14 | 0.26 | 6 | 0 | output |
| fmt_empty_header | 200 x 300, '>' without a name | 0.2 | python | 0.07 | 14 | 1 | truth DIFF, PAL2NAL failed, engines =; error by design: '>' without a name |
| fmt_empty_header | 200 x 300, '>' without a name | 0.2 | cpp | 0.08 | 13 | 1 | truth DIFF, PAL2NAL failed, engines =; error by design: '>' without a name |
| fmt_empty_header | 200 x 300, '>' without a name | 0.2 | PAL2NAL v14 | 0.02 | 0 | 0 | no output |
| fmt_ids_spaces | 200 x 300, descriptions with spaces | 0.3 | python | 0.10 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_ids_spaces | 200 x 300, descriptions with spaces | 0.3 | cpp | 0.10 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_ids_spaces | 200 x 300, descriptions with spaces | 0.3 | PAL2NAL v14 | 0.22 | 6 | 0 | output |
| fmt_ids_tabs | 200 x 300, tabs in headers | 0.2 | python | 0.10 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_ids_tabs | 200 x 300, tabs in headers | 0.2 | cpp | 0.10 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_ids_tabs | 200 x 300, tabs in headers | 0.2 | PAL2NAL v14 | 0.26 | 6 | 0 | output |
| fmt_ids_unicode_desc | 200 x 300, non-ASCII description | 0.2 | python | 0.10 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_ids_unicode_desc | 200 x 300, non-ASCII description | 0.2 | cpp | 0.09 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_ids_unicode_desc | 200 x 300, non-ASCII description | 0.2 | PAL2NAL v14 | 0.24 | 6 | 0 | output |
| fmt_ids_nonascii | 200 x 300, non-ASCII IDs | 0.2 | python | 0.07 | 14 | 1 | truth DIFF, PAL2NAL DIFF, engines =; error by design: non-ASCII ID |
| fmt_ids_nonascii | 200 x 300, non-ASCII IDs | 0.2 | cpp | 0.09 | 14 | 1 | truth DIFF, PAL2NAL DIFF, engines =; error by design: non-ASCII ID |
| fmt_ids_nonascii | 200 x 300, non-ASCII IDs | 0.2 | PAL2NAL v14 | 0.26 | 6 | 0 | output |
| fmt_ids_10000_chars | 200 x 300, 10,000-char IDs | 4.2 | python | 0.15 | 19 | 0 | truth =, PAL2NAL =, engines = |
| fmt_ids_10000_chars | 200 x 300, 10,000-char IDs | 4.2 | cpp | 0.14 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_ids_10000_chars | 200 x 300, 10,000-char IDs | 4.2 | PAL2NAL v14 | 0.31 | 23 | 0 | output |
| fmt_ids_10000_chars_clustal_out | 200 x 300, 10,000-char IDs, clustal out | 4.2 | python | 0.26 | 98 | 0 | PAL2NAL =, engines = |
| fmt_ids_10000_chars_clustal_out | 200 x 300, 10,000-char IDs, clustal out | 4.2 | cpp | 0.14 | 13 | 0 | PAL2NAL =, engines = |
| fmt_ids_10000_chars_clustal_out | 200 x 300, 10,000-char IDs, clustal out | 4.2 | PAL2NAL v14 | 0.36 | 23 | 0 | output |
| fmt_ids_10000_chars_clustal_in | 20 x 300, 10,000-char IDs, CLUSTAL in | 1.2 | python | 0.10 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_ids_10000_chars_clustal_in | 20 x 300, 10,000-char IDs, CLUSTAL in | 1.2 | cpp | 0.10 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_ids_10000_chars_clustal_in | 20 x 300, 10,000-char IDs, CLUSTAL in | 1.2 | PAL2NAL v14 | 0.05 | 0 | 0 | output |
| fmt_dot_gap | 200 x 300, '.' gaps | 0.2 | python | 0.10 | 13 | 0 | truth =, PAL2NAL DIFF, engines =; documented difference: '.' is a gap for COAL, PAL2NAL misaligns the row |
| fmt_dot_gap | 200 x 300, '.' gaps | 0.2 | cpp | 0.09 | 13 | 0 | truth =, PAL2NAL DIFF, engines =; documented difference: '.' is a gap for COAL, PAL2NAL misaligns the row |
| fmt_dot_gap | 200 x 300, '.' gaps | 0.2 | PAL2NAL v14 | 0.24 | 6 | 0 | output |
| fmt_lower | 200 x 300, lower-case DNA | 0.2 | python | 0.10 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_lower | 200 x 300, lower-case DNA | 0.2 | cpp | 0.09 | 13 | 0 | truth =, PAL2NAL =, engines = |
| fmt_lower | 200 x 300, lower-case DNA | 0.2 | PAL2NAL v14 | 0.24 | 6 | 0 | output |
| fmt_lower_nogap | 200 x 300, lower-case DNA, -nogap | 0.2 | python | 0.12 | 16 | 0 | PAL2NAL =, engines = |
| fmt_lower_nogap | 200 x 300, lower-case DNA, -nogap | 0.2 | cpp | 0.10 | 12 | 0 | PAL2NAL =, engines = |
| fmt_lower_nogap | 200 x 300, lower-case DNA, -nogap | 0.2 | PAL2NAL v14 | 0.26 | 6 | 0 | output |
| fmt_rna | 200 x 300, RNA | 0.2 | python | 0.10 | 14 | 0 | truth =, PAL2NAL =, engines = |
| fmt_rna | 200 x 300, RNA | 0.2 | cpp | 0.10 | 12 | 0 | truth =, PAL2NAL =, engines = |
| fmt_rna | 200 x 300, RNA | 0.2 | PAL2NAL v14 | 0.21 | 6 | 0 | output |
| fmt_iupac_all | 200 x 300, N in every 4-fold codon | 0.2 | python | 0.12 | 16 | 0 | truth =, PAL2NAL =, engines = |
| fmt_iupac_all | 200 x 300, N in every 4-fold codon | 0.2 | cpp | 0.10 | 12 | 0 | truth =, PAL2NAL =, engines = |
| fmt_iupac_all | 200 x 300, N in every 4-fold codon | 0.2 | PAL2NAL v14 | 0.24 | 6 | 0 | output |
| fmt_all_X | 200 x 300, protein all X | 0.2 | python | 0.11 | 16 | 0 | truth =, PAL2NAL =, engines = |
| fmt_all_X | 200 x 300, protein all X | 0.2 | cpp | 0.09 | 14 | 0 | truth =, PAL2NAL =, engines = |
| fmt_all_X | 200 x 300, protein all X | 0.2 | PAL2NAL v14 | 0.12 | 6 | 0 | output |
| fmt_empty_dna_record | 200 x 300, DNA record 7 empty | 0.2 | python | 0.10 | 13 | 1 | truth DIFF, PAL2NAL failed, engines =; error by design: a gene without DNA |
| fmt_empty_dna_record | 200 x 300, DNA record 7 empty | 0.2 | cpp | 0.09 | 14 | 1 | truth DIFF, PAL2NAL failed, engines =; error by design: a gene without DNA |
| fmt_empty_dna_record | 200 x 300, DNA record 7 empty | 0.2 | PAL2NAL v14 | 0.18 | 7 | 0 | no output |
| fmt_empty_dna_record_keepgoing | 200 x 300, DNA record 7 empty, --keep-going | 0.2 | python | 0.12 | 16 | 2 | truth DIFF, engines =; by design: the empty gene is left out (exit 2) |
| fmt_empty_dna_record_keepgoing | 200 x 300, DNA record 7 empty, --keep-going | 0.2 | cpp | 0.09 | 13 | 2 | truth DIFF, engines =; by design: the empty gene is left out (exit 2) |
| clustal_200k_cols | 5 x 200,000, CLUSTAL in with # marks | 4.2 | python | 0.52 | 36 | 0 | truth =, engines = |
| clustal_200k_cols | 5 x 200,000, CLUSTAL in with # marks | 4.2 | cpp | 0.19 | 37 | 0 | truth =, engines = |
| clustal_200k_cols_blockonly | 5 x 200,000, CLUSTAL in, -blockonly | 4.2 | python | 0.63 | 38 | 0 | engines = |
| clustal_200k_cols_blockonly | 5 x 200,000, CLUSTAL in, -blockonly | 4.2 | cpp | 0.21 | 41 | 0 | engines = |
| clustal_200k_cols_clustal_out | 5 x 200,000, CLUSTAL in, clustal out (marks shown) | 4.2 | python | 0.63 | 35 | 0 | engines = |
| clustal_200k_cols_clustal_out | 5 x 200,000, CLUSTAL in, clustal out (marks shown) | 4.2 | cpp | 0.36 | 57 | 0 | engines = |
| clustal_200k_cols_one_block | 5 x 200,000, CLUSTAL in, one 200,000-char block | 4.1 | python | 0.61 | 37 | 0 | engines = |
| clustal_200k_cols_one_block | 5 x 200,000, CLUSTAL in, one 200,000-char block | 4.1 | cpp | 0.19 | 37 | 0 | engines = |
| line_100MB_utr | 1 x 400 + 100 MB 3' UTR on one line | 100.0 | python | 2.92 | 543 | 0 | truth =, engines = |
| line_100MB_utr | 1 x 400 + 100 MB 3' UTR on one line | 100.0 | cpp | 1.02 | 306 | 0 | truth =, engines = |
| fmt_mixed_lf_then_cr | 2,000 x 400, LF for 300 records then CR | 3.3 | python | 0.11 | 16 | 1 | truth DIFF, PAL2NAL DIFF, engines =; FINDING 4: PAL2NAL reads this file |
| fmt_mixed_lf_then_cr | 2,000 x 400, LF for 300 records then CR | 3.3 | cpp | 0.10 | 14 | 1 | truth DIFF, PAL2NAL DIFF, engines =; FINDING 4: PAL2NAL reads this file |
| fmt_mixed_lf_then_cr | 2,000 x 400, LF for 300 records then CR | 3.3 | PAL2NAL v14 | 3.01 | 37 | 0 | output |
| fmt_mixed_lf_then_crlf | 2,000 x 400, LF for 300 records then CRLF | 3.3 | python | 0.26 | 21 | 0 | truth =, PAL2NAL =, engines =; FINDING 4: PAL2NAL reads this file |
| fmt_mixed_lf_then_crlf | 2,000 x 400, LF for 300 records then CRLF | 3.3 | cpp | 0.12 | 14 | 0 | truth =, PAL2NAL =, engines =; FINDING 4: PAL2NAL reads this file |
| fmt_mixed_lf_then_crlf | 2,000 x 400, LF for 300 records then CRLF | 3.3 | PAL2NAL v14 | 3.05 | 30 | 0 | output |
| many_100k | 100,000 x 400, gaps | 158.6 | python | 3.40 | 444 | 0 | truth =, PAL2NAL =, engines = |
| many_100k | 100,000 x 400, gaps | 158.6 | cpp | 1.52 | 35 | 0 | truth =, PAL2NAL =, engines = |
| many_100k | 100,000 x 400, gaps | 158.6 | PAL2NAL v14 | 195.45 | 1,626 | 0 | output |
| many_1M | 1,000,000 x 400, gaps | 1,587.8 | python | 34.15 | 1,070 | 0 | truth =, engines = |
| many_1M | 1,000,000 x 400, gaps | 1,587.8 | cpp | 16.88 | 201 | 0 | truth =, engines = |
| many_2M_short | 2,000,000 x 30, gaps | 285.8 | python | 20.64 | 2,073 | 0 | truth =, engines = |
| many_2M_short | 2,000,000 x 30, gaps | 285.8 | cpp | 18.10 | 401 | 0 | truth =, engines = |
| many_1M_j1 | 1,000,000 x 400, -j 1 | 1,587.8 | python | 84.99 | 347 | 0 | truth =, engines = |
| many_1M_j1 | 1,000,000 x 400, -j 1 | 1,587.8 | cpp | 14.84 | 208 | 0 | truth =, engines = |
| many_1M_j2 | 1,000,000 x 400, -j 2 | 1,587.8 | python | 45.92 | 473 | 0 | truth =, engines = |
| many_1M_j2 | 1,000,000 x 400, -j 2 | 1,587.8 | cpp | 15.35 | 200 | 0 | truth =, engines = |
| many_1M_j4 | 1,000,000 x 400, -j 4 | 1,587.8 | python | 32.46 | 532 | 0 | truth =, engines = |
| many_1M_j4 | 1,000,000 x 400, -j 4 | 1,587.8 | cpp | 14.46 | 201 | 0 | truth =, engines = |
| many_1M_j8 | 1,000,000 x 400, -j 8 | 1,587.8 | python | 32.05 | 678 | 0 | truth =, engines = |
| many_1M_j8 | 1,000,000 x 400, -j 8 | 1,587.8 | cpp | 15.64 | 201 | 0 | truth =, engines = |
| many_1M_j16 | 1,000,000 x 400, -j 16 | 1,587.8 | python | 34.33 | 1,077 | 0 | truth =, engines = |
| many_1M_j16 | 1,000,000 x 400, -j 16 | 1,587.8 | cpp | 17.41 | 207 | 0 | truth =, engines = |
| many_1M_j32 | 1,000,000 x 400, -j 32 (2x cores) | 1,587.8 | python | 33.85 | 1,875 | 0 | truth =, engines = |
| many_1M_j32 | 1,000,000 x 400, -j 32 (2x cores) | 1,587.8 | cpp | 18.14 | 208 | 0 | truth =, engines = |
| many_1M_stdout | 1,000,000 x 400, stdout | 1,587.8 | python | 32.40 | 1,107 | 0 | truth =, engines = |
| many_1M_stdout | 1,000,000 x 400, stdout | 1,587.8 | cpp | 17.01 | 200 | 0 | truth =, engines = |
| many_1M_clustal | 1,000,000 x 400, clustal out | 1,587.8 | python | 53.59 | 1,088 | 0 | engines = |
| many_1M_clustal | 1,000,000 x 400, clustal out | 1,587.8 | cpp | 23.39 | 225 | 0 | engines = |
| many_1M_paml | 1,000,000 x 400, paml out | 1,587.8 | python | 33.09 | 1,080 | 0 | engines = |
| many_1M_paml | 1,000,000 x 400, paml out | 1,587.8 | cpp | 16.98 | 200 | 0 | engines = |
| many_1M_codon | 1,000,000 x 400, codon out | 1,587.8 | python | 109.38 | 1,087 | 0 | engines = |
| many_1M_codon | 1,000,000 x 400, codon out | 1,587.8 | cpp | 25.50 | 232 | 0 | engines = |
| many_2M_short_j1 | 2,000,000 x 30, -j 1 | 285.8 | python | 34.93 | 648 | 0 | truth =, engines = |
| many_2M_short_j1 | 2,000,000 x 30, -j 1 | 285.8 | cpp | 15.85 | 399 | 0 | truth =, engines = |
| many_2M_short_j8 | 2,000,000 x 30, -j 8 | 285.8 | python | 19.22 | 1,296 | 0 | truth =, engines = |
| many_2M_short_j8 | 2,000,000 x 30, -j 8 | 285.8 | cpp | 16.25 | 408 | 0 | truth =, engines = |
| many_100k_shuffled | 100,000 x 400, DNA records shuffled | 158.6 | python | 6.00 | 443 | 0 | truth =, engines = |
| many_100k_shuffled | 100,000 x 400, DNA records shuffled | 158.6 | cpp | 1.95 | 35 | 0 | truth =, engines = |
| files_1 | 100,000 x 400, 1 DNA file | 158.6 | python | 3.43 | 446 | 0 | truth =, engines = |
| files_1 | 100,000 x 400, 1 DNA file | 158.6 | cpp | 1.56 | 34 | 0 | truth =, engines = |
| files_10 | 100,000 x 400, 10 DNA files | 158.6 | python | 3.52 | 421 | 0 | truth =, engines = |
| files_10 | 100,000 x 400, 10 DNA files | 158.6 | cpp | 1.74 | 35 | 0 | truth =, engines = |
| files_1000 | 100,000 x 400, 1,000 DNA files | 158.6 | python | 2.27 | 179 | 1 | engines DIFF, TRACEBACK; FINDING 3: open-file limit (Python engine with workers) |
| files_1000 | 100,000 x 400, 1,000 DNA files | 158.6 | cpp | 2.42 | 45 | 0 | truth =, engines DIFF; FINDING 3: open-file limit (Python engine with workers) |
| files_1000_clustal | 100,000 x 400, 1,000 DNA files, clustal out | 158.6 | python | 2.90 | 180 | 1 | engines DIFF, TRACEBACK; FINDING 3: open-file limit (Python engine with workers) |
| files_1000_clustal | 100,000 x 400, 1,000 DNA files, clustal out | 158.6 | cpp | 2.37 | 69 | 0 | engines DIFF; FINDING 3: open-file limit (Python engine with workers) |
| stdin_dna | 100,000 x 400, DNA from stdin | 158.6 | python | 3.48 | 449 | 0 | truth =, engines = |
| stdin_dna | 100,000 x 400, DNA from stdin | 158.6 | cpp | 1.63 | 34 | 0 | truth =, engines = |
| stdin_pep | 100,000 x 400, protein from stdin | 158.6 | python | 3.37 | 448 | 0 | truth =, engines = |
| stdin_pep | 100,000 x 400, protein from stdin | 158.6 | cpp | 1.59 | 34 | 0 | truth =, engines = |
| stdin_gzip | 100,000 x 400, gzip -dc | coal pep - | 158.6 | python | 3.79 | 446 | 0 | truth =, engines = |
| stdin_gzip | 100,000 x 400, gzip -dc | coal pep - | 158.6 | cpp | 1.90 | 34 | 0 | truth =, engines = |
| stdin_dna_stdout | 100,000 x 400, DNA from stdin, stdout | 158.6 | python | 3.50 | 443 | 0 | truth =, engines = |
| stdin_dna_stdout | 100,000 x 400, DNA from stdin, stdout | 158.6 | cpp | 1.73 | 37 | 0 | truth =, engines = |
| files_70000 | 70,000 x 50, 70,000 DNA files | 4.3 | python | 0.32 | 18 | 1 | engines =; FINDING 3: open-file limit |
| files_70000 | 70,000 x 50, 70,000 DNA files | 4.3 | cpp | 0.26 | 10 | 1 | engines =; FINDING 3: open-file limit |
| files_5000 | 5,000 x 50, 5,000 DNA files | 0.3 | python | 0.21 | 13 | 1 | engines =; FINDING 3: open-file limit (both engines, clean error) |
| files_5000 | 5,000 x 50, 5,000 DNA files | 0.3 | cpp | 0.23 | 10 | 1 | engines =; FINDING 3: open-file limit (both engines, clean error) |
| files_5000_strict_j1 | 5,000 x 50, 5,000 DNA files, -j 1 | 0.3 | python | 0.21 | 13 | 1 | engines =; FINDING 3: open-file limit (both engines, clean error) |
| files_5000_strict_j1 | 5,000 x 50, 5,000 DNA files, -j 1 | 0.3 | cpp | 0.24 | 8 | 1 | engines =; FINDING 3: open-file limit (both engines, clean error) |
| files_5000_default_ulimit | 5,000 x 50, 5,000 DNA files, ulimit -n 1024 | 1.1 | python | 1.00 | 147 | 1 | engines =; FINDING 3: open-file limit (both engines, clean error) |
| files_5000_default_ulimit | 5,000 x 50, 5,000 DNA files, ulimit -n 1024 | 1.1 | cpp | 1.57 | 17 | 1 | engines =; FINDING 3: open-file limit (both engines, clean error) |
| files_5000_j1_default_ulimit | 5,000 x 50, 5,000 DNA files, -j 1 | 1.1 | python | 0.95 | 136 | 1 | engines =; FINDING 3: open-file limit (both engines, clean error) |
| files_5000_j1_default_ulimit | 5,000 x 50, 5,000 DNA files, -j 1 | 1.1 | cpp | 1.50 | 13 | 1 | engines =; FINDING 3: open-file limit (both engines, clean error) |
| files_70000_ulimit_65536 | 70,000 x 50, 70,000 DNA files, ulimit -n 65536 | 15.7 | python | 16.53 | 70 | 1 | engines DIFF, TRACEBACK; FINDING 2 (Python, fixed) / FINDING 3 (both) |
| files_70000_ulimit_65536 | 70,000 x 50, 70,000 DNA files, ulimit -n 65536 | 15.7 | cpp | 201.79 | 669 | 1 | engines DIFF; FINDING 2 (Python, fixed) / FINDING 3 (both) |
| files_60000_ulimit_65536 | 60,000 x 50, 60,000 DNA files, ulimit -n 65536 | 13.4 | python | 92.64 | 8,083 | 0 | truth =, engines =; FINDING 3: 129 KB per open file (Python) |
| files_60000_ulimit_65536 | 60,000 x 50, 60,000 DNA files, ulimit -n 65536 | 13.4 | cpp | 152.12 | 663 | 0 | truth =, engines =; FINDING 3: 129 KB per open file (Python) |
| opt_100k_fasta | 100,000 x 400, -output fasta | 158.6 | python | 3.48 | 451 | 0 | truth =, PAL2NAL =, engines = |
| opt_100k_fasta | 100,000 x 400, -output fasta | 158.6 | cpp | 1.55 | 34 | 0 | truth =, PAL2NAL =, engines = |
| opt_100k_fasta | 100,000 x 400, -output fasta | 158.6 | PAL2NAL v14 | 184.41 | 1,626 | 0 | output |
| opt_100k_clustal | 100,000 x 400, -output clustal | 158.6 | python | 5.46 | 422 | 0 | PAL2NAL =, engines = |
| opt_100k_clustal | 100,000 x 400, -output clustal | 158.6 | cpp | 1.77 | 59 | 0 | PAL2NAL =, engines = |
| opt_100k_clustal | 100,000 x 400, -output clustal | 158.6 | PAL2NAL v14 | 193.38 | 1,626 | 0 | output |
| opt_100k_paml | 100,000 x 400, -output paml | 158.6 | python | 3.42 | 449 | 0 | PAL2NAL =, engines = |
| opt_100k_paml | 100,000 x 400, -output paml | 158.6 | cpp | 1.54 | 35 | 0 | PAL2NAL =, engines = |
| opt_100k_paml | 100,000 x 400, -output paml | 158.6 | PAL2NAL v14 | 182.94 | 1,626 | 0 | output |
| opt_100k_codon | 100,000 x 400, -output codon | 158.6 | python | 11.55 | 426 | 0 | PAL2NAL =, engines = |
| opt_100k_codon | 100,000 x 400, -output codon | 158.6 | cpp | 2.39 | 69 | 0 | PAL2NAL =, engines = |
| opt_100k_codon | 100,000 x 400, -output codon | 158.6 | PAL2NAL v14 | 237.77 | 1,798 | 0 | output |
| opt_100k_nogap | 100,000 x 400, -nogap | 158.6 | python | 3.95 | 451 | 0 | PAL2NAL =, engines = |
| opt_100k_nogap | 100,000 x 400, -nogap | 158.6 | cpp | 1.69 | 33 | 0 | PAL2NAL =, engines = |
| opt_100k_nogap | 100,000 x 400, -nogap | 158.6 | PAL2NAL v14 | 202.60 | 1,156 | 0 | output |
| opt_100k_nomismatch | 100,000 x 400, -nomismatch | 158.6 | python | 4.09 | 446 | 0 | PAL2NAL =, engines = |
| opt_100k_nomismatch | 100,000 x 400, -nomismatch | 158.6 | cpp | 1.66 | 33 | 0 | PAL2NAL =, engines = |
| opt_100k_nomismatch | 100,000 x 400, -nomismatch | 158.6 | PAL2NAL v14 | 182.28 | 1,626 | 0 | output |
| opt_100k_nogap_nomismatch_clustal | 100,000 x 400, -nogap -nomismatch clustal | 158.6 | python | 3.64 | 434 | 0 | PAL2NAL =, engines = |
| opt_100k_nogap_nomismatch_clustal | 100,000 x 400, -nogap -nomismatch clustal | 158.6 | cpp | 1.55 | 36 | 0 | PAL2NAL =, engines = |
| opt_100k_nogap_nomismatch_clustal | 100,000 x 400, -nogap -nomismatch clustal | 158.6 | PAL2NAL v14 | 178.07 | 1,156 | 0 | output |
| opt_100k_blockonly | 100,000 x 400, CLUSTAL in with marks, -blockonly | 164.6 | python | 5.14 | 395 | 0 | PAL2NAL =, engines = |
| opt_100k_blockonly | 100,000 x 400, CLUSTAL in with marks, -blockonly | 164.6 | cpp | 2.18 | 144 | 0 | PAL2NAL =, engines = |
| opt_100k_blockonly | 100,000 x 400, CLUSTAL in with marks, -blockonly | 164.6 | PAL2NAL v14 | 152.46 | 1,246 | 0 | output |
| opt_100k_marks_clustal_out | 100,000 x 400, CLUSTAL in with marks, clustal out | 164.6 | python | 6.02 | 358 | 0 | PAL2NAL =, engines = |
| opt_100k_marks_clustal_out | 100,000 x 400, CLUSTAL in with marks, clustal out | 164.6 | cpp | 2.27 | 168 | 0 | PAL2NAL =, engines = |
| opt_100k_marks_clustal_out | 100,000 x 400, CLUSTAL in with marks, clustal out | 164.6 | PAL2NAL v14 | 179.18 | 1,648 | 0 | output |
| opt_100k_strict | 100,000 x 400, --strict | 158.6 | python | 3.16 | 452 | 0 | truth =, engines = |
| opt_100k_strict | 100,000 x 400, --strict | 158.6 | cpp | 1.42 | 36 | 0 | truth =, engines = |
| opt_100k_keepgoing | 100,000 x 400, --keep-going | 158.6 | python | 3.17 | 448 | 0 | truth =, engines = |
| opt_100k_keepgoing | 100,000 x 400, --keep-going | 158.6 | cpp | 1.42 | 34 | 0 | truth =, engines = |
| opt_100k_report | 100,000 x 400, --report | 158.6 | python | 3.66 | 451 | 0 | truth =, engines = |
| opt_100k_report | 100,000 x 400, --report | 158.6 | cpp | 2.47 | 38 | 0 | truth =, engines = |
| opt_100k_json | 100,000 x 400, --json | 158.6 | python | 3.66 | 445 | 0 | engines = |
| opt_100k_json | 100,000 x 400, --json | 158.6 | cpp | 2.30 | 37 | 0 | engines = |
| opt_100k_progress | 100,000 x 400, --progress | 158.6 | python | 3.13 | 444 | 0 | truth =, engines DIFF; --progress lines are timing-based (expected) |
| opt_100k_progress | 100,000 x 400, --progress | 158.6 | cpp | 1.43 | 36 | 0 | truth =, engines DIFF; --progress lines are timing-based (expected) |
| opt_100k_maxwarn0 | 100,000 x 400, --max-warnings 0 | 158.6 | python | 3.08 | 452 | 0 | truth =, engines = |
| opt_100k_maxwarn0 | 100,000 x 400, --max-warnings 0 | 158.6 | cpp | 1.41 | 36 | 0 | truth =, engines = |
| opt_100k_quiet | 100,000 x 400, -q | 158.6 | python | 3.12 | 451 | 0 | truth =, engines = |
| opt_100k_quiet | 100,000 x 400, -q | 158.6 | cpp | 1.39 | 35 | 0 | truth =, engines = |
| opt_1M_report | 1,000,000 x 400, --report | 1,587.8 | python | 35.42 | 1,075 | 0 | truth =, engines = |
| opt_1M_report | 1,000,000 x 400, --report | 1,587.8 | cpp | 24.52 | 210 | 0 | truth =, engines = |
| opt_1M_progress_strict | 1,000,000 x 400, --progress --strict | 1,587.8 | python | 29.68 | 1,075 | 0 | truth =, engines DIFF; --progress lines are timing-based (expected) |
| opt_1M_progress_strict | 1,000,000 x 400, --progress --strict | 1,587.8 | cpp | 14.93 | 205 | 0 | truth =, engines DIFF; --progress lines are timing-based (expected) |
| codontable_1 | 20,000 x 400, -codontable 1 | 31.7 | python | 1.13 | 135 | 0 | truth =, PAL2NAL =, engines = |
| codontable_1 | 20,000 x 400, -codontable 1 | 31.7 | cpp | 0.36 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_1 | 20,000 x 400, -codontable 1 | 31.7 | PAL2NAL v14 | 28.73 | 327 | 0 | output |
| codontable_1_nogap | 20,000 x 400, -codontable 1 -nogap | 31.7 | python | 0.94 | 152 | 0 | PAL2NAL =, engines = |
| codontable_1_nogap | 20,000 x 400, -codontable 1 -nogap | 31.7 | cpp | 0.39 | 16 | 0 | PAL2NAL =, engines = |
| codontable_1_nogap | 20,000 x 400, -codontable 1 -nogap | 31.7 | PAL2NAL v14 | 30.79 | 232 | 0 | output |
| codontable_2 | 20,000 x 400, -codontable 2 | 31.7 | python | 0.84 | 153 | 0 | truth =, PAL2NAL =, engines = |
| codontable_2 | 20,000 x 400, -codontable 2 | 31.7 | cpp | 0.36 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_2 | 20,000 x 400, -codontable 2 | 31.7 | PAL2NAL v14 | 27.07 | 328 | 0 | output |
| codontable_2_nogap | 20,000 x 400, -codontable 2 -nogap | 31.7 | python | 0.95 | 145 | 0 | PAL2NAL =, engines = |
| codontable_2_nogap | 20,000 x 400, -codontable 2 -nogap | 31.7 | cpp | 0.38 | 17 | 0 | PAL2NAL =, engines = |
| codontable_2_nogap | 20,000 x 400, -codontable 2 -nogap | 31.7 | PAL2NAL v14 | 29.36 | 232 | 0 | output |
| codontable_3 | 20,000 x 400, -codontable 3 | 31.7 | python | 0.82 | 148 | 0 | truth =, PAL2NAL =, engines = |
| codontable_3 | 20,000 x 400, -codontable 3 | 31.7 | cpp | 0.37 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_3 | 20,000 x 400, -codontable 3 | 31.7 | PAL2NAL v14 | 27.61 | 328 | 0 | output |
| codontable_3_nogap | 20,000 x 400, -codontable 3 -nogap | 31.7 | python | 0.94 | 147 | 0 | PAL2NAL =, engines = |
| codontable_3_nogap | 20,000 x 400, -codontable 3 -nogap | 31.7 | cpp | 0.38 | 16 | 0 | PAL2NAL =, engines = |
| codontable_3_nogap | 20,000 x 400, -codontable 3 -nogap | 31.7 | PAL2NAL v14 | 30.82 | 232 | 0 | output |
| codontable_4 | 20,000 x 400, -codontable 4 | 31.7 | python | 0.86 | 143 | 0 | truth =, PAL2NAL =, engines = |
| codontable_4 | 20,000 x 400, -codontable 4 | 31.7 | cpp | 0.37 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_4 | 20,000 x 400, -codontable 4 | 31.7 | PAL2NAL v14 | 30.39 | 328 | 0 | output |
| codontable_4_nogap | 20,000 x 400, -codontable 4 -nogap | 31.7 | python | 0.94 | 145 | 0 | PAL2NAL =, engines = |
| codontable_4_nogap | 20,000 x 400, -codontable 4 -nogap | 31.7 | cpp | 0.38 | 16 | 0 | PAL2NAL =, engines = |
| codontable_4_nogap | 20,000 x 400, -codontable 4 -nogap | 31.7 | PAL2NAL v14 | 30.99 | 232 | 0 | output |
| codontable_5 | 20,000 x 400, -codontable 5 | 31.7 | python | 0.83 | 150 | 0 | truth =, PAL2NAL =, engines = |
| codontable_5 | 20,000 x 400, -codontable 5 | 31.7 | cpp | 0.37 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_5 | 20,000 x 400, -codontable 5 | 31.7 | PAL2NAL v14 | 26.40 | 328 | 0 | output |
| codontable_5_nogap | 20,000 x 400, -codontable 5 -nogap | 31.7 | python | 0.94 | 148 | 0 | PAL2NAL =, engines = |
| codontable_5_nogap | 20,000 x 400, -codontable 5 -nogap | 31.7 | cpp | 0.40 | 16 | 0 | PAL2NAL =, engines = |
| codontable_5_nogap | 20,000 x 400, -codontable 5 -nogap | 31.7 | PAL2NAL v14 | 28.34 | 232 | 0 | output |
| codontable_6 | 20,000 x 400, -codontable 6 | 31.7 | python | 0.86 | 149 | 0 | truth =, PAL2NAL =, engines = |
| codontable_6 | 20,000 x 400, -codontable 6 | 31.7 | cpp | 0.36 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_6 | 20,000 x 400, -codontable 6 | 31.7 | PAL2NAL v14 | 30.29 | 328 | 0 | output |
| codontable_6_nogap | 20,000 x 400, -codontable 6 -nogap | 31.7 | python | 0.94 | 151 | 0 | PAL2NAL =, engines = |
| codontable_6_nogap | 20,000 x 400, -codontable 6 -nogap | 31.7 | cpp | 0.38 | 16 | 0 | PAL2NAL =, engines = |
| codontable_6_nogap | 20,000 x 400, -codontable 6 -nogap | 31.7 | PAL2NAL v14 | 33.62 | 232 | 0 | output |
| codontable_9 | 20,000 x 400, -codontable 9 | 31.7 | python | 0.86 | 142 | 0 | truth =, PAL2NAL =, engines = |
| codontable_9 | 20,000 x 400, -codontable 9 | 31.7 | cpp | 0.36 | 20 | 0 | truth =, PAL2NAL =, engines = |
| codontable_9 | 20,000 x 400, -codontable 9 | 31.7 | PAL2NAL v14 | 26.31 | 328 | 0 | output |
| codontable_9_nogap | 20,000 x 400, -codontable 9 -nogap | 31.7 | python | 0.94 | 150 | 0 | PAL2NAL =, engines = |
| codontable_9_nogap | 20,000 x 400, -codontable 9 -nogap | 31.7 | cpp | 0.38 | 16 | 0 | PAL2NAL =, engines = |
| codontable_9_nogap | 20,000 x 400, -codontable 9 -nogap | 31.7 | PAL2NAL v14 | 28.59 | 232 | 0 | output |
| codontable_10 | 20,000 x 400, -codontable 10 | 31.7 | python | 0.89 | 142 | 0 | truth =, PAL2NAL =, engines = |
| codontable_10 | 20,000 x 400, -codontable 10 | 31.7 | cpp | 0.34 | 18 | 0 | truth =, PAL2NAL =, engines = |
| codontable_10 | 20,000 x 400, -codontable 10 | 31.7 | PAL2NAL v14 | 107.14 | 363 | 0 | output |
| codontable_10_nogap | 20,000 x 400, -codontable 10 -nogap | 31.7 | python | 0.94 | 144 | 0 | PAL2NAL =, engines = |
| codontable_10_nogap | 20,000 x 400, -codontable 10 -nogap | 31.7 | cpp | 0.38 | 16 | 0 | PAL2NAL =, engines = |
| codontable_10_nogap | 20,000 x 400, -codontable 10 -nogap | 31.7 | PAL2NAL v14 | 109.97 | 258 | 0 | output |
| codontable_11 | 20,000 x 400, -codontable 11 | 31.7 | python | 0.82 | 144 | 0 | truth =, PAL2NAL =, engines = |
| codontable_11 | 20,000 x 400, -codontable 11 | 31.7 | cpp | 0.36 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_11 | 20,000 x 400, -codontable 11 | 31.7 | PAL2NAL v14 | 28.50 | 328 | 0 | output |
| codontable_11_nogap | 20,000 x 400, -codontable 11 -nogap | 31.7 | python | 0.97 | 151 | 0 | PAL2NAL =, engines = |
| codontable_11_nogap | 20,000 x 400, -codontable 11 -nogap | 31.7 | cpp | 0.38 | 16 | 0 | PAL2NAL =, engines = |
| codontable_11_nogap | 20,000 x 400, -codontable 11 -nogap | 31.7 | PAL2NAL v14 | 31.21 | 233 | 0 | output |
| codontable_12 | 20,000 x 400, -codontable 12 | 31.7 | python | 0.87 | 144 | 0 | truth =, PAL2NAL =, engines = |
| codontable_12 | 20,000 x 400, -codontable 12 | 31.7 | cpp | 0.37 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_12 | 20,000 x 400, -codontable 12 | 31.7 | PAL2NAL v14 | 31.27 | 328 | 0 | output |
| codontable_12_nogap | 20,000 x 400, -codontable 12 -nogap | 31.7 | python | 0.98 | 145 | 0 | PAL2NAL =, engines = |
| codontable_12_nogap | 20,000 x 400, -codontable 12 -nogap | 31.7 | cpp | 0.40 | 16 | 0 | PAL2NAL =, engines = |
| codontable_12_nogap | 20,000 x 400, -codontable 12 -nogap | 31.7 | PAL2NAL v14 | 33.48 | 231 | 0 | output |
| codontable_13 | 20,000 x 400, -codontable 13 | 31.7 | python | 0.86 | 141 | 0 | truth =, PAL2NAL =, engines = |
| codontable_13 | 20,000 x 400, -codontable 13 | 31.7 | cpp | 0.37 | 18 | 0 | truth =, PAL2NAL =, engines = |
| codontable_13 | 20,000 x 400, -codontable 13 | 31.7 | PAL2NAL v14 | 31.16 | 328 | 0 | output |
| codontable_13_nogap | 20,000 x 400, -codontable 13 -nogap | 31.7 | python | 1.05 | 153 | 0 | PAL2NAL =, engines = |
| codontable_13_nogap | 20,000 x 400, -codontable 13 -nogap | 31.7 | cpp | 0.39 | 16 | 0 | PAL2NAL =, engines = |
| codontable_13_nogap | 20,000 x 400, -codontable 13 -nogap | 31.7 | PAL2NAL v14 | 33.59 | 232 | 0 | output |
| codontable_14 | 20,000 x 400, -codontable 14 | 31.7 | python | 0.85 | 145 | 0 | truth =, PAL2NAL =, engines = |
| codontable_14 | 20,000 x 400, -codontable 14 | 31.7 | cpp | 0.36 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_14 | 20,000 x 400, -codontable 14 | 31.7 | PAL2NAL v14 | 26.65 | 328 | 0 | output |
| codontable_14_nogap | 20,000 x 400, -codontable 14 -nogap | 31.7 | python | 0.95 | 149 | 0 | PAL2NAL =, engines = |
| codontable_14_nogap | 20,000 x 400, -codontable 14 -nogap | 31.7 | cpp | 0.42 | 18 | 0 | PAL2NAL =, engines = |
| codontable_14_nogap | 20,000 x 400, -codontable 14 -nogap | 31.7 | PAL2NAL v14 | 30.09 | 232 | 0 | output |
| codontable_15 | 20,000 x 400, -codontable 15 | 31.7 | python | 0.86 | 146 | 0 | truth =, PAL2NAL =, engines = |
| codontable_15 | 20,000 x 400, -codontable 15 | 31.7 | cpp | 0.37 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_15 | 20,000 x 400, -codontable 15 | 31.7 | PAL2NAL v14 | 29.57 | 328 | 0 | output |
| codontable_15_nogap | 20,000 x 400, -codontable 15 -nogap | 31.7 | python | 0.92 | 151 | 0 | PAL2NAL =, engines = |
| codontable_15_nogap | 20,000 x 400, -codontable 15 -nogap | 31.7 | cpp | 0.38 | 16 | 0 | PAL2NAL =, engines = |
| codontable_15_nogap | 20,000 x 400, -codontable 15 -nogap | 31.7 | PAL2NAL v14 | 31.73 | 232 | 0 | output |
| codontable_16 | 20,000 x 400, -codontable 16 | 31.7 | python | 0.85 | 143 | 0 | truth =, PAL2NAL =, engines = |
| codontable_16 | 20,000 x 400, -codontable 16 | 31.7 | cpp | 0.36 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_16 | 20,000 x 400, -codontable 16 | 31.7 | PAL2NAL v14 | 31.49 | 328 | 0 | output |
| codontable_16_nogap | 20,000 x 400, -codontable 16 -nogap | 31.7 | python | 0.96 | 146 | 0 | PAL2NAL =, engines = |
| codontable_16_nogap | 20,000 x 400, -codontable 16 -nogap | 31.7 | cpp | 0.40 | 16 | 0 | PAL2NAL =, engines = |
| codontable_16_nogap | 20,000 x 400, -codontable 16 -nogap | 31.7 | PAL2NAL v14 | 33.35 | 231 | 0 | output |
| codontable_21 | 20,000 x 400, -codontable 21 | 31.7 | python | 0.87 | 147 | 0 | truth =, PAL2NAL =, engines = |
| codontable_21 | 20,000 x 400, -codontable 21 | 31.7 | cpp | 0.36 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_21 | 20,000 x 400, -codontable 21 | 31.7 | PAL2NAL v14 | 27.09 | 328 | 0 | output |
| codontable_21_nogap | 20,000 x 400, -codontable 21 -nogap | 31.7 | python | 0.95 | 149 | 0 | PAL2NAL =, engines = |
| codontable_21_nogap | 20,000 x 400, -codontable 21 -nogap | 31.7 | cpp | 0.39 | 16 | 0 | PAL2NAL =, engines = |
| codontable_21_nogap | 20,000 x 400, -codontable 21 -nogap | 31.7 | PAL2NAL v14 | 29.14 | 232 | 0 | output |
| codontable_22 | 20,000 x 400, -codontable 22 | 31.7 | python | 0.89 | 149 | 0 | truth =, PAL2NAL =, engines = |
| codontable_22 | 20,000 x 400, -codontable 22 | 31.7 | cpp | 0.40 | 18 | 0 | truth =, PAL2NAL =, engines = |
| codontable_22 | 20,000 x 400, -codontable 22 | 31.7 | PAL2NAL v14 | 31.83 | 328 | 0 | output |
| codontable_22_nogap | 20,000 x 400, -codontable 22 -nogap | 31.7 | python | 0.95 | 147 | 0 | PAL2NAL =, engines = |
| codontable_22_nogap | 20,000 x 400, -codontable 22 -nogap | 31.7 | cpp | 0.38 | 16 | 0 | PAL2NAL =, engines = |
| codontable_22_nogap | 20,000 x 400, -codontable 22 -nogap | 31.7 | PAL2NAL v14 | 34.26 | 232 | 0 | output |
| codontable_23 | 20,000 x 400, -codontable 23 | 31.7 | python | 0.84 | 151 | 0 | truth =, PAL2NAL =, engines = |
| codontable_23 | 20,000 x 400, -codontable 23 | 31.7 | cpp | 0.36 | 19 | 0 | truth =, PAL2NAL =, engines = |
| codontable_23 | 20,000 x 400, -codontable 23 | 31.7 | PAL2NAL v14 | 27.58 | 328 | 0 | output |
| codontable_23_nogap | 20,000 x 400, -codontable 23 -nogap | 31.7 | python | 0.98 | 148 | 0 | PAL2NAL =, engines = |
| codontable_23_nogap | 20,000 x 400, -codontable 23 -nogap | 31.7 | cpp | 0.38 | 16 | 0 | PAL2NAL =, engines = |
| codontable_23_nogap | 20,000 x 400, -codontable 23 -nogap | 31.7 | PAL2NAL v14 | 30.04 | 232 | 0 | output |
| keepgoing_100k_one_empty | 100,000 x 400, DNA record 777 empty, --keep-going | 158.6 | python | 3.77 | 463 | 2 | truth DIFF, engines =; by design: the empty gene is left out (exit 2) |
| keepgoing_100k_one_empty | 100,000 x 400, DNA record 777 empty, --keep-going | 158.6 | cpp | 1.64 | 38 | 2 | truth DIFF, engines =; by design: the empty gene is left out (exit 2) |
| keepgoing_100k_one_empty_paml | 100,000 x 400, DNA record 777 empty, --keep-going paml | 158.6 | python | 4.24 | 459 | 2 | engines =; by design: the empty gene is left out (exit 2) |
| keepgoing_100k_one_empty_paml | 100,000 x 400, DNA record 777 empty, --keep-going paml | 158.6 | cpp | 1.67 | 37 | 2 | engines =; by design: the empty gene is left out (exit 2) |
| strict_100k_one_empty | 100,000 x 400, DNA record 777 empty (error) | 158.6 | python | 2.97 | 400 | 1 | engines =; error by design: a gene without DNA |
| strict_100k_one_empty | 100,000 x 400, DNA record 777 empty (error) | 158.6 | cpp | 0.65 | 25 | 1 | engines =; error by design: a gene without DNA |
| resume_cpp_output_fasta | 1,000,000 x 400, -output fasta, kill/resume rounds |  | cpp | 13.06 (ref) |  | 16 runs | OK: byte-identical after every kill/resume, no leftovers |
| resume_python_output_fasta | 1,000,000 x 400, -output fasta, kill/resume rounds |  | python | 28.84 (ref) |  | 13 runs | OK: byte-identical after every kill/resume, no leftovers |
| resume_cpp_output_clustal | 1,000,000 x 400, -output clustal, kill/resume rounds |  | cpp | 15.07 (ref) |  | 13 runs | OK: byte-identical after every kill/resume, no leftovers |
| resume_python_output_clustal | 1,000,000 x 400, -output clustal, kill/resume rounds |  | python | 49.82 (ref) |  | 7 runs | OK: byte-identical after every kill/resume, no leftovers |
| resume_cpp_output_fasta_nogap_report_rep_resume.json | 1,000,000 x 400, -output fasta -nogap --report rep_resume.json, kill/resume rounds |  | cpp | 20.7 (ref) |  | 11 runs | OK: byte-identical after every kill/resume, no leftovers |
| resume_python_output_paml_keep-going | 1,000,000 x 400, -output paml --keep-going, kill/resume rounds |  | python | 33.62 (ref) |  | 11 runs | OK: byte-identical after every kill/resume, no leftovers |
| resume_cpp_output_fasta | 1,000,000 x 400, -output fasta, corruption rounds |  | cpp | 13.86 (ref) |  | 12 runs | OK: every damaged checkpoint refused as documented, output byte-identical |
| resume_python_output_fasta | 1,000,000 x 400, -output fasta, corruption rounds |  | python | 34.82 (ref) |  | 12 runs | OK: 1 corruption round(s) killed before the first checkpoint (nothing to corrupt; normal run, output byte-identical) |
| resume_cpp_output_clustal | 1,000,000 x 400, -output clustal, corruption rounds |  | cpp | 15.52 (ref) |  | 12 runs | OK: 1 corruption round(s) killed before the first checkpoint (nothing to corrupt; normal run, output byte-identical) |
| resume_cpp_output_fasta | 1,000,000 x 400, -output fasta, corruption rounds |  | cpp | 14.14 (ref) |  | 12 runs | OK: every damaged checkpoint refused as documented, output byte-identical |
| resume_python_output_fasta | 1,000,000 x 400, -output fasta, corruption rounds |  | python | 31.58 (ref) |  | 12 runs | OK: 1 corruption round(s) killed before the first checkpoint (nothing to corrupt; normal run, output byte-identical) |
| resume_cpp_output_clustal | 1,000,000 x 400, -output clustal, corruption rounds |  | cpp | 19.74 (ref) |  | 12 runs | OK: every damaged checkpoint refused as documented, output byte-identical |
| resume_python_output_clustal | 1,000,000 x 400, -output clustal, corruption rounds |  | python | 61.85 (ref) |  | 12 runs | OK: every damaged checkpoint refused as documented, output byte-identical |
| wrongcode_100k | 100,000 x 400 built with code 2, run with code 1 | 164.7 | python | 21.81 | 996 | 0 | truth =, engines = |
| wrongcode_100k | 100,000 x 400 built with code 2, run with code 1 | 164.7 | cpp | 14.10 | 352 | 0 | truth =, engines = |
| wrongcode_100k_maxwarn0 | 100,000 x 400 wrong code, --max-warnings 0 | 164.7 | python | 37.57 | 1,046 | 0 | truth =, engines = |
| wrongcode_100k_maxwarn0 | 100,000 x 400 wrong code, --max-warnings 0 | 164.7 | cpp | 39.18 | 352 | 0 | truth =, engines = |
| wrongcode_100k_report | 100,000 x 400 wrong code, --report | 164.7 | python | 22.52 | 1,001 | 0 | truth =, engines = |
| wrongcode_100k_report | 100,000 x 400 wrong code, --report | 164.7 | cpp | 17.05 | 352 | 0 | truth =, engines = |
| wrongcode_100k_nomismatch | 100,000 x 400 wrong code, -nomismatch | 164.7 | python | 24.00 | 1,157 | 0 | engines = |
| wrongcode_100k_nomismatch | 100,000 x 400 wrong code, -nomismatch | 164.7 | cpp | 14.47 | 369 | 0 | engines = |
| wrongcode_100k_strict | 100,000 x 400 wrong code, --strict | 164.7 | python | 22.36 | 983 | 0 | truth =, engines = |
| wrongcode_100k_strict | 100,000 x 400 wrong code, --strict | 164.7 | cpp | 13.71 | 353 | 0 | truth =, engines = |
| wrongcode_1M | 1,000,000 x 400 wrong code | 1,648.8 | python | 349.76 | 6,392 | 0 | truth =, engines =; FINDING 9: memory per warning |
| wrongcode_1M | 1,000,000 x 400 wrong code | 1,648.8 | cpp | 217.56 | 3,377 | 0 | truth =, engines =; FINDING 9: memory per warning |
| wrongcode_1M_maxwarn0_report | 1,000,000 x 400 wrong code, --max-warnings 0 --report | 1,648.8 | python | 520.42 | 6,372 | 0 | truth =, engines =; FINDING 9: memory per warning |
| wrongcode_1M_maxwarn0_report | 1,000,000 x 400 wrong code, --max-warnings 0 --report | 1,648.8 | cpp | 650.24 | 3,300 | 0 | truth =, engines =; FINDING 9: memory per warning |
| unrelated_100k | 100,000 x 400 unrelated DNA | 164.7 | python | 537.41 | 14,399 | 0 | truth =, engines =; FINDING 9: memory per warning |
| unrelated_100k | 100,000 x 400 unrelated DNA | 164.7 | cpp | 451.06 | 5,495 | 0 | truth =, engines =; FINDING 9: memory per warning |
| unrelated_100k_strict_keepgoing | 100,000 x 400 unrelated, --strict --keep-going | 164.7 | python | 17.93 | 518 | 1 | engines =; FINDING 9: memory per warning |
| unrelated_100k_strict_keepgoing | 100,000 x 400 unrelated, --strict --keep-going | 164.7 | cpp | 5.80 | 76 | 1 | engines =; FINDING 9: memory per warning |
| unrelated_1M | 1,000,000 x 400 unrelated DNA | 1,648.8 | python | 161.09 | 15,182 | -9 | killed by the memory watchdog at 16 GB; FINDING 9: killed by the 16 GB memory watchdog |
| unrelated_1M | 1,000,000 x 400 unrelated DNA | 1,648.8 | cpp | 1136.85 | 17,182 | -9 | killed by the memory watchdog at 16 GB; FINDING 9: killed by the 16 GB memory watchdog |
| unrelated_1M_report_maxwarn0 | 1,000,000 x 400 unrelated, --report --max-warnings 0 | 1,648.8 | python | 485.03 | 15,049 | -9 | killed by the memory watchdog at 16 GB; FINDING 9: killed by the 16 GB memory watchdog |
| unrelated_1M_report_maxwarn0 | 1,000,000 x 400 unrelated, --report --max-warnings 0 | 1,648.8 | cpp | 427.93 | 2,684 | -9 | killed by the memory watchdog at 16 GB; FINDING 9: killed by the 16 GB memory watchdog |
| normalised_ids_1M | 1,000,000 x 400, IDs g1 vs g1.1, DNA reversed | 1,646.8 | python | 48.47 | 1,058 | 0 | truth =, engines = |
| normalised_ids_1M | 1,000,000 x 400, IDs g1 vs g1.1, DNA reversed | 1,646.8 | cpp | 16.31 | 414 | 0 | truth =, engines = |
| half_ids_1M | 1,000,000 x 400, 50% of DNA IDs differ (by order) | 1,649.8 | python | 38.97 | 1,539 | 0 | truth =, engines = |
| half_ids_1M | 1,000,000 x 400, 50% of DNA IDs differ (by order) | 1,649.8 | cpp | 20.59 | 606 | 0 | truth =, engines = |
| half_ids_1M_strict | 1,000,000 x 400, 50% IDs differ, --strict | 1,649.8 | python | 10.49 | 265 | 1 | engines =; error by design: partial ID match with --strict |
| half_ids_1M_strict | 1,000,000 x 400, 50% IDs differ, --strict | 1,649.8 | cpp | 3.04 | 201 | 1 | engines =; error by design: partial ID match with --strict |
| concurrent_1000_distinct | 2,000 x 400, 1,000 processes, 1,000 -o files in one folder | 2.4 | cpp | 27.9 (all) | | 0 | 1,000 outputs = reference, no temp files left |
| concurrent_1000_distinct | 2,000 x 400, 1,000 processes (400 at a time), 1,000 -o files | 2.4 | python | 52.9 (all) | | 0 | 1,000 outputs = reference, no temp files left |
| concurrent_1000_same | 2,000 x 400, 1,000 processes, one -o file | 2.4 | cpp | 29.5 (all) | | 536 x 1 | file correct; FINDING 1 ("no such file or directory") |
| concurrent_1000_same | 2,000 x 400, 1,000 processes (400 at a time), one -o file | 2.4 | python | 49.7 (all) | | 495 x 1 | file correct; FINDING 1 |
| concurrent_same_o | 100,000 x 400, 4 runs with different options into one -o file, 3 trials | 158.6 | cpp | | | 0 after lost output | FINDING 1: exit 0 with another run's output |
| concurrent_same_o | 100,000 x 400, 4 runs with different options into one -o file, 3 trials | 158.6 | python | | | 0 with an empty file | FINDING 1: exit 0 with an empty output file |
| nohup_sighup | 1,000,000 x 400 under nohup, SIGHUP after 4 s (before the fix) | 1,587.8 | python | | | 129 | FINDING 5: died; fixed in this PR (now 0, output complete) |
| nohup_sighup | 1,000,000 x 400 under nohup, SIGHUP after 4 s (before the fix) | 1,587.8 | cpp | | | 129 | FINDING 5: died; fixed in this PR (now 0, output complete) |
| orphans_head1 | 1,000,000 x 400, `\| head -1`, -j 8 and default | 1,587.8 | python, cpp | | | 141 | no worker left behind, stderr empty |
| sigterm_workers / sigkill_workers | 1,000,000 x 400 -j 8, SIGTERM / SIGKILL to the main process | 1,587.8 | python | | | 143 / 137 | no worker left 4 s later, no partial output left |
| minrepro_mixed_cr | 2 genes, 65,600-byte first record then a lone-CR record | 0.07 | python, cpp | 0.1 | 13 | 1 | FINDING 4 (PAL2NAL: correct output) |
| minrepro_quadratic_20k | 1 x 20,000 poly-A after 20,000 codons of poly-A | 0.1 | python / cpp | 0.15 / 1.54 | 22 / 15 | 0 | identical outputs; FINDING 10 |
| hint_cost_200k | 1 x 200,000, last residue wrong, default / -q | 0.8 | python | 5.65 / 0.40 | 155 / 66 | 0 | FINDING 11 |
| hint_cost_200k | 1 x 200,000, last residue wrong, default / -q | 0.8 | cpp | 1.25 / 0.16 | 63 / 27 | 0 | FINDING 11 |
| devfull / fsize | output to /dev/full, ulimit -f 2048 (WSL 1 enforces neither) | 2.4 | python, cpp | | | 0 | not testable on this machine (see the end) |
| readonly_dir | -o into a read-only folder; --report there; --resume there | 2.4 | python, cpp | | | 1 | `cannot write ...: permission denied (choose a folder you can write to)` |
| readonly_cwd / readonly_as_TMPDIR | clustal to stdout from a read-only cwd / TMPDIR | 2.4 | python, cpp | | | 0 | ran (temp folder falls back) |
| tmpdir_unset / missing / is_file / unwritable | clustal -o with TMPDIR unset, /nonexistent, a file, /proc | 2.4 | python, cpp | | | 0 | identical outputs |
| head1_fasta/clustal/paml/codon/json | 100,000 x 400 `\| head -1` | 158.6 | python, cpp | | | 141 | stderr empty, no traceback |
| out_is_input / out_is_folder / out_symlink_loop / input_is_folder / input_is_dev_zero / two_stdin | 2,000 x 400 | 2.4 | python, cpp | | | 1 / 1 / 1 / 1 / 1 / 2 | one-line errors, identical in both engines |
| win_spaces_o | Windows, 100,000 x 400 in `dir with spaces\sub dir`, -o `out put.fa` | 158.6 | python / cpp | 5.3 / 1.7 | | 0 | truth =, engines = |
| win_forward_slashes / win_stdin / win_spaces_clustal_report | Windows, 100,000 x 400 | 158.6 | python / cpp | 5.6 / 1.8; 5.5 / 2.4; 10.6 / 3.2 | | 0 | truth = |
| win_closed_pipe | Windows, 100,000 x 400, reader closes after one line | 158.6 | python / cpp | 5.8 / 2.0 | | 141 | no traceback |
| win_memlimit_256MB | Windows Job Object 256 MB, 1,000,000 x 400 | 1,587.8 | python | 15.4 | 269 (job) | 1 | FINDING 6: `MemoryError` traceback (also with -j 1) |
| win_memlimit_512MB | Windows Job Object 512 MB, 1,000,000 x 400 | 1,587.8 | python | 37.5 | 537 (job) | 1 | `coal: error: a worker process died (out of memory?)` |
| win_memlimit_512MB_j1 | Windows Job Object 512 MB, 1,000,000 x 400, -j 1 | 1,587.8 | python | 86.0 | 314 (job) | 0 | output complete (1,230,888,890 bytes) |
| win_memlimit_1024MB | Windows Job Object 1 GB, 1,000,000 x 400 | 1,587.8 | python | 41.9 | 551 (job) | 0 | output complete |
| win_memlimit_256MB..1024MB | Windows Job Object 256 MB to 1 GB, 1,000,000 x 400 | 1,587.8 | cpp | 13.7-16.2 | 228 (job) | 0 | output complete in every job |
| mcp_codon_align | MCP tool, 100,000 x 400 | 158.6 | server + engine | 2.8 | 121 | ok | all_ok, 100,000 OK |
| mcp_codon_align | MCP tool, 1,000,000 x 400, fasta / clustal + save_as folder | 1,587.8 | server + engine | 33.2 / 32.5 | 803 | ok | all_ok; FINDING 18 (run folders kept) |
| mcp_codon_align | MCP tool, 100,000 x 400 wrong code, default / strict | 164.7 | server + engine | 15.6 / 15.8 | 925 | ok | needs_attention, suggested_codontable 2 |
| mcp_explain_report | MCP tool on the 134 MB report of the 1,000,000 run | 133.7 | server | 2.0 | 802 | ok | plain-English summary |
| mcp_errors | MCP tool: missing file; save_as onto an existing file | | server | 0.01 / 0.34 | | tool error | clear messages, nothing overwritten |
| gui_100k | Gradio GUI (Windows), 100,000 x 400 by file path, code 5 then code 1 | 158.6 | app | about 15 / 6 | | ok | page responsive during the run; preview capped at 200 x 400 with the rest summarised; FINDING 7 at start-up |
| gui_extra_install | `pip install -e "apps/coal-app[gui]"` on Python 3.14 | | pip | | | 1 | FINDING 8 (pydantic-core has no 3.14 wheel) |
| pal2nal_alone | PAL2NAL v14 on 100,000 x 400 wrong code / -nomismatch / unrelated DNA | 164.7 | PAL2NAL v14 | 803 / 697 / 1,136 | 1,989 / 1,204 / 10,484 (maxrss) | 0 | output = truth (= COAL) |
| baseline_unittest | `python -m unittest discover -s tests` before any change | | both | 215 | | 1 | 84 errors, 6 failures, all WSL 1 `forkserver` (FINDING 19); all pass with `fork` |

## Findings

Severity order: crashes, lost output and wrong exit status first, then unbounded memory or super-linear time, then
rough edges. No run wrote a wrong alignment: every alignment that was written was byte-identical to the simulated
truth and to PAL2NAL v14 where it could run, and the two engines never differed in output (finding 1 is a run that
exits 0 with an empty or another run's file; the pieces that were written were right). Two of the findings (2 and 5)
are fixed in this PR; the rest are reported.

### Crashes, lost output and wrong exit status

**1. Two runs writing the same `-o` file at the same time: one can exit 0 with an empty file or another run's
output.** Every `-o` run writes its rows to the fixed names `<out>.coal-partial`, `<out>.coal-checkpoint.json` and
`<out>.coal-entries.jsonl` next to the output (`src/coal/checkpoint.py:56-60`, `cpp/src/main.cpp:648`) so that
`--resume` can find them, and renames the partial file into place at the end. Two runs with the same `-o` therefore
share one partial file: each truncates and writes it, the first to finish renames it away, and the others end with
`coal: error: cannot write same.fa: no such file or directory` (misleading: the folder is fine), or, worse, finish
"successfully" after the file they were writing was already renamed by the other run. Four runs of the 100,000 x 400
input with different options into one file, three trials per engine (`concurrent_same_o` rows):

| engine, trial | exit codes fasta / paml / clustal / -nogap | what `same.fa` held afterwards |
|---|---|---|
| coal-cpp 1 | 1 / **0** / 1 / **0** | the `-nogap` output: the paml run exited 0 but its output is gone |
| coal-cpp 2 | 1 / **0** / 1 / **0** | the `-nogap` output (same) |
| coal-cpp 3 | **0** / **0** / 1 / **0** | the `-nogap` output: two runs exited 0 with no output |
| Python 1 | 1 / **0** / 1 / 1 | **an empty file** (0 bytes) after an exit 0 |
| Python 2 | 0 / 1 / 1 / 1 | the fasta output (correct) |
| Python 3 | 1 / 1 / 1 / 1 | no file (correct: every run failed) |

With 1,000 identical processes writing one file (`concurrent_1000_same`), 536 of 1,000 (coal-cpp) and 495 of 1,000
(Python, 400 at a time) failed with the same "no such file or directory"; the surviving file was correct only because
every process wrote the same bytes. 1,000 processes writing 1,000 different files in one folder were fine (no temp
file collisions, no leftovers). What should happen: the pre-rename file should have a per-run name (as the
`.coal_<random>.tmp` files used without `--resume` do) or be created with `O_EXCL` and refused with "another COAL run
is writing `same.fa`"; a run must never exit 0 unless the bytes it wrote are in the output file. Smallest repro:

```
coal pep.fasta nuc.fasta -o same.fa -output paml & coal pep.fasta nuc.fasta -o same.fa -nogap & wait; echo $?
```

(any input that takes more than about a second; here `tests/tools/stress_gen.py d --nseq 100000 --ncodons 400 --gaps`).

**2. Python engine: more than 65,535 DNA files is an `OverflowError` traceback.** `DnaIndex` kept the file number
of every record in an `array("H")` (`src/coal/seqio.py:282`); record 1 of file 65,536 raised
`OverflowError: unsigned short is greater than maximum` in `_add`, a 26-line traceback with no `coal: error:` line
(`files_70000_ulimit_65536`, 70,000 one-record files, `ulimit -n 65536`). coal-cpp uses an `int` and failed later
and cleanly on the open-file limit (finding 3). **Fixed in this PR** (one line: `array("I")`, 2 more bytes per
sequence); with the fix the Python engine now fails like coal-cpp, with `coal: error: cannot open 'nuc_65518.fasta':
Too many open files`. Repro before the fix:

```
python tests/tools/stress_gen.py d --nseq 70000 --ncodons 50 --files 70000
cd d && ulimit -n 65536 && COAL_ENGINE=python coal pep.fasta nuc_*.fasta -o out.fa
```

**3. Both engines keep every DNA file open; the Python engine needs 129 KB per file and fails with 1,000 files
when it has workers.** `DnaIndex.get` opens each DNA file once and keeps the handle (`seqio.py:341-343`;
`seqio.cpp:479` does the same with an `ifstream`). Python 3.14 raised `io.DEFAULT_BUFFER_SIZE` to 128 KB, so each
handle costs 129 KB (measured with tracemalloc): 1,000 files = 130 MB, 60,000 files = 7.7 GB (`files_60000`: Python
8,083 MB peak against 663 MB for coal-cpp, for a 15.7 MB input). On Python 3.13 and older the buffer is 8 KB (0.5 GB
at 60,000 files). Memory therefore grows with the number of input files, not only with the number of sequences as
`docs/usage.md` says. The open-file limit is hit as well: with the default `ulimit -n 1024`, 1,000 DNA files work in
coal-cpp but fail in the Python engine as soon as it starts workers (the pool's pipes take the last descriptors):
`coal: error: cannot write out.fa: too many open files`, preceded by interpreter noise (`Exception ignored in atexit
callback ... AssertionError: can only join a child process`, and with clustal output `Exception ignored in atfork
callback ... SystemExit: 143`). `-j 1` works, 5,000 files fail in both engines with the clean
`coal: error: cannot open 'nuc_01018.fasta': Too many open files`. What should happen: close a file when the next one
is opened (records are fetched in protein order, which is file order for round-robin and sorted inputs), or an LRU
of a few hundred handles, and a smaller buffer (`buffering=65536` or the raw file); and the limit (about `ulimit -n`
minus 10 files) belongs in `docs/usage.md` until then. Repro: `stress_gen.py d --nseq 100000 --ncodons 400 --gaps
--files 1000`, then `COAL_ENGINE=python coal d/pep.fasta d/nuc_*.fasta -o out.fa` on a machine with 16 cores.

**4. A DNA FASTA whose line endings switch from LF to lone CR after the first 64 kB is refused with a wrong
message; PAL2NAL reads it.** Both readers decide the line-end style from the first 64 kB
(`seqio.py:250-259 _has_lone_cr`, `seqio.cpp:349-365`). A record that begins after that with classic-Mac `\r` ends
is not seen as a record: its `>` line is appended to the previous sequence's data (and loses its header) or, when
the first `\r` line is itself a header, everything after it becomes that header's description. The error is then
`coal: error: number of sequences differs (protein: 2000, DNA: 301): ... 1699 proteins have no DNA record`, which
sends the user looking for missing records in a file that has them all. PAL2NAL v14 reads the same file and its
output equals the truth (`fmt_mixed_lf_then_cr`; LF then CRLF is fine: `fmt_mixed_lf_then_crlf`). Smallest repro
(`tests/test_stress_findings.py::MixedLineEndings`, an expected failure for both engines): a two-gene DNA file whose
first record is 65,600 bytes with LF and whose second record is `>b\rATG...\r`: `coal: error: b: protein and DNA are
inconsistent (no match found)`; the same file with a 1 kB first record works. Such files come from concatenating
exports of different tools; what should happen is what the protein reader already does (Python's universal
newlines): treat a lone `\r` as a line end wherever it is, or at least say "line 123: a lone CR after LF line ends".

**5. `nohup coal ...` died when the terminal closed (exit 129), in both engines.** `install_signal_handlers`
(`src/coal/cli.py:159-161`, `cpp/src/sysio.cpp:73-75`) installed COAL's SIGHUP handler whether or not the signal was
already ignored, which is exactly what `nohup` (and `setsid`) arrange; a plain Python process under `nohup` ignores
the hang-up and finishes, COAL exited 129 with no output (`nohup_sighup`, 1,000,000 x 400, `kill -HUP` 4 s in).
Long runs started with `nohup` from an SSH session are the common case for an aligner of this size. **Fixed in this
PR** (the handler is installed only when SIGHUP is not ignored, the usual idiom; a plain SIGHUP still cleans up and
exits 129, and `tests/test_stress_findings.py::SighupUnderNohup` checks both engines by holding them on standard
input while the ignored signal arrives).

**6. Python engine out of memory in the main process: a raw `MemoryError` traceback.** Inside a Windows Job Object
limited to 256 MB, the 1,000,000 x 400 input made the Python engine die with a 15-line traceback ending in
`MemoryError` (`win_memlimit_256MB_python`, with and without `-j 1`; the engine needs about 314 MB for this input in
one process). A worker running out of memory is handled (`coal: error: a worker process died (out of memory?)`,
`win_memlimit_512MB_python`), the main process is not. coal-cpp finished the same input in 228 MB inside the 256 MB
job. What should happen: `coal: error: out of memory ...` with exit 1, like the other failures (the output file was
not left behind, which is right).

**7. `python -m coal_app` on Windows crashes at start-up when stdout is not UTF-8.** Gradio's MCP banner prints an
emoji; with a cp1252 console or a pipe, `launch()` dies with `UnicodeEncodeError: 'charmap' codec can't encode
character '\U0001f528'` before the server answers (gradio 5.49.1, Python 3.13.5, `coal_app/ui.py:768`). Setting
`PYTHONUTF8=1` (or `sys.stdout.reconfigure(encoding="utf-8", errors="replace")` before `launch`, or `--no-mcp`)
avoids it. The frozen Windows app was not tested.

**8. The `gui` extra cannot be installed on Python 3.14.** `gradio[mcp]==5.49.1` pins pydantic-core 2.33, which
has no Python 3.14 wheel and whose PyO3 refuses to build for 3.14 (`error: the configured Python interpreter version
(3.14) is newer than PyO3's maximum supported version (3.13)`), so `pip install -e "apps/coal-app[gui]"` fails on
this machine's Ubuntu (Python 3.14.4). The `mcp` extra (mcp 2.2.0) installs and works. The GUI was tested with
Windows Python 3.13 instead.

### Unbounded memory and super-linear time

**9. Memory grows with the number of per-codon warnings, not with the number of sequences: a million genes of
wrong DNA cannot be run on a 32 GB machine, in either engine.** `docs/usage.md` says a sequence with warnings keeps
its report entry in memory, about 0.6 KB (coal-cpp) or 1.4 KB (Python). That is per *sequence*; the entry holds
every per-codon warning, and a gene whose DNA is unrelated to its protein has one per codon (400 here):

| input (x 400 codons) | PAL2NAL v14 | coal-cpp | Python engine |
|---|---:|---:|---:|
| 100,000 genes, wrong genetic code (about 13 mismatches each) | 803 s, 2.0 GB | 14 s, **353 MB** | 22 s, **996 MB** |
| 1,000,000 genes, wrong genetic code | not run | 218 s, **3.4 GB** | 350 s, **6.4 GB** |
| 100,000 genes of unrelated DNA (every codon mismatches) | 1,136 s, 10.5 GB | 451 s, **5.5 GB** | 537 s, **14.4 GB** |
| 1,000,000 genes of unrelated DNA | not run | **killed at 16.4 GB** after 321,000 genes (19 min) | **killed at 16.6 GB** after 78,000 genes (2.7 min) |

(the clean 100,000-gene run: 35 MB / 444 MB; the last row was stopped by a watchdog at 16 GB to keep the machine
usable; the Python engine was about 150 KB per gene, coal-cpp 55 KB; PAL2NAL's numbers are its own stand-alone runs.)
The output of every finished run was byte-identical to PAL2NAL's and to the simulated truth, and `--strict
--keep-going` on the unrelated input (every gene refused, exit 1) needed only 76 / 518 MB, so it is the kept
warnings. `--max-warnings 0 --report` on the million wrong-code genes wrote 13.4 million warning lines (520 s / 650 s)
at the same memory. What should happen: keep per sequence only what the summary and the report need (counts, the
first few warnings) and stream the rest, as the report writer already streams `--report`; the documented per-sequence
figures are then true again. Repro: `stress_gen.py d --nseq 100000 --ncodons 400 --mismatch unrelated`, then
`coal d/pep.fasta d/nuc.fasta -o out.fa` (`--table 2` + `-codontable 1` for the wrong-code variant).

**10. coal-cpp's exact-fit search is quadratic in the gene length on low-complexity DNA.** `find()` in
`cpp/src/core.cpp:187-205` tries every start position and verifies each one codon by codon; when a long stretch of
the DNA almost fits the protein, each start costs up to the whole gene before failing. A poly-A gene placed after an
equally long poly-A stretch (`--repeat poly --repeat-prefix 1`; the Python engine finds the same, correct answer
through its translated-frame search, `core.py:436-480`):

| residues | coal-cpp | Python engine |
|---:|---:|---:|
| 10,000 | 0.42 s | 0.11 s |
| 20,000 | 1.5 s | 0.15 s |
| 50,000 | 8.5 s | 0.24 s |
| 200,000 | 138 s | 0.73 s |
| 1,000,000 | **not finished in 30 min** (killed) | 4.0 s |

A (GPP)n gene behind a (GPP)n stretch: 3.3 s against 0.22 s at 50,000 residues. Real inputs this long and this
repetitive are rare, but the same scan runs for every exact-fit search, so any long gene whose DNA has a near-copy of
its start pays for it. What should happen: find the first run of ordinary residues in the translated frames and
verify candidates with one memory compare, as the Python engine does (`_plain_run` / `frames[st % 3].startswith`);
`tests/test_stress_findings.py::QuadraticExactFitScan` (COAL_SLOW) measures the ratio. Repro: `stress_gen.py d --nseq
1 --ncodons 200000 --repeat poly --repeat-prefix 1 --no-start-met`, then `coal-cpp d/pep.fasta d/nuc.fasta`.

**11. The Python engine's fallback on one very long gene is 15-20x slower than with `-q`, because of the
genetic-code hint.** One wrong residue at the end of a 1,000,000-residue gene: 30.2 s and 694 MB; the same run with
`-q` (which skips the hint) 1.7 s and 266 MB; coal-cpp 5.7 s / 175 MB (0.5 s with `-q`). At 200,000 residues: 5.65 s
against 0.40 s (coal-cpp 1.25 s against 0.16 s). `-nomismatch`, which prints no mismatch warnings, still pays it
(34.5 s). `docs/benchmarks.md` already names the hint as the cost of the fallback on long genes; the numbers here
say it is 95 % of the run and 2.6x the memory, so re-testing the gene under 17 codes could stop after the first few
thousand codons, or be skipped when `-q`/`-nomismatch` means nobody reads it.

**12. Many worker processes buy no time after four, but cost about 45 MB each.** 1,000,000 x 400 with the Python
engine: `-j 1` 85 s / 347 MB, `-j 2` 46 s / 473 MB, `-j 4` 32 s / 532 MB, `-j 8` 32 s / 679 MB, `-j 16` (default here)
34 s / 1,077 MB, `-j 32` 34 s / 1,875 MB (PSS). Output byte-identical in all. The default is "all cores", so a
32-core node spends 1.9 GB for the speed of 4 processes; `docs/benchmarks.md` says four is the ceiling, so the default
could stop there. coal-cpp's threads gave nothing on this machine: 14.5-18 s for `-j 1` to `-j 32` (the docs promise
about 1.2x from two threads on; here the run is bound by reading 1.2 GB and writing 1.2 GB on WSL 1's file system).

### Rough edges

**13. The Python engine takes 109 s for the codon format of 1,000,000 x 400 (3.7 GB of output), 3.3x its fasta
time, while coal-cpp takes 25.5 s (1.5x).** Clustal: 54 s against 23 s. Fine, but the interleaved writers are where
the Python engine's time goes at this size (`many_1M_codon`, `many_1M_clustal`).

**14. Error messages that list every input file.** With thousands of DNA files the count and "no sequences" errors
start with every path joined by ", " (`seqio.py:206`, `pipeline.py:766`): a 1,022-file message is a 16 kB line. Name
the first few and the count.

**15. The first checkpoint comes late.** Checkpoints are written "about every 2 s of processing", but processing
starts after the index and pairing passes: a 1,000,000 x 400 Python run interrupted 12 s into its 35 s had no
checkpoint yet (the C++ engine had one after 6 s of 13 s). Everything else about `--resume` held: 36 kill/resume
rounds (SIGTERM and SIGKILL at random points, fasta, clustal, paml, `-nogap --report`, `--keep-going`, both engines)
gave byte-identical output with no leftovers; truncated, non-JSON, out-of-range (`completed`, `sink_bytes`) and
partial-file-deleted checkpoints were all refused with the documented notice and a correct normal run; the other
engine's checkpoint was refused with the documented error (exit 1).

**16. A SIGKILL during the second pass of clustal/codon output leaves the fixed-width scratch copy in `TMPDIR`.**
One of the kill-and-resume rounds (Python engine, 1,000,000 x 400, clustal, SIGKILL during the final pass) left
`/tmp/.coal_rows_v91_2rfn`, 743 MB, which nothing removes afterwards: it is not one of the files `--resume` knows
about (`<out>.coal-partial`, `.coal-stage`), and SIGKILL cannot run the cleanup. `--resume` then redoes the pass
correctly, so the output is right, but a few such kills fill `/tmp`. A stale-file sweep at start-up (files named
`.coal_rows_*` older than a day in the temp folder), or writing the scratch copy next to the output under a name
the checkpoint records, would cover it.

**17. `--progress` lines differ between engines** (they are timing-based: `coal: progress: 4% ... about 39 s left`
against `18% ... about 9 s left`), so stderr is not byte-identical with `--progress`; everything else was.

**18. The MCP tool keeps every run's output in its run folder for 24 h.** Three `codon_align` calls on the
1,000,000 x 400 input left 4.4 GB in `/tmp/coal_app_runs`; `save_as` copies the alignment rather than moving it.
`explain_report` on the 134 MB report of that run needs 800 MB. Otherwise the tools were fine: 100,000 sequences in
2.8 s, 1,000,000 in 32-33 s (fasta and clustal, `save_as` into a folder), the wrong-code input in 16 s with the right
hint (`suggested_codontable: 2`), missing files and `save_as` onto an existing file refused with clear messages.

**19. Python 3.14's default `forkserver` start method cannot start on WSL 1, and COAL shows a raw traceback.**
Every multi-process run died with `EOFError` from `multiprocessing/forkserver.py` (84 errors in the baseline suite).
This is the platform's fault (fork and spawn work), but the engine could catch the pool start-up failure and say
`coal: error: could not start worker processes (...); run with -j 1`.

**20. The 100 MB single line is fine but not free:** one record of 400 codons plus a 100 MB 3' UTR on one line
takes 2.9 s / 543 MB in the Python engine and 1.0 s / 306 MB in coal-cpp (three copies of the line).

## What I did not get to

- **A full disk.** WSL 1 does not enforce `ulimit -f`, lets a normal user's writes to `/dev/full` succeed (and
  root's too), has no size-limited tmpfs and no loop devices, so neither "output to a full disk" nor "`TMPDIR` full"
  could be produced; on Windows a small volume needs administrator rights. Writing through a full device was
  therefore not tested at all. The read-only folder, unwritable/missing `TMPDIR` and the output-is-a-folder /
  input / symlink-loop cases were.
- **`ulimit -v` on Linux.** Not enforced by WSL 1; the 1 GB input under 256 MB was run on Windows inside a Job
  Object instead (Python engine: `MemoryError` traceback at 256 MB, clean worker-death error at 512 MB with
  workers, success in 314 MB with `-j 1`; coal-cpp: 228 MB).
- **1,000 concurrent Python-engine processes.** 1,000 `coal-cpp` processes ran at once; the Python engine was run
  400 at a time (1,000 Python processes need about 13 GB of RAM on this 32 GB machine with Windows holding 20 GB).
- **PAL2NAL on the mismatched 100,000-gene inputs** ran stand-alone with a 30-minute cap (see the table); PAL2NAL
  on 1,000,000 sequences (62 minutes in `docs/benchmarks.md`) was not run.
- **The genetic-code probes of the app, the hosted Space, the frozen Windows app and the GUI's file upload path**
  (only the "Use file paths" boxes were used), and the GUI under Linux (the `gui` extra does not install on Python
  3.14).
- **UNC paths on Windows** (`\\wsl$\...` was not reachable from the Windows Python) and paths longer than 260
  characters.
- **Sanitizer builds** (`-fsanitize=address,undefined`) of coal-cpp on the large inputs: the nightly fuzz covers
  small random inputs; a sanitizer run of the 1,000,000-sequence and 1,000,000-residue cases was not done.
- **numpy acceleration** of the Python engine (`core.py:_numpy`): numpy was not installed, so every Python-engine
  number is the pure-Python path.
- **Timings on a real Linux kernel.** Everything ran on WSL 1 (slow file opens, no `forkserver`); absolute times are
  comparable within this report, not with `docs/benchmarks.md`.


## Re-test on v0.5.0

v0.5.0 (tag `v0.5.0`, 22 commits after the tested `c6f5728`) changes only version strings, one cosmetic vector
initialiser in `cpp/src/report.cpp`, packaging, the website and docs: `src/coal`, `cpp/src/core.cpp`,
`cpp/src/seqio.cpp`, `cpp/src/sysio.cpp` and the app are byte-identical apart from the version. Every finding was
rerun against an untouched v0.5.0 build (`~/work/COAL050`, its own venv and `coal-cpp`; `group050.jsonl`) and
reproduced unchanged:

| finding | v0.4.2 (`c6f5728`) | v0.5.0 |
|---|---|---|
| 1 concurrent runs into one `-o` | exit 0 with another run's / an empty output | same (cpp: `-nogap` run's output survives a paml run that exited 0; python: 2 of 2 trials) |
| 2 > 65,535 DNA files | `OverflowError` traceback (python) | same |
| 3 1,000 DNA files, 16 workers | python "too many open files", cpp fine | same; 60,000 files 8.1 GB |
| 4 LF then lone CR after 64 kB | both engines: count error | same |
| 5 `nohup` + SIGHUP | exit 129, no output | same |
| 8 gui extra on 3.14 | pydantic-core build fails | same pins |
| 9 memory per warning | 20,000 unrelated genes: 3.1 GB (python) / 1.15 GB (cpp) | same (3.14 / 1.15 GB) |
| 10 quadratic exact-fit scan | 50,000 residues: cpp 8.5 s | 9.2 s |
| 11 hint cost | 1M-residue gene: python 30 s / 694 MB, `-q` 1.7 s | 35.8 s / 691 MB, `-q` 1.6 s |
| clean inputs | 1M x 400: cpp 16.9 s / 201 MB, python 34 s / 1.07 GB | 16.9 s / 201 MB, 37.7 s / 1.07 GB |
| `--resume` | byte-identical | byte-identical (2 rounds per engine) |

## Fixes (branch `stress-fixes-0.5.0`, on v0.5.0)

Each fix is in both engines unless said, with a test in `tests/test_stress_findings.py` (17 tests, 1.8 s; the
quadratic one under `COAL_SLOW`), and every differential suite (`test_differential`, `test_cpp_vs_python`,
`test_cpp_vs_reference`, `test_realdata`, `test_engine_fuzz`) plus 2,300 extra engine-fuzz cases (seeds 7, 11, 23,
31) agree between the engines after them. The full suite shows only the two WSL 1 artefacts (`/proc` permission
wording, the pool test under `fork`).

| finding | fix | commit |
|---|---|---|
| 1 | `<out>.coal-lock` (process id) next to the output while a run writes it; a second run on the same file is refused with `coal: error: another COAL run (process N) is writing out.fa: wait for it to finish, or write to another file (if that process is gone, remove out.fa.coal-lock)`, exit 1; a lock of a dead process is replaced; removed by the normal end, an error and the signal handlers | `a7b40ae` |
| 2 | `DnaIndex.fno` is an unsigned int (python) / `uint32_t` (cpp, which wrapped silently) | `e2e7699`, `a7b40ae` |
| 3 | an open-file budget: the soft limit raised towards the hard limit when the files need it, at most limit minus headroom files open, least recently used closed first; python read buffers sized to about 64 MB in all. 600 files round robin with a hard limit of 160 descriptors work in both engines | `a7b40ae` |
| 4 | `\n`, `\r\n` and a lone `\r` end a line wherever they are in a DNA file (as the protein reader) | `6416b17` |
| 5 | an already ignored SIGHUP stays ignored | `19ab856` |
| 6 | `coal: error: out of memory: ...` (exit 1, report code `out_of_memory`) for a `MemoryError` / `std::bad_alloc` in the main process | `6416b17` |
| 7 | the GUI's console output uses `errors="replace"` on a non-UTF-8 console (`tolerant_console`); the MCP stdio server is untouched | `6416b17` |
| 9 | partial: with `-q` and no `--report`/`--json` the run report keeps no per-sequence entries (20,000 unrelated genes: cpp 103 MB instead of 1,149 MB; python 1.16 GB instead of 3.15 GB, the rest being warnings in flight between workers). The default mode still keeps every message for the summary and the hint; `docs/usage.md` now gives the per-gene cost and the follow-up (keep only what the summary and hint need, stream the rest to the entries file) | `0eb8791` |
| 10 | coal-cpp `find()`: a pattern of ordinary residues is a KMP string search in the translated frames; mixed patterns search their first run of residues and probe the last residue before the full check. 200,000-residue repeat: 138 s to 0.16 s; 1,000,000: not finished in 30 min to 0.8 s | `a7b40ae` |
| 11 | genes of more than 100,000 residues are left out of the genetic-code hint's sample: the 1M-residue gene with one wrong codon takes 1.6 s (python) / 0.5 s (cpp) instead of 30 s / 5.7 s | `99bd164` |
| 14 | messages name at most three input files and the count | `99bd164` |

Not fixed, by choice: 12 (the default worker count; a documented ceiling), 13, 15, 16, 17, 18, 19 (environment),
20, and the default-mode memory of 9 (a redesign of the run report's per-sequence store, both engines).
