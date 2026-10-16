# Benchmarks

```
python benchmarks/run_bench.py --quick          # smoke test, seconds
python benchmarks/run_bench.py                  # full grid (PAL2NAL may take minutes on the largest sizes)
python benchmarks/run_bench.py --sizes 1000x500 100x10000 --repeats 5
python benchmarks/run_bench.py --coal-jobs 2 4 0                     # also time COAL with -j 2, -j 4, all cores
python benchmarks/run_bench.py --no-pal2nal --sizes 500000x500 --coal-jobs 0   # scale test, COAL only
```

- `gen_data.py` makes synthetic datasets that are correct by construction: N descendants of a random
  ancestral gene of L codons, ~10% codons mutated, ~5% deleted (gaps). Cached in `benchmarks/data/` (git-ignored).
- Each tool is run as a command-line process (times include Perl/Python start-up). Median of `--repeats` runs.
- Peak memory needs `pip install psutil`; otherwise those columns are blank.
- `outputs_identical` = COAL FASTA output byte-identical to PAL2NAL `-output fasta`.
- Results: `results/<name>.csv` + `results/<name>.txt` (platform, versions). Only compare numbers from the same machine.
- Large datasets take disk space in `benchmarks/data/` (500,000 × 500 is about 1 GB) and a few seconds to generate.
