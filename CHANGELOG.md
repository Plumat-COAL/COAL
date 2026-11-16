# Changelog

## Unreleased

### Fixed

- **The development tree is 0.6.1.dev0** (coalr 0.6.1.9000): `main` still said 0.6.0 after the changes below
  were merged, so a build of it reported, and `coal` accepted as, the published release. The six version places,
  CMake, coalr's DESCRIPTION, the vendored `coal_version.h`, `CITATION.cff` (no `date-released` in a development
  tree) and the website's version spans name the development version or the release it becomes, and the
  dev-version banner is back at the top of README.md, docs/install.md and docs/tutorials/README.md. A release
  version with entries under "Unreleased" now fails `tests/test_feature_combinations.py`, and
  `scripts/release_version.py` refuses its tag (also when coalr's DESCRIPTION has another version).
  (Final review of 2026-11-16, H1.)
- **coalr install lines name a release**: docs/r.md, README.md, coalr's README, tutorial 4 and
  packaging/orthologr/README.md install `coalr_0.6.1.tar.gz` of the v0.6.1 release (or `install_github()` at the
  tag `v0.6.1`); the unpinned `install_github()` is shown only as the development version, and until 0.6.1 is
  published a note gives the v0.6.0 tag. `tests/test_docs_site.py` keeps the lines at the release version.
  (Final review of 2026-11-16, H1.)
- `release.yml`: the dry-run version check is `scripts/release_version.py --dry-run` (tested) instead of the same
  inline Python in two jobs, and a release that already exists (a re-run, or one made on GitHub's release page)
  gets the prerelease flag of a new one (`gh release edit --prerelease=<true|false>`; it kept whatever it had).
  (Final review of 2026-11-16, L5.)
- The release checklist is one list (docs/development.md "Releasing COAL"; CONTRIBUTING.md points to it) and a
  release commit made by it passes the tests: `tests/test_user_files_and_lows.py` finds the post-0.6.0 entries
  under the first section, whatever it is called; the list names the report example's version strings, coalr's
  NEWS heading (`scripts/release_version.py` also refuses a tag whose NEWS.md does not start with
  `# coalr <version>`) and what a prerelease (`rc1`) sets: COAL's places, `CITATION.cff` and the site at
  `X.Y.Zrc1`, coalr left at its development version, which `tests/test_r_bindings.py` and coalr's own test accept
  for a prerelease. Checked by making the 0.6.1 and 0.6.1rc1 commits on a copy and running both test suites.
  The release notes list `coalr_<version>.tar.gz` / `coalr_src.tar.gz` (not for a prerelease), and the `files`
  job fails a release (not a prerelease) without coalr. (Audit of #100.)
- The Bioconda recipe's notes no longer call 0.6.0 the next release. `ClustersMsa.test_mafft` no longer skips
  with "needs coal-cpp and Clustal Omega" where MAFFT is installed and Clustal Omega is not. (Audit of #100.)
- `tests/test_orthologr_standin.py` gives its subprocesses an absolute `PYTHONPATH` (`tests/subprocess_env.py`):
  with `PYTHONPATH=src` its run in `examples/clean` failed with "No module named coal_compat". (L4.)
- `coal pipeline`: a `.fasta` file of your own in `DIR/codon_alignments/` or `DIR/protein_alignments/` is no longer
  deleted by the next run into `DIR`. A cluster run records the names of the per-cluster files it wrote
  (`codon_alignments_files` / `protein_alignments_files` in `work/outputs.json`) and the next run removes those
  only; a run with `--keep-alignments` / `--keep-protein-alignments` into a folder that holds a file of yours is
  refused with exit 1 before anything changes, since one of its files could have that name. A folder recorded by
  0.6.0 (no list) is cleared as before. (Final review of 2026-11-15, L3.)
- `coal pipeline`: a `DIR/work/` folder that no COAL run made (not empty, no `work/outputs.json`, no COAL
  `report.json` of an older version) is no longer written into: the run is refused with exit 1 and a message
  before anything changes; `--overwrite` writes into it (files there with COAL's own names are replaced; nothing
  else is touched). A `work` that is a file is
  refused, also with `--overwrite`. (L2.)
- `coal pipeline --report FILE`: an existing `FILE` (outside `DIR` or in it) that is not a COAL report is no longer
  replaced without a word: exit 1 before anything changes, unless `--overwrite`. A COAL report there (`"schema":
  "coal.<name>/<n>"`, as `coal pipeline` and `coal --report` write) is replaced as before. (L1.)
- `coal pipeline --digits N` (and `digits=N` from Python) now also sets the digits of `dnds.kaks`
  (`--kaks-format`), which always had 6; without `--digits` it keeps KaKs_Calculator's 6, so a run without it
  writes the same files as before. `coal-cpp codon-dnds` has `--kaks-digits N` for `--kaks-out` (default 6).
  `digits` outside 1 to 17 is refused by the Python API too (it reached coal-cpp). (L6.)
- `coal-cpp kaks` / `dnds`: the refusal of a gzip-compressed input names the commands that decompress it for
  this engine, the Python package's `coal kaks` and `coal dnds` (it named `coal` and `coal pipeline`). (L5.)
- Windows, MAFFT (AUDIT-91 notes): a `mafft` on `PATH` that is a batch file of your own (not MAFFT's package, so
  it would run through `cmd.exe`) gives way to the pinned `mafft-win\mafft.bat` package in COAL's cache when that
  is installed and unchanged, with a note once; `COAL_MAFFT` still picks any batch file. A batch file of your own
  is no longer refused because of the folder COAL was started from: it runs in the run's folder, and now with an
  environment without Windows' hidden per-drive current folders, so that folder never reaches it (its own folder,
  the run's folder and `TMP` / `TEMP` are checked as before). MAFFT's package in a folder whose path has a `;`
  (which would split `PATH`) is refused with a message; non-ASCII letters in its folder are tested, mocked and in
  CI's real download (now under a cache folder with such letters). When that real-download test cannot reach
  mafft.cbrc.jp, its skip is shown as a warning on the CI run and in the job summary, with the reason.
- `coal pipeline`, follow-ups of the audit of #96: a cluster run's list of alignment files is empty from its start
  and written again when the engine ends with an error or an interrupt too, so such a run no longer leaves a
  folder whose every `.fasta` counts as COAL's (while the engine writes there is no list, so a run killed then
  is cleared and resumed as before); an error writing the list is a warning, not a traceback. A `.FASTA` (any
  case but the engine's `.fasta`) in such a folder counts as yours. A `work/` with only an empty `logs/`, or
  whose `work/outputs.json` begins as COAL writes it but cannot be read, is no longer refused as someone else's
  `work/` (output-named files in `DIR` are still refused when no record or report names them). `--report` naming
  a folder, a device or a named pipe is refused before the run (a folder failed at its end, a pipe blocked it). report.json
  `parameters` has `kaks_digits` (the digits of `dnds.kaks`). docs/pipeline.md says what a deleted `work/` means
  for the alignment folders; docs/alignment.md no longer lists the current folder among a batch file's checks.
- Texts: `CITATION.cff` has `date-released: 2026-11-15` for 0.6.0; `docs/python-api.md` no longer says its
  examples were run with 0.6.0.dev0; coalr's NEWS no longer lists a development-version note for 0.6.0;
  `scripts/release_notes.py` no longer gives a wrong size for the 0.6.0 section.
- **GY: a ridge start from each of the first start's two runs** (review of main after #97, M1). The ridge start of a
  diverged pair (t at its bound) takes kappa and omega t from where the first start stopped. Since the first start
  can be run a second time with longer steps (below), that was the second run's stop, and from there the ridge
  start could reach a lower part of the ridge than from the first run's stop, which 0.6.0 used: `GCTGATCCC` against
  `GATGCCGCG` under code 4 had lnL -12.5851 (kappa 17) where 0.6.0 had -12.5748 (kappa 1e-6). A pair whose first
  start was replaced by its second run now gets a ridge start from each stop. That pair is back at -12.5748, and
  `GTTAAGTAC` against `GGTAATGAG` under code 6, `saturated` with lnL -15.8668 in 0.6.0 and since, is `ok` with lnL
  -15.3300. Nothing else changes: at 17 digits no GY row of 4,388 yeast, 2,862 *Salmonella* / *Shigella* and 3,142
  *E. coli* / *Salmonella* pairs differs, and of the 25,200 random pairs four do, each with a higher likelihood (by
  0.010 to 0.24; two go from `saturated` to `ok`, none to `not_converged`, of which there are still none). Against
  0.6.0 no pair of the 25,200 has a likelihood lower by more than 1e-4 any more.
- GY, inside: the optimiser's stopping rule near a bound, the lengthening of a step and the choice among the starts
  are functions of their own (`probe_zone_ends`, `lengthen_step`, `choose_start`, `apply_saturation`), and the
  choice among the starts has a test with made-up starts (`tests/cpp/kaks_gy_choose_start.cpp`): the rule that an
  unconverged start displaces a usable one only when it is more than 1e-4 higher, which no pair found so far
  reaches, is now held in place. The same bytes at 17 digits as before this change on the fixtures, the three sets
  above and the 25,200 random pairs.
- **GY: a coarse scan of kappa and omega for small pairs and for pairs next to the plateau of saturation** (audits
  of #101 and #103). Fixed starts are not a global search: `GTGGACAAC` against `TCCGGGAAG` under code 11 was
  `saturated` with lnL -17.3590 (every start ended on the ridge of saturation, where the likelihood barely depends
  on kappa and omega), and its maximum is -16.5761 at t 7.96, kappa 1e-6 (`ok` now). After the ordinary starts and
  their continuation, which are as before, a pair of 40 codons or fewer, a pair whose result did not converge, and a
  pair whose likelihood is within 0.05 lnL per codon of the plateau of saturation (sum log pi_a pi_b, where kappa
  and omega have no gradient left) has the likelihood, maximised over t, evaluated at 36 points of a grid in (kappa,
  omega), and is started from the two highest (a pair of more than 40 codons only from a point higher than its
  result). The new choice is taken only if the likelihood reported for it, after the flat-likelihood rule, is higher
  than without the scan. More than 1e-3 below the global maximum of an independent likelihood, without and with the
  scan (this change's samples and its audit's, each with its own search): random pairs of 1 to 40 codons 67 and 4 of
  600, 75 and 8 of 642; of 60 to 200 codons 14 and 0 of 150, 11 and 2 of 150; of 300 to 800 codons 5 and 1 of 60, 9
  and 1 of 60. Most of the remaining maxima have omega or kappa at a bound, not all. No pair of more than 40 codons
  in the real sets below is within 0.18 per codon of the plateau (the smallest: 0.182, yeast), so the scan leaves
  them alone and GY's CPU time there is the same as before (two intermediate versions cost 8 to 12 %). Against main
  at 17 digits: of 37,800 random pairs (seeds 1, 2, 97, 98, 11, 4242; nine codes) 6,837 have a higher likelihood (up
  to 3.4), none is lower, none is `not_converged`; 1,099 go from `saturated` to `ok`, 363 to `warning`, 281 from
  `ok` and 138 from `warning` to `saturated`. No GY row of the 4,388 yeast, 2,862 *Salmonella* / *Shigella* and
  3,142 *E. coli* / *Salmonella* pairs changes except the identical pairs (next). docs/kaks.md no longer says that
  COAL's likelihood is the higher one against codeml as if guaranteed: it is on every pair compared, and the search
  is local.
- **GY: the likelihood of identical sequences** is the sum of log pi over the codons (the likelihood at t = 0); it
  was printed twice that (for 4 codons about 10 too low). Only the lnL column of `--ml-columns` changes: 665 of
  the 37,800 random pairs and 5 of the real pairs above.
- GY, inside: the continuation of an unconverged start is a function (`continue_starts`) with a test through a
  second hook (`tests/cpp/kaks_gy_choose_start.cpp`: at most two continuations, the other unconverged starts
  dropped, one that ends lower replaced by its start), and the scan's gate (`scan_wanted`) has one too, with cases
  either side of 40 codons and of 0.05 per codon; the hooks are compiled into every build of the library
  (COAL_TEST_HOOKS applies only to coal-cpp, and no CI job builds a hooked library; the comment says so); the
  second ridge start is not made where it would repeat the first; `ridge_omega`, `first_run_stop_*` names. The
  same bytes at 17 digits as main on the fixtures, the three real sets and the 37,800 random pairs.
- **GY: a start of the optimiser that crawled to its step limit is run again with longer steps.** On a nearly flat
  likelihood (two random codons with every site different: `AACTGG` against `GGTACA`) the first steps have the
  gradient's own size, 1e-4 in log kappa and log omega, and no curvature for the BFGS update to scale them with, so
  the iteration used its 100 steps to gain 3e-7 log units and the pair was `not_converged`. Such a start (limit
  reached, not converged, t not at its bound) is now run a second time with these steps lengthened while the
  likelihood improves, and the second run is taken if it converges with a likelihood that is not lower. That pair is
  now `saturated` (t at its bound, lnL -13.5807 where the crawl stopped at -13.6046). Every start that converges or
  ends on the ridge of saturation is computed exactly as before: same bytes at 12 digits for the 20 codeml fixtures,
  128 simulated pairs, the test sets, and 4,388 yeast, 3,685 *E. coli* and 159 *Salmonella* / *Shigella* pairs; of the
  audit's 568 random and high-divergence pairs, 5 rows change (the pair above, and four saturated pairs whose third or
  fourth start had crawled: lnL equal to 1e-8, with the kept estimate moving on the flat likelihood: the Ka of
  `rand_L2_11` in the fourth digit, 5.3289 to 5.3304, and one dS in the fifth). "Not lower" means not lower by more
  than 1e-7, the tolerance within which two likelihoods count as the same maximum. The nine other methods are
  byte-identical.

- **GY: a converged maximum is no longer displaced by a start that ran out of steps** (audit of #95, M1). The four
  starts of a diverged pair were compared by likelihood alone, so a start that stopped at its step limit with a higher
  likelihood made the pair `not_converged` although three starts had converged (`TGGTTACGCTTT` against `TACAGCTTTCGC`:
  `saturated` in 0.6.0, `not_converged` after the change above). Now a start is usable if it converged or ended on the
  ridge of saturation (t at its bound, converged or not); the highest usable start wins (the earlier or the converged
  one within 1e-7); a start that is not usable is continued from where it stopped, twice at most, when its likelihood
  is higher than that of the best usable start, and only an unconverged point more than 1e-4 above the best usable
  start, or the absence of any usable start, gives `not_converged`. That pair is `saturated` with lnL -24.0821 (0.6.0:
  -24.1755). A start that never took a step the second run would lengthen is not run a second time (it would repeat
  itself step for step). On 25,200 random pairs of 1 to 40 codons under nine genetic codes: `not_converged` 103 in
  0.6.0, 47 after the change above, none now, and no pair goes from `ok` or `saturated` to `not_converged`.
- **GY: kappa or omega near a bound** (audits of #95, L5, and of #97, M1). Near 0 the likelihood is linear in kappa or
  omega, so the gradient in the logarithm fades before the bound: the iteration stopped about a factor 10 short of it,
  at a point that depended on the path (omega 1.1e-6 to 1.5e-5 instead of 1e-6, up to 1e-5 log units below the
  likelihood at the bound), and an iteration near or on the bound could not see a maximum further inside. When the
  iteration stops with a parameter below 1e-3 or above 99, both ends of that zone are now tried: the bound, taken if
  the likelihood there is not lower, and the inner edge (1e-3 or 99), taken if the likelihood there is higher; the
  iteration goes on from the better one. Later starts are seeded from where the first one stopped by itself, never
  from a bound it was moved to (a start on the bound stalls there: one random pair of 25,200 lost 6.6e-3 log units
  that way in the first version of this change). An iteration that stops outside the zones is as before: of 4,388
  yeast pairs no GY row changes at 12 digits; of 2,862 *Salmonella* / *Shigella* pairs 6 and of 3,142 *E. coli* /
  *Salmonella* pairs 8 (9 different pairs: omega or kappa now exactly 1e-6, lnL higher by 7e-7 to 6e-6, dS the same to
  six digits, no status change); one of the 20 codeml fixtures (`real0_8`: omega 1.143e-6 to 1e-6). On the 25,200
  random pairs, against the change above: no likelihood lower by more than 1e-4 (three flat pairs lower by 1.4e-7 to
  1.8e-7, the tolerance of the flat-likelihood rule, which now calls them `saturated`); 39 flat pairs of 2 to 4 codons
  go from `ok` or `warning` to `saturated` and one back. Against 0.6.0 one pair of the 25,200 had a lower likelihood
  (3 codons, saturated in both, 0.010 lower since the change above): `GCTGATCCC` / `GATGCCGCG`, put right by the
  ridge start from each run (above).
- `coal pipeline clusters --keep-protein-alignments` (or `--keep-alignments`): when the engine step ends while an
  alignment folder it writes is not there (`os.makedirs` failed: permission denied, disk full), the run now records
  an empty list of files for it in `work/outputs.json`. It recorded none, which means "every `.fasta` here is the
  run's", so a `protein_alignments/` you made afterwards had every `.fasta` in it removed by the next run into `DIR`;
  now they are yours (kept, and a run writing that folder is refused as for any file of yours there). The warning
  for a list that could not be written now says what follows: the next run takes every `.fasta` file there for this
  run's (it said the next run keeps them).
- `coal pipeline` no longer writes through a symbolic link in `DIR/work/`. `work/outputs.json.tmp` (the record),
  `state.json.tmp`, `inputs/sources.json.tmp`, the logs (`work/logs/*.log`), the decompressed and `--prefix-ids`
  copies of the inputs, `proteins.sorted.faa.tmp`, `*.kept.faa.tmp`, and the files `coal-cpp translate` and DIAMOND
  write there (`*.faa`, `*.genes.tsv`, the hit tables, `subject.dmnd` / `query.dmnd`) were opened as they were
  found: a link planted with such a name made the next run write into its target, outside `DIR` (the record, the
  logs, the proteins), and the renamed record became that link. They are now created anew like the tables in
  `DIR` (a link or hard link with the name is removed, never followed; `O_CREAT | O_EXCL | O_NOFOLLOW`), a log is
  appended to only when it is no link and has no other name, and `--resume` no longer reuses a stage whose output
  is a link. `coal-cpp translate` and DIAMOND write into a folder of their stage, `work/NAME.staging` (made anew
  for the stage, mode 0700, removed after it and by the next run after a kill), and their files are renamed into
  place (a link put at a name after COAL removed it, while DIAMOND searched, was followed); the
  decompressed copy of a `.gz` input is created with `O_EXCL | O_NOFOLLOW`. A folder COAL writes into that is a
  symbolic link (`DIR/work`, `work/logs`, `work/inputs`, an alignment folder) is refused with exit 1 before
  anything changes, also with `--overwrite` (it was written through: `O_NOFOLLOW` guards only the last part of a
  path); to keep the working files elsewhere, give `-o` a folder there. The programs' command lines in `report.json`
  name those staging folders (the same paths in every run). (Found by the audits of #99 and #102.)

### Added

- **coalr is a release asset**: `.github/workflows/release.yml` has a `coalr` job (`scripts/sync_r_sources.py
  --check`, DESCRIPTION at the version, `R CMD build`, `R CMD INSTALL`, the installed package's version and one
  example run) that adds `coalr_<version>.tar.gz` and the version-free `coalr_src.tar.gz` to the release; both
  are in the release's exact file list, in `release-smoke.yml`'s check (with the tarball's DESCRIPTION and
  contents) and in `tests/test_docs_site.py`. Not for a prerelease: `0.6.1rc1` is no R package version
  (`scripts/release_version.py` `r_version()` is empty for it), so a prerelease has no coalr file.
- The checksums and the exact file list are a job of their own (`files`), which runs on a dry run too
  (dispatched with `publish` off), so a dry run shows the list a release would publish; the `release` job
  publishes exactly the files that job checked.
- `coal_dnds_pairs()` / `coal_dnds_clusters()` (coalr): attribute `"alignment"` with the protein-alignment scoring
  used; docs/r.md, coalr's README and the help pages say that the default scoring is EMBOSS needle's (as in `coal
  pipeline`) and that `preset = "pwalign"` gives orthologr's alignments, with which dS can differ by more than 10%
  for some pairs (118 of 2,958 in tutorial 4). (Final review of 2026-11-16, L3.)
- `nightly.yml` job `msa-tools`: Clustal Omega and MAFFT installed with no fallback, and the cluster-mode MSA tests
  must not skip for want of them (ci.yml's install stays best effort). (Final review of 2026-11-16, L4.)
- **GY is validated on every genetic code** (NCBI 1 to 6, 9 to 16, 21 to 23; the audits of 0.6.0 covered the
  standard code and three pairs of code 2). No change was needed: per code, the stop codons are out of the
  chain's states and of the F3x4 frequencies, synonymous / nonsynonymous follows the code, and columns with a
  stop codon of the code are skipped. `benchmarks/kaks/validate_gy.py codes`: on 40 random pairs per code (680)
  the reported lnL equals the likelihood recomputed with numpy from the model's definition under that code to
  1e-11 relative, with S, dS, dN from the rates, and is a maximum; against codeml with the matching `icode`
  (the eleven tables codeml has), COAL's lnL is never lower on 132 pairs and omega, kappa and t agree to codeml's
  printed digits. Tests: `tests/test_kaks_gy.py::TestGeneticCodes` (20 pairs under ten non-standard codes with
  codeml's stored output, `tests/kaks/make_gy_code_fixtures.py`; the likelihood from the definition under
  codes 2, 4, 11, 12, 14, 16, 21, 22, 23; the stop codons of codes 1, 2 and 4). docs/kaks.md, "Genetic codes".

### Internal

- `coal pipeline`: the record of whose files are in the output folder (`work/outputs.json`: what an earlier run
  wrote, the checks before a run changes anything, the clean-up at its start, the alignment-folder lists around the
  engine) is its own module, `src/coal/outputs_record.py` (class `OutputsRecord`), out of `dnds_pipeline.py`
  (4,447 lines, 86 methods in `_Run`), with unit tests that need no engine (`tests/test_outputs_record.py`: the
  user's `.fasta` kept, a killed run resumed, a foreign `work/` and `--report` refused, tampered lists, 0.6.0
  records, `.FASTA`). Behaviour unchanged but for the fix above: the same files, byte for byte, on `examples/dnds`
  and the cluster fixture. Comments in `dnds_pipeline.py` cite the pull request of a change instead of review
  documents that are not in the repository.
- `coal.outputs_record`: `OutputsRecord.check_before_start()` is the sequence of checks `_Run.start()` makes
  (`--report` path, the run's inputs, earlier outputs and alignment folders, output names, `work/`, `--report FILE`),
  so the unit tests run that same sequence; `open_new` / `open_append` / `remove_file` write without following
  links. The alignment-folder options of `_Run` come from `outputs_record.ALIGNMENT_DIRS`. Tests for four rules no
  test caught when broken on purpose (a record's `read_in_place` voids the report's list, a broken link at a
  recorded name, a `work/logs` link, a named pipe as `work/outputs.json`).

## 0.6.0 — 2026-11-15

**In short.** This release makes COAL more than a PAL2NAL replacement, and changes its licence.

- **Ka/Ks (dN/dS) built in:** `coal dnds` (a tidy table with standard errors and a status per pair) and
  `coal kaks` (KaKs_Calculator's format), with ten methods: YN00 (default), NG86, LWL85, MLWL, LPB93, MLPB, MYN,
  Comeron, Li, and GY, the maximum-likelihood estimate of Goldman & Yang (1994) as PAML codeml computes it for a
  pair (`coal dnds --ml-columns` adds kappa, t and lnL; in coalr `methods = "GY"`, feature `gy94_ml`). Checked
  against KaKs_Calculator 2.0, PAML's yn00 and codeml, and against simulated truth (docs/kaks.md).
- **One step from coding sequences to a dN/dS table:** `coal pipeline pairs` (the orthologs of two species, found
  with DIAMOND) and `coal pipeline clusters` (gene families: DIAMOND DeepClust, then Clustal Omega or MAFFT), with
  a JSON report, `coal pipeline install` for the external tools, and `coal cite` for what to cite
  (docs/pipeline.md). Also `coal translate`, and COAL's own protein aligner `coal-cpp-align`, which the platform
  wheels ship next to `coal-cpp`.
- **R package `coalr` 0.6.0** (packaging/r/coalr, docs/r.md): COAL's C++ core inside R (Ka/Ks, codon alignment,
  the pair and cluster steps of the pipeline; no Python needed), and the whole pipeline through reticulate. The
  orthologr fork's `dNdS(engine = "coal")` uses it.
- **App:** the web GUI and the Windows app have **Ka/Ks** and **dN/dS pipeline** tabs; the MCP server has the tools
  `dnds_pairs` and `pipeline_pairs`.
- **Windows:** `coal pipeline install mafft` (x86-64) installs MAFFT's own Windows package, pinned by its
  SHA-256, so `coal pipeline clusters --aligner mafft` runs on native Windows.
- **Tutorials** on real, small, public genomes, from coding sequences to an interpreted dN/dS table
  (docs/tutorials/); CI reruns four of them.
- **Licence: GPL-3.0-only from this release** (0.5.0 and earlier stay MIT), because the Ka/Ks engine is derived
  from GPL-3 code. Running COAL puts no condition on your data or scripts; see the first entry under "Changed".
- **For users of 0.5.0, the existing commands behave as before** except for what "Changed" lists: the `-o` lock
  (two runs onto one output file are refused), tools in COAL's cache need their checksum record (re-install them
  once: `coal pipeline install <tool> --force`), `pal2nal` without arguments is a usage error, and `coal --help`
  ends with the list of the other commands.
- **What each install gives:** the platform wheels have everything; the source package has the Python engine
  only (codon alignment), and `coal kaks`, `coal dnds` and `coal pipeline` say so there (docs/install.md).

Much of what follows was found by the stress test of 2026-10-23, the full verification of 2026-10-24, two
independent reviews (2026-11-12 and 13), the first CI runs on Windows and macOS, and a review of the project as a
new user would meet it (2026-11-13). Fixes to features that were never in a release are part of that feature's
entry under "Added", not under "Fixed".

### Added

- **Ka/Ks method `GY`: maximum likelihood under the codon model of Goldman & Yang (1994)** (`coal kaks -m GY`,
  `coal dnds -m GY`, `coal pipeline ... -m GY`, coalr `methods = "GY"`; aliases `GY-HKY`, KaKs_Calculator's name
  of its GY rows, and `GY94`). The branch length t, kappa and omega are estimated from the codon pairs with F3x4
  codon frequencies, as PAML's codeml does for a pair (`runmode = -2`, `CodonFreq = 2`), and dN and dS are those
  of Yang & Nielsen (2000), so dN/dS is omega; every genetic code. The simulation study of 2026-11-14 showed YN
  and MYN biased by 10 to 40 % from a dS of about 1 while codeml's likelihood stayed within 8 %: `GY` is that
  estimate, in process. COAL's own code (`cpp/src/kaks_gy.cpp`; only the symmetric eigen solver is shared with
  YN): the rate matrix is decomposed once per (kappa, omega), t is maximised out by Newton steps, and kappa and
  omega by a BFGS iteration with the exact gradient, from fixed starting points (four for a diverged pair, one
  of them on the ridge of saturation), so the result is deterministic, independent of the thread count, of the
  order of the two sequences and of the order of the columns. **Validation** (`benchmarks/kaks/validate_gy.py`,
  docs/kaks.md "GY: maximum likelihood"): on 639 pairs (320 simulated, 160 yeast, 159 *Salmonella* /
  *Shigella*) COAL's log-likelihood is never below codeml's (beyond 1e-4) and is higher for 20 (codeml at
  a bound or at a lower maximum); where both are at the same maximum, dN, dS, omega, kappa and t are codeml's to
  its printed precision. Against simulated truth (6,400 pairs of 300 codons) the median error of dS and of dN/dS
  is within 6 % at every dS from 0.1 to 3, and dS ± 1.96 SE holds the truth for 91 to 97 % of the pairs. **Output:** `coal dnds --ml-columns` (new) adds the columns `kappa`,
  `t` and `lnL` after `note`; `coal kaks` fills ML-Score, AICc and Model (`HKY`) for GY rows, which come after
  MLPB as in KaKs_Calculator. **Standard errors** of dN and dS from the observed information and the delta
  method (codeml's `getSE = 1`; checked to five digits against an independent numpy computation, and within
  about 15 % of codeml's own, whose Hessian is cruder); `NA` for a saturated pair and for the distance that an
  omega at its bound makes 0. **Statuses:** an estimate of dS or dN above 10
  substitutions per site (always the case when t ends at its bound of 50) is `NA` with status `saturated`, and
  the other estimate is kept; a pair whose likelihood is flat in t up to the bound (one codon, every site
  different) is `saturated` without an estimate; identical sequences give dN = dS = 0. **Speed:** 23 to 47 CPU seconds per 1,000 pairs of 300 to 400 codons, 15 to 25 times YN's time and 3 to 4
  times faster than codeml. No other method's
  output changes (byte-identical at `--digits 17`, both output formats, all nine methods, on 4,388 yeast and
  3,685 *E. coli* pairs). `-m GY` used to be refused ("does not implement natively"); MS, MA and the gamma
  variants still are. `coal cite` lists Goldman & Yang (1994) for it. Tests: `tests/test_kaks_gy.py` (20 pairs
  with codeml's stored output, three under the mitochondrial code; the likelihood recomputed from the model's
  definition; the standard errors against an independent numpy computation; edge cases; formats), `tests/test_kaks_properties.py` (GY in the property tests: P1 to P8 and
  P10), coalr `test-gy.R`.
  **Small probabilities** (audit of #92, MEDIUM-1): a pattern probability below 1e-8 (codons two or three
  substitutions apart at a t below about 5e-4: thousands of nearly identical codons with one multiply changed
  codon) is computed by a series of non-negative terms (uniformisation) instead of the eigen sum, whose absolute
  error of 1e-16 was the whole value there: such a pair got status `ok` at the lower bound of t with a dN 100
  times too low (60,000 codons, one triple change: dN 2.0e-7; now 2.4e-5, NG 2.0e-5, and the likelihood equals
  an independent matrix exponential to about 1e-12 relative, where codeml's own is rounding noise).
- **MAFFT on native Windows: `coal pipeline install mafft` installs it** (x86-64; also `--install-tools` and
  `coal.msa.install_mafft()`). It downloads MAFFT's own Windows package, the "all-in-one"
  `mafft-7.526-win64-signed.zip` from mafft.cbrc.jp (MAFFT 7.526, the version COAL pins elsewhere), pinned by size
  (27,476,513 bytes) and SHA-256 (`7c3a508c8512af04007d68019c6fbd4ff00d891244b58f87ec9aa0e53fbae8a7`; MAFFT
  publishes no checksum: computed on 2026-11-15 from two identical downloads), both checked before anything is
  unpacked. Every member name is checked first (one folder `mafft-win/`; no absolute, drive, `..`, backslash or
  link member; the DIAMOND zip's rules for each part of a name; no two names that differ only in case); the tree
  goes into a new folder that replaces the old install, and its manifest records the hash of every file and the
  package (it is unpacked from the same open file whose size and hash were checked). `find_mafft()` then gives
  `<cache>\tools\mafft\mafft-win\mafft.bat`, but a run does not start that batch file: `mafft.bat` expands its
  folder, the current folder, `TMP` and `PATH` again unquoted inside `cmd.exe`, where `&` or `^` in a path would
  run as a command (independent audit of this change). COAL starts the package's `bash.exe` with
  `/usr/bin/mafft` directly, with what the batch file sets (`PATH`, `MAFFT_BINARIES`, a temporary folder of the
  run's own as `TMPDIR` / `TMP` / `TEMP` and as the current folder, the input given by its name there; a run
  folder or install on a network (UNC) path is refused, as MSYS2 misreads such paths), with no standard input,
  and ends its whole
  process tree (`taskkill /T`) on a timeout or interrupt. A `.bat` of your own in `COAL_MAFFT` still runs
  through `cmd.exe`, with a command line COAL builds (`coal.msa.batch_command_line`: the program's full path and
  any argument that is not plain in quotes; `%`, `!`, `"` and a trailing backslash refused), and is refused while
  its folder, the current folder or `TMP` / `TEMP` has one of `& | < > ^ % ! "` in its path. An interrupted
  install's leftovers are removed by the next one, and a failed swap puts the old install back. MAFFT is BSD; the package bundles MSYS2 components and
  MAFFT's page gives GPL for the package; COAL does not redistribute it. CI's Windows C++ job downloads it
  (`tests/test_windows_mafft_real.py`, `COAL_TEST_REAL_MAFFT=1`) with the cache folder, the current folder and
  `TMP` / `TEMP` in paths with spaces, `&`, `^` and `( )`, and runs `coal pipeline clusters --aligner mafft`
  on `examples/dnds`. Clustal Omega on Windows is unchanged (no automatic install; the message's MAFFT option now
  names this install); MAFFT on Windows machines that are not x86-64 gets the message of what works there.
- **Ka/Ks (dN/dS) engine: `coal kaks` and `coal dnds`.** In-process C++ (no external program runs), **derived
  from the source code of KaKs_Calculator 2.0 and of PAML's yn00** (both GPL version 3; attribution in
  `cpp/src/kaks.cpp`, `NOTICE` and docs/kaks.md; this is why COAL is GPL-3.0-only from this release). Methods:
  YN00 (default), NG86, LWL85, MLWL, LPB93, MLPB and MYN, for AXT or aligned FASTA (every pair of an MSA, or
  consecutive pairs), any of COAL's genetic codes, on all threads. `coal kaks -i in.axt -o out.kaks -m YN` takes
  KaKs_Calculator 2.0's `-i -o -m -c` and writes its format (22 columns, NA rules) for these seven methods; it is
  not a drop-in for every KaKs_Calculator call: MS/MA and the gamma variants are refused (GY: its own entry above), and the default
  method is YN where KaKs_Calculator's is MA. `coal dnds` writes a tidy TSV with standard errors and a per-pair
  status (ok, warning for internal stop codons, saturated, no_data, error), so one bad pair never stops a run.
  NG is KaKs_Calculator's NG, whose site counting is not the textbook Nei-Gojobori one (docs/kaks.md).
  Validation on 999 orthologr pairs and 600 simulated pairs (docs/kaks.md): against a KaKs_Calculator 2.0.1 built
  to print 15 digits, Ka and Ks differ by at most 1.4e-13 (YN, MYN) or not at all (NG, LWL family) wherever both
  give a number; against PAML 4.10.10 yn00 built to print 12 digits, run one pair at a time on the columns
  without gaps, ambiguous bases or stop codons and compared with `coal dnds --yn-flavour paml`, dN and dS differ
  by at most 5e-12, which is yn00's printing precision. Every row that differs from the stock KaKs_Calculator's
  text is classified (docs/kaks.md lists the deliberate differences). C++ API `<coal/kaks.hpp>` in `coal::core`;
  the library has no global mutable state and can be called from many threads.

  Settled before the release (states that no release ever had are not listed as fixes):

  - **`--yn-flavour paml` agrees with yn00 at the K80 boundary 2P + Q = 1** (deep test of 2026-11-14, D01).
    When the transition and transversion proportions of the 4-fold (or non-degenerate) sites satisfy 2P + Q = 1
    exactly, yn00 finds K80 undefined while COAL, forming Q in KaKs_Calculator's order, found a residue of 1e-17
    and a kappa of 60 there. The paml flavour now forms Q as yn00's `GetKappa` does. The S. cerevisiae /
    C. glabrata pair of the finding: before S 194.22, N 462.78, dN 0.5024, dS 2.1598 (`ok`); after S 128.0,
    N 529.0, dN 0.4123, as yn00 (dS: yn00 prints 3.6389, COAL `saturated`, as for every pair whose Ks
    correction is undefined). All pairs with `--yn-flavour paml`: 10 of 4,388 Scer / Cgla rows change, 8 of
    4,208 Scer / Klac, 1 of 3,685 E. coli K-12 / Sakai, 0 of 2,862 Salmonella / Shigella; the 10 Scer / Cgla
    pairs disagreed with yn00 in S and N before and agree to yn00's printed digits after (weighting 1 and 0). The
    300-pair samples against yn00 (Scer / Cgla, Scer / Klac, Salmonella / Shigella, K-12 / Sakai, human / mouse,
    A. thaliana / lyrata): every value agrees to yn00's printed digits where both give one (Scer / Cgla before:
    1 pair of 300 differed). The default flavour (KaKs_Calculator's numbers) is unchanged: 0 rows differ on all
    four sets.
  - **Exit status 2 when no pair has an estimate because of the input.** `coal dnds` / `coal kaks` end with 2
    (not 0) when every result is `error` or `no_data`, with one line naming the usual cause: the input looks like
    protein sequences, or the sequences of the pairs differ in length (not aligned), and the way to a codon
    alignment. When only some pairs differ in length the exit status is 0 and a warning says how many. A run with
    any estimate, `saturated` or `not_converged` row exits 0 (docs/kaks.md "Exit status").
  - **`--min-codons N`, `--max-ds X`** (`coal dnds`, `coal codon-dnds`, `coal pipeline`; `min_codons`, `max_ds`
    in coalr): flags for estimates that rest on few codons or a saturated dS. They never drop a row and never
    change a number: the status becomes `warning` and the note says `fewer than N codons compared` / `dS above X`.
    Off by default. docs/kaks.md and docs/pipeline.md have a new section "Filtering before interpretation"
    (short alignments, coverage and identity, dS saturation, very small dS, starting thresholds).
  - **gzip-compressed input** (recognised by its first bytes): the `coal` command decompresses it for
    `coal kaks`, `coal dnds` and `coal translate`; `coal-cpp` itself says "is gzip-compressed: decompress it first
    (gunzip -k FILE)" instead of "control character". Write errors name the file and the cause
    (`write error on 'out.tsv': No space left on device`). The codon alignment (`coal PEP DNA` with either
    engine, `coal-cpp PEP DNA`, the `pal2nal` drop-in) does not decompress: a gzip-compressed protein alignment,
    DNA file or standard input is "`d.fa.gz is gzip-compressed: decompress it first (zcat FILE.gz > FILE, or pipe
    it: zcat FILE.gz | coal ... -): sequence files must be plain text`" (`input_format`; the same bytes in both
    engines) instead of "control character" (protein) or "sequence data before the first '>' header" (DNA).
    bzip2 and xz files are named the same way ("`d.fa.bz2 is bzip2-compressed: decompress it first (bzcat
    FILE.bz2 > FILE, or pipe it: bzcat FILE.bz2 | coal ... -): ...`", "`is xz-compressed`" with `xzcat`), and
    each of the three also after a UTF-8 byte-order mark; bzip2 by `BZh`, the block size and the magic of its first
    block (or of an empty stream), so a text that starts with "BZh" is not taken for one (both engines, the same
    bytes, tests/test_engine_review_fixes.py; follow-up audit 2026-11-15). `coal-cpp kaks` / `dnds` refuse gzip
    only, with their own words, as before.
  - **YN and MYN are faster on close sequences, with the same bits.** The transition-probability matrix is not
    computed for a pair in which no codon differs at two or three positions (nothing reads it then): checked
    byte for byte at 17 digits on the validation sets with GCC and Clang; about half the CPU time on a set of
    closely related pairs, 10 to 20% (YN) on the 999 orthologr pairs (docs/kaks.md "Benchmark").
  - **YN and MYN are faster for every pair, with the same bits** (profiling of 2026-11-15, #85). The 61 × 61
    eigendecomposition behind their pathway weighting (most of their time) runs four independent columns at a
    time, sums only over each codon's neighbours (the terms left out are exactly +0) and reuses one workspace per
    pair; on x86 with GCC or Clang it also has an AVX2 copy (`target("avx2")`) chosen at run time when the CPU has
    AVX2, as the aligner's row loops already were, with the baseline code as the fallback. Both copies perform the
    same operations in the same order (no FMA, no reassociation), so the numbers are the same bits. No build flag
    changed: the wheels, coalr and the Bioconda recipe stay baseline x86-64 / generic aarch64. CPU time: about 10
    to 14% less (YN, MYN) on the AVX2 path on 2,862 *Salmonella* / *Shigella* pairs (`-j 1`, a shared 2-core
    machine), 2 to 6% less on the baseline path; 5 to 9% for a `coal pipeline pairs` run across three data sets;
    less on other data. Every output of the
    pipeline (all nine methods, `--digits 17`, both aligners' alignment files) was byte-identical before and after
    on five data sets. `tests/test_kaks_dispatch.py` checks that the AVX2 copy runs where the CPU has it and that
    both copies give the same bytes (testing options `--kaks-baseline`, `--kaks-isa`, not in the help).
  - **Attribution corrected.** The documentation said the Ka/Ks engine was an independent implementation with no
    copied code. It is derived from KaKs_Calculator 2.0 and PAML's yn00 (GPL v3). The sources, `NOTICE`, the R
    package and docs/kaks.md now say so, with the upstream copyright notices and the list of changes; COAL's licence
    changed accordingly (first entry). `coal kaks` is described as KaKs_Calculator's format for seven of its methods
    with the default YN, not as a drop-in; `-m NG` is documented as KaKs_Calculator's site convention, which is not
    the textbook NG86 one.
  - **Ka/Ks: YN and MYN keep the estimate that has sites.** With no nonsynonymous sites (`GCT` against `GCA`) or no
    synonymous sites (Met/Trp only) both Ka and Ks were NA; the one that has sites is now estimated, with the
    numbers KaKs_Calculator's YN and PAML give (docs/kaks.md, difference 8).
  - **Ka/Ks: MYN gives the defined value where KaKs_Calculator prints `nan`** (a kappa class of weight 0 added as
    0 × NaN; TN93 dividing 0 by 0 when the sites have no T-C or no A-G pair; site base frequencies of a sequence
    without such sites). An identical short pair was reported as "Ka saturated"; it is Ka = 0 now (differences 9
    and 10). Checked on 4,000 random short pairs against the reference and against the reference with those two
    lines corrected.
  - **Ka/Ks: no data race.** Fisher's test took log-factorials from `lgamma()`, which writes the global `signgam`;
    the engine has its own log-gamma now (ThreadSanitizer: one race before, none after; results unchanged).
  - **`coal kaks` / `coal dnds` (and `coal translate`, `coal codon-dnds`) refuse an output file that is the input**
    (exit 2) instead of emptying it, and **an AXT record with one sequence line is that record's error**, no longer
    cut in half and compared with itself.
  - **Ka/Ks: what the comparison with KaKs_Calculator does and does not show** (docs/kaks.md, "Differences",
    rewritten; nothing changes on real genes: the output for the 999 orthologr and 600 simulated pairs is
    byte-identical). The statements "no value that both tools print differs by more than 1e-9" and "the corrected
    reference gives COAL's numbers" were made from 4,000 random pairs in which no codon differs at more than one
    position, and are now scoped to them. On 3,000 short pairs with codons drawn independently YN and MYN print
    different numbers in 33 fields (23 since the next entry but one: where the iteration does not converge COAL
    now prints none), each now explained and tested (differences 11 to 13): KaKs_Calculator drops
    the differences of codons that differ at several positions once a NaN is in its iteration and prints a Ka
    that is too small (COAL is right: it gives yn00's unweighted value, or the corrected reference's); an
    iteration that does not converge is not an estimate in either program; a saturated kind of site changes the
    path of the iteration within its tolerance. KaKs_Calculator's MYN also depends on the order of the pairs in a
    file (49 of 3,000 rows), so the comparison now runs it one pair per process.
    **Not said before: the MYN fix above ("gives the defined value where KaKs_Calculator prints `nan`") also
    changed MYN numbers that COAL used to print** on such short pairs (13 Ka and 3 Ks among the 3,000; for
    example Ka 0.0781 became 0.1673, where YN gives 0.1672): COAL had followed the reference into the lost codons.
  - **Ka/Ks: a YN or MYN iteration that does not converge gives no estimate: status `not_converged`, NA.** It was
    `ok` with a note, then `warning`, with the values of the last round printed, as KaKs_Calculator and PAML do.
    Those values are not reproducible: for the one-codon pair `AAA` / `ATG` the same source printed Ka 3.1688 on
    Linux, 3.1736 on Windows and macOS, and "saturated" with another compiler and C library (CI of 2026-11-13).
    The wandering iteration amplifies a last-bit difference of `exp` / `log`, or a fused multiply-add, into the
    third digit. Now `coal dnds`, `coal kaks`, the pipeline and coalr give `NA` for Ka, Ks, Ka/Ks, the standard
    errors and the weighted difference counts (in `coal kaks` also Fisher's p and Divergence-Time), the status is
    `not_converged` (new; between `warning` and `saturated`, and also where the last round was saturated, since
    that differs by platform too), and the note keeps the last round for inspection. The site counts are kept.
    `cluster_summary.tsv` has a new column `n_not_converged` (after `n_ok`, which does not count these pairs), the
    pipeline report has `counts.pairs.not_converged` (docs/pipeline.schema.json), the summaries on stderr name the
    count, and the C++ `Result` has `last_ka` / `last_ks`. Among the real-gene validation sets this changes one
    row (a simulated pair, MYN: Ka 0.58117 and Ks 1.43957 became NA); the 999 orthologr pairs' output is
    byte-identical. Among 3,000 random short pairs: 123 YN and 87 MYN rows. docs/kaks.md, difference 12.
  - **Ka/Ks: the same digits on every compiler.** The engine is compiled without floating-point contraction
    (a pragma in `kaks.cpp`, `-ffp-contract=off` / `/fp:precise` in CMake), so a build for a CPU with FMA (every
    arm64 Mac) no longer differs from an x86-64 one in the last digits. Eight builds (GCC 13 and 16, Clang 18,
    `-O2` / `-O3`, with and without `-march=x86-64-v3`) are byte-identical on about 67,000 rows. What remains is
    the C library's `exp` and `log`: with another one, no row of the 999 + 600 real-gene pairs changes at the
    default precision, and the last one to three of 15 digits change in 8 to 12% of them. No measurable cost in
    speed. docs/kaks.md, "Platform reproducibility", says what is guaranteed and what is not.
  - **Ka/Ks: Fisher's test with a cell that is zero up to rounding.** When the weighted differences add up to the
    number of sites (Sd = S = 2, or 2.0000000000000004), the test was made or left out by the last bit (p = 0 or
    NA, by platform). The cell now counts as zero (p = 0, as KaKs_Calculator prints). No row of the real-gene sets
    changes.
  - **Ka/Ks, MYN: an infinite kappa is `saturated`, not `no_data`.** A pair whose kappa_TC or kappa_AG is
    estimated as infinite (TN93's transition term of one class of sites exactly at its limit) had status
    `no_data` ("no synonymous sites; no nonsynonymous sites"), and `coal kaks` on such pairs alone exited 2,
    suggesting that the input was not a codon alignment. Now `saturated`, with a note naming the kappa;
    Substitutions keeps its count (docs/kaks.md, difference 15). 79 of 240,000 simulated pairs, 3 real pairs among
    the deep test's genomes. In `coal kaks`'s format the Substitution-Rate-Ratio of such a row is `NA`
    (KaKs_Calculator prints `inf` for the infinite kappa, e.g. `2.30174:inf:1:1:1:1`; COAL's output has no nan /
    inf: property fuzz P3, seed 8585); the dnds table, the C++ API (`Result::rates` keeps the inf) and every other
    field are unchanged. One *Salmonella* / *Shigella* MYN row of the deep test changes (p2067).
  - **Ka/Ks, MLWL: the site formula at kappa exactly 2 or 0.5 no longer depends on the order of the codon
    columns.** kappa = 2 ts / tv is often exactly at MLWL's thresholds, and the floating-point sums put it on
    either side according to the column order, a site or more of N and S apart (property fuzz P6, seed 4242;
    KaKs_Calculator does the same). The comparison is now made on exact sums (exact kappa 2 is "not above 2"); the
    printed values are the same doubles. One MLWL row of the deep test's real genes changes (*S. cerevisiae* /
    *C. glabrata* p2754: Ka 0.249535 -> 0.256528, Ks 1.81172 -> 1.65011, S 64.667 -> 71, the values KaKs_Calculator
    prints for 38 of 40 random orders of its columns); 10 of the 3,000 random short pairs with independent codons
    (docs/kaks.md, difference 18).
  - **Ka/Ks, LWL family: no substitution counts from a saturated correction** (review of 2026-11-15, M1). In a
    saturated class A, B and K are the 0 that stands for the undefined correction, and Sd / Nd were summed with
    those zeros, LPB / MLPB's S = Sd / Ks and N = Nd / Ka divided them (a code-11 pair of 11 codons: LPB S + N =
    42.5 sites). Now Sd (`S-Substitutions`) is `NA` when Ks is saturated and Nd when Ka is, and so are LPB /
    MLPB's S and N on that side and Fisher's p; LWL / MLWL's fold-class site counts, the Fold-S/N-Substitutions
    and the other side's count are kept (docs/kaks.md, difference 19; difference 5 no longer says that
    KaKs_Calculator prints `-nan` for a saturated Ks: it prints the quotient). All 9 methods at 17 digits, both
    formats: *Salmonella* / *Shigella* 530 of 25,758 rows change (LWL 129, LPB 129, MLWL 136, MLPB 136: every
    saturated row of the four methods, nothing else), *S. cerevisiae* / *C. glabrata* 8,694 of 39,492. coalr:
    re-synced; `coal_info()$features` gains `"lwl_saturated_counts_na"` and, for two earlier vendored changes that
    had no name, `"myn_kappa_inf_saturated"` and `"mlwl_exact_kappa"`.
  - **Ka/Ks: a multiple-hit correction whose argument is 0 in exact arithmetic is `saturated`** whatever the
    rounding (an argument at or below 1e-10 counts as 0): a rounding residue printed dN or dS of 19 to 28 as an
    estimate (fuzz finding K03; NG, the LWL family, the YN/MYN distances, default Comeron and Li; the kappa
    estimates and the gestimator / seqinr flavours keep their references' tests). Comeron's "fewer than
    `--comeron-min-sites` sites" no longer depends on the order of the codon columns (difference 14). Standard
    code: 22 rows of 4,388 S. cerevisiae / C. glabrata pairs and 2 Comeron rows of the Salmonella / Shigella pairs
    change; nothing else.
  - **Ka/Ks, genetic codes other than the standard one:** LWL/LPB's Arg special case and MLWL/MLPB's Ile/Arg rules
    follow the selected code (K01, K02; with code 2, ATA to ATG is synonymous). Standard code unchanged.
  - **Ka/Ks, YN/MYN: a row without any estimate no longer prints `Sd` / `Nd` that depend on the order of the
    sequences** (now `NA` where the pathways were weighted; K06).
  - **Accuracy against simulated truth, and filtering advice from it** (docs/kaks.md "Accuracy against simulated
    truth"; `benchmarks/kaks/sim_truth/` reproduces it). YN/MYN are unbiased up to dS about 0.5 and 10 to 40% off
    beyond; the counting methods are biased in both directions by codon usage; reliability limits by gene length,
    and how well the standard errors are calibrated. The filtering advice is now `--max-ds 1` (0.5 with strong
    codon bias; 2 only as an outer bound), at least 100 codons for genes read one by one, and dS of at least about
    0.05 for ratios; docs/pipeline.md and the tutorials follow it (tutorial 2 shows both the outer bounds and the
    stricter values, and its `--min-codons` / `--max-ds` example uses 100 and 1). `--yn-weighting 0` is less
    biased at omega below 1 and dS 1 to 1.5 (the default is unchanged).
  - **`coal kaks -j1`, `coal dnds -j1`: attached values of one-letter options are read**, as the main command
    does (`-j4`, `-j=4`, `-mYN`, `-ofile`, chained flags `-qj4`), in `coal kaks`, `coal dnds`, `coal translate`,
    `coal codon-dnds` and `coal-cpp-align` (which also takes `--jobs=4` now). They were "unrecognized argument"
    while `-j 1` worked. The value of an option is never taken apart (`-o -q1` still writes a file called `-q1`).
- **Ka/Ks: the methods `Comeron` and `Li`**, the two dN/dS methods of the R package orthologr that do not need
  KaKs_Calculator (`Comeron`, alias `C95`: Comeron 1995 as libsequence's gestimator computes it, orthologr's
  default; `Li`: Li 1993 as `seqinr::kaks()` computes it). In `coal kaks` (22 columns, what the method defines),
  `coal dnds` (Li with standard errors), `coal pipeline` (`-m Comeron`), `<coal/kaks.hpp>` (`Method::Comeron`,
  `Method::Li`) and coalr (`coal_kaks()`, `coal_dnds_pairs()`, `coal_dnds_clusters()`); standard genetic code
  only. Both are ports: `cpp/src/kaks_comeron.cpp` is derived from libsequence (Copyright (C) 2003-2009 Kevin
  Thornton, GPL version 3 or later; orthologr's copy, modified by Sarah Scharfenberg and Hajk-Georg Drost under
  GPL version 2 or later) and `cpp/src/kaks_li.cpp` from `kaks.c` of seqinr 4.2-44 (`GPL (>= 2)`); the file
  headers, `NOTICE`, coalr's `inst/COPYRIGHTS` and `Authors@R` say so. Two flavours each: COAL's policy by
  default (terminal stop dropped, internal stops skipped, `NA` with status `saturated` where the correction is
  undefined), and the reference's own reading of the pair (`--comeron-flavour gestimator`, `--li-flavour
  seqinr`; in the pipeline and coalr with `--terminal-stops compare` for orthologr 0.4.3, which compares the two
  stop codons). Validation (docs/kaks.md, "Comeron and Li"): on 14,153 pairs (999 + 984 real, 600 simulated,
  570 fixtures, 11,000 random short ones) every number that gestimator's code or `seqinr::kaks` returns is the
  same double in the reference flavours, and every 999 / 9.999999 is `NA`; with orthologr 0.4.3's own
  alignments of its 1,000-gene example coalr gives orthologr's table for all 984 pairs (at the 4 digits
  gestimator prints). From the coding sequences alone it does not: orthologr 0.4.3 aligns proteins without a
  substitution matrix, COAL with BLOSUM62, and 540 of the 984 protein alignments differ. `Li` is not `LPB`
  (other site classes and pathway weights; the numbers differ). libsequence's count for codons that differ at
  two positions (it passes its arguments in another order than its function takes them) is kept in both
  flavours and documented.
  With the merge of the review fixes of 2026-11-13: `kaks_comeron.cpp` and `kaks_li.cpp` are compiled without
  floating-point contraction like `kaks.cpp` (pragma and build flags; the vendored copies in coalr too), so the
  two methods print the same bytes whether or not the target has FMA: GCC 13 and Clang 18, default and
  `-march=x86-64-v3`, with and without the build-system flag, `-O0`: 0 of 169,828 rows differ on the 14,153
  validation pairs, where 3,045 of 28,306 reference-flavour rows differed at 17 digits with contraction
  (docs/kaks.md, "Platform reproducibility"). The status `not_converged` (YN and MYN only) is in every status
  list, count and summary that Comeron and Li rows go through; they never have it.
- **One step from coding sequences to dN/dS: `coal pipeline` (docs/pipeline.md).**
  `coal pipeline pairs --query q_cds.fa --subject s_cds.fa -o out` translates both files, finds orthologs with
  DIAMOND (reciprocal best hits, or best hits; or takes a pairs file), aligns each pair's proteins (COAL's global
  aligner by default; Clustal Omega, MAFFT or DIAMOND's own alignment as options), builds the codon alignment and
  computes Ka/Ks (YN00 by default, several methods allowed). `coal pipeline clusters --cds all_cds.fa -o out`
  clusters with DIAMOND DeepClust (or takes its two-column clusters file), aligns each cluster with Clustal Omega
  (two members: the global aligner), and computes Ka/Ks for every pair of a cluster, with a per-cluster summary
  (counts by status; median, mean, quartiles, minimum and maximum of dN, dS and dN/dS). Clusters above
  `--max-cluster-size` (200) are computed against the representative, sampled reproducibly (`--pair-sample`,
  `--seed`) or skipped. Output: `dnds.tsv`, `cluster_summary.tsv`, `genes.tsv`, `failed.tsv` and a JSON report
  (`docs/pipeline.schema.json`); optionally KaKs_Calculator's format and the codon alignments. One bad gene, pair or
  cluster never stops a run (exit 2 when something failed). No process or temporary file per pair: the new native
  steps `coal-cpp translate` and `coal-cpp codon-dnds` (C++ API `<coal/dnds.hpp>`) do translation, pair alignment,
  codon threading and Ka/Ks in one threaded process fed through a pipe; output order does not depend on `-j`.
  `coal pipeline install diamond` installs DIAMOND into COAL's cache, like Clustal Omega (pinned by SHA-256). On orthologr's
  1,000-gene example (999 pairs): 2.5 s on 2 cores, every value equal to the steps run one by one and to
  KaKs_Calculator 2.0 on the same codon alignments; 19,980 pairs in 20.5 s with 41 MB.

  Settled before the release (states that no release ever had are not listed as fixes):

  - **Every FASTA record is a gene, and the run says so when that looks wrong.** Input with several transcripts
    per gene (found from `gene:`, `[gene=`, `locus_tag=`, `gene_symbol:` header fields, or from shared ID stems
    such as `AT1G01010.1` / `.2`) gets one warning with counts. `--longest-isoform` keeps the longest CDS per gene
    before the search (`--gene-map TSV` when the headers do not name the gene); dropped records are in `genes.tsv`
    with status `dropped`, and `report.json` records the choice (docs/pipeline.md "One record is one gene").
    Second review (2026-11-13): stable gene identifiers are read before gene symbols (`locus_tag`, `GeneID:`,
    Ensembl `gene:`, FlyBase `parent=FBgn`, Phytozome `locus=`, then `[gene=` / `gene_symbol:`), and records that
    share a symbol but are clearly different loci (other contig, locus tag or GeneID: several strains in one
    file, two copies of `tuf`) are not merged, with a warning. More ID styles: TransDecoder `_i1.p1`, MAKER
    `-RA` / `-mRNA-1`, maize `_T001`, RAP-DB `-01`, Phytozome proteins `.1.p`; RefSeq `NM_` / `XM_` accessions get
    a warning that the ID holds no gene. `--gene-map` wins for the records it names and the others fall back to
    the headers and IDs (it used to switch detection off for every file); a map that names no record of the
    input, or gives a record two genes, is an error (exit 1). **`--longest-isoform` never drops more than half
    of a file's records on the evidence of ID suffixes alone**: it stops with exit 1 and asks for `--gene-map`
    or the new `--isoform-rule suffix` (`map`: use the map only). docs/pipeline.md has the table of header
    styles (one test file each) and how the heuristics fail (`HSP70.1 / HSP70.2`, `scaffold7.1`). The progress
    and summary lines count "records" until `--longest-isoform` has made genes of them, and a repeated ID is
    counted as not usable in both.
  - **A file is complete or absent** (second review, 2026-11-13). Every table and alignment file of the pipeline
    is written as `NAME.tmp`, created anew, and renamed when it is complete; before, only `dnds.tsv` was, and a
    stage that failed on a full disk left a cut-off `codon_alignments.fasta`, `protein_alignments.fasta` and
    `dnds.kaks` under their final names. After a failed or interrupted stage none of them is there, the report
    lists only the files that exist, and the `.tmp` files are removed (the next run removes what a killed run
    left). A symbolic link in the output folder under an output's name is replaced, never written through
    (`dnds.kaks -> another file` used to overwrite that file). `coal-cpp codon-dnds` has `--new-files` for this
    and writes each per-cluster file through `NAME.fasta.tmp`. docs/pipeline.md "Output".
  - **A hint for genes of another genetic code** (second review, 2026-11-13). Mitochondrial genes in a nuclear
    CDS set got `internal_stop` warnings and estimates from mistranslated codons, and nothing said why. When at
    least two genes of a set have two or more internal stop codons under `-c` and none under another table COAL
    supports, the run prints one warning per set (how many genes, which table, three examples, the two-run
    recipe); `report.json` has it as `genetic_code_hints`. docs/pipeline.md "Genes of another genetic code":
    mitochondrial and plastid genes, the recipe, and that `Comeron` / `Li` are standard-code only. A per-gene
    code option (`--code-map`) was evaluated and not built: the engine has one code per run; the design is in
    docs/development.md.
  - **Smaller things from the second review** (2026-11-13). `--diamond-args`, `--deepclust-args` and
    `--realign-args` refuse the DIAMOND options COAL sets itself (`--out`, `--db`, `--query`, `--outfmt`,
    `--threads`, ...; exit 2 before anything runs: `--diamond-args '--out x.tsv'` used to send the hits elsewhere
    and end with an unrelated error), all three are documented (`--realign-args` was not), and `realign.tsv`
    always has its header line. The hint about IDs of `--pairs` / `--clusters` that match no record comes when
    half of one side is unknown or no pair can be made (it needed 90% of one side). Messages name the file the
    user gave, not its decompressed copy (`coal pipeline`, `coal kaks`, `coal dnds`). The default `-j` is one
    thread per CPU the process may use (its affinity: `taskset`, SLURM, cgroups), and `coal-cpp` is never given
    more threads than that either (`report.json` `threads`: `requested`, `usable_cpus`, `engine`, `tools`).
    DIAMOND's version is read as the first dotted number of `diamond version` (`2.1` is 2.1.0, not "too old"),
    and a version that cannot be read gives one warning. A gzip input has `sha256_uncompressed` next to `sha256`
    in the report.
  - **`--min-identity P`**, the third flag next to `--min-codons` and `--max-ds` (`coal pipeline`,
    `coal codon-dnds`; `min_identity` in coalr): a pair whose protein identity (the column `prot_identity`) is
    below P percent gets status `warning` and the note `protein identity below P%`; no row is dropped, no number
    changes. docs/pipeline.md says why: a cluster can hold pairs that are not related at all, and they get
    status `ok`.
  - **Documentation from the second review:** docs/pipeline.md "How much the aligner matters" (150 *Arabidopsis*
    pairs: 120, 90 and 75 of them have the same values as the default aligner with MAFFT, DIAMOND's alignment and
    Clustal Omega; up to 0.43 in dN/dS for single pairs); docs/kaks.md says what `P-Value(Fisher)` tests and
    that `SE_dN` / `SE_dS` are standard errors of dN and of dS (NG, YN, MYN, Li), with none for the ratio; README
    example 5 (the orthologr fork) passes `sensitivity_mode = "default"`, so that it finds the 12 pairs of
    example 3 (orthologr's default, DIAMOND `--fast`, finds 8 on the example files).
  - **Cluster mode does not depend on the order of the input.** DIAMOND DeepClust is order-dependent, so it is
    given the proteins sorted by ID: the same sequences give the same `clusters.tsv` whatever the order of the
    `--cds` files and records. (Development builds before 2026-11-13 clustered in input order; their clusters can
    differ from this release's.) docs/pipeline.md "What the search and clustering defaults do" documents the
    best-hit rule, DIAMOND's sensitivity, a worked comparison with orthologr's fork, and why default DeepClust
    clusters hold many saturated pairs.
  - **Provenance in `report.json`:** start and end time (UTC), platform, Python and COAL versions, the thread
    counts as resolved, every input file with size and SHA-256 (`--hash-inputs auto|always|never`; above 2 GiB the
    hash is skipped with the reason), the exact DIAMOND and coal-cpp command lines, the aligner's command
    template, the citation keys (`cite`), the run's warnings. The tables are unchanged by it.
  - **Robustness:** gzip-compressed inputs are decompressed into `work/inputs/`; a full disk (any operating-system
    error) ends the run with one line, exit 1 and a report that is complete or absent; a `--pairs` / `--clusters`
    file whose IDs match (almost) nothing gets one line with example IDs of both sides; DIAMOND older than 2.1.0
    stops cluster mode with the reason (pair mode warns; tested with 2.2.8); `-j` above the CPU count is clamped
    for DIAMOND and the aligners; a tool ended by a signal is reported as "killed by signal 9 (SIGKILL; out of
    memory?)"; the same file as `--query` and `--subject` is a warning.
  - **`coal pipeline ... --keep-protein-alignments`** writes the protein alignments the codon alignments are made
    from: `protein_alignments.fasta` in pair mode (records 1+2, 3+4, ...), one file per cluster in
    `protein_alignments/` in cluster mode, with the same records, order and file names as `--keep-alignments`
    writes, so a codon row is three times its protein row, column for column. Until now the protein alignments
    could not be obtained from the pipeline (review of 2026-11-12). The folder is recorded in `work/outputs.json`
    (`protein_alignments_dir`), so the next run removes that run's files and nothing else; a `protein_alignments/`
    that no COAL run wrote is never written into. The engine step has `coal codon-dnds --protein-out FILE` and
    `--protein-dir DIR`; the report's `outputs` has `protein_alignments` (docs/pipeline.md, the schema).
  - **`coal pipeline`: one run per output folder** (`<out_dir>.coal-lock`; a second run is refused), **an earlier
    run's files are removed when a run starts** (no `dnds.kaks` or `codon_alignments.fasta` of another run next to
    a new `dnds.tsv`), **an interrupt ends the engine and is reported** (SIGINT/SIGTERM to the Python process kill
    `coal-cpp codon-dnds` and the aligners; `report.json` gets `"status": "interrupted"`, a new report field, and
    the exit code is 130 or 128 + signal), **a protein file is explained before DIAMOND runs**, and **a non-ASCII
    sequence ID is that gene's failure** (`id_not_ascii`), not the end of the run (docs/pipeline.md).
  - **`coal pipeline` no longer deletes its inputs or files it did not write.** The clean-up at the start of a run
    (previous entry) removed every file with an output's name, also when it was this run's input or the user's own:
    `--pairs out/orthologs.tsv -o out` failed with "cannot open" after deleting the pairs and the earlier
    `dnds.tsv`, and `clusters --clusters proj/clusters.tsv -o proj` deleted the clusters file and every
    `proj/codon_alignments/*.fasta`. Now a run removes only what an earlier run is known to have written (the
    list each run keeps in `work/outputs.json`, or the outputs named in an older `report.json`), never one of its
    inputs, and it never writes onto an input: `--pairs out/orthologs.tsv` is read and left as it is, any other
    input under a name the run writes is refused before anything changes, and so is `--keep-alignments` into a
    `codon_alignments/` folder that holds files COAL did not write (docs/pipeline.md).
  - **`coal pipeline`: one rule for text, on every system** (docs/pipeline.md, "Text encoding"). Every table is
    read and written as UTF-8 with `\n` line ends, never with the system's default encoding or line ends; sequence
    IDs are bytes that are not re-encoded. Bytes that are not valid UTF-8 are shown as `\xhh` (an ID in Latin-1:
    `g\xe8ne`) by `coal-cpp translate` and by the Python layer alike: before, `coal-cpp` wrote `?` for each such
    byte and Python U+FFFD, so two different IDs could look the same. `coal-cpp-align` writes its alignment to
    standard output in binary mode (`\n` on Windows, as its `-o` file). The Windows CI failure itself was the
    test reading `failed.tsv` with the default encoding (cp1252); the tests now run the pipeline under an ASCII
    locale, with a cp1252 console and with Python's `-X warn_default_encoding` as an error, and compare every
    output byte with the UTF-8 run.
  - **Deep test of 2026-11-14 (D01-D08).** **The same IDs in two `--cds` files no longer drop a species
    without a word** (D05). `clusters --cds a.fa b.fa` with IDs repeated across the files (BRAKER / AUGUSTUS
    `g1.t1` in every species) kept the first file's record of each ID and left out the other's: 300 + 300
    Salmonella / Shigella CDS named `gene0`..`gene299` in both gave 289 clusters and 12 pairs instead of about
    300 cross-species pairs, exit 2, and 300 failures saying "the ID occurs more than once in the file" (each
    occurred once per file). Now such a run stops before anything runs (exit 1) and names both files, the count
    and examples; the new `--prefix-ids` puts each file's name and `|` before its IDs (`q_sameid|gene0`; the
    report's `inputs` give each file's `id_prefix`), which gives 600 usable records, 291 clusters and 352 pairs
    on the same files. `coal-cpp translate` with several files says "the ID occurs in two input files, A and B"
    for such a record (within one file the wording is as before).
    **The genetic-code hint no longer names a wrong code** (D02). It named the code that made the most flagged
    genes stop-free, and codes that also read TAA / TAG "clean" frameshifted pseudogenes: it fired in 6 of the 9
    real pair datasets of the deep test and was wrong in each (S. cerevisiae mitochondrial genes: 14, should be
    3; C. glabrata / K. lactis mitochondria: 4, should be 3; human and mouse Ensembl: 4, should be 2; E. coli K-12,
    A. lyrata NCBI and Ensembl mouse all-transcript files: 14 for pseudogenes and NMD transcripts). Now
    pseudogenes (`[pseudo=true]`, `*pseudogene` biotype) and frameshift-length CDS are not evidence, genes
    scattered over nuclear chromosomes are not either when the headers name the molecule, the code whose extra
    stop or unused codons the genes avoid wins (AGA / AGG for 2, CGA / CGC for 3) before the one with the fewest
    unneeded stop reassignments, and when the sequences cannot tell the codes apart all candidates are named and
    none is recommended (report: `clean_under` null, new `candidates` and `evidence`). After, on the same 9
    datasets: Scer / Cgla and Scer / Klac: code 3 for both sets (17, 8 and 7 genes); human / mouse whole Ensembl
    sets: code 2 (12 genes each); no hint for Salmonella / Shigella, K-12 / Sakai, the A. thaliana / lyrata and
    human / mouse subsets, A. thaliana / lyrata whole genomes, the human / mouse all-transcript files. Also
    A. gossypii: 3 (was 4); Y. pestis (78 genes, almost all `[pseudo=true]`): none (was 14); D. melanogaster
    Ensembl: "4 or 5" (was 14; 5 is right). A synthetic set of four nuclear genes given two internal TGA each
    now gets the candidates 3, 4, 5, 9, 10, 13 and 21 instead of "code 4".
    **A byte order mark or classic Mac line ends no longer lose genes** (D03). `coal-cpp translate` skipped a
    UTF-8 byte order mark and split lone `\r` line ends, the Python-side readers did not: with a BOM the first
    query gene was "not in the query file" (300 Salmonella CDS: 1 pair error, exit 2), with `\r` line ends every
    pair was "the CDS is empty" (300 errors), the first line of a `--gene-map` was silently ignored (4 records
    dropped instead of 5, exit 0), and the first member of a `--clusters` file was dropped with an invisible U+FEFF
    in its cluster ID (exit 2). Now the FASTA index, the header reader of `--longest-isoform`, the `--pairs`,
    `--clusters` and `--gene-map` readers and `coal --batch` job files read such files as coal-cpp does; on the
    deep test's files every output is byte-identical to the run on the plain file (exit 0).
    **Messages** (D04, D06, D07, D08; no number changes). A `coal --batch` job line that starts with `coal` (as
    the README's wording suggested) is refused before any job runs, naming the fix, instead of 400 errors "cannot
    open 'coal'"; the README says "without the program name". A pipeline run with no pair says so ("no pairs:
    DIAMOND found no reciprocal best hit ..." / "none of the N clusters has two usable members"; exit status
    unchanged), a `--pairs` file whose columns look swapped and one with repeated pairs get one warning each (the
    pairs are still computed as given), and `$COAL_CLUSTALO` / `$COAL_MAFFT` / `$COAL_DIAMOND` set to something
    that is not a file is named as the cause instead of "not found (PATH, $COAL_CLUSTALO, COAL's cache)".
    `coal dnd x.fa` (and `coal-cpp dnd x.fa`) says "did you mean `coal dnds`?". `coal pipeline install mafft`
    installs MAFFT 7.526 (pinned, as Clustal Omega 1.2.4 is; it took whatever Bioconda had) and prints "MAFFT
    7.526 installed". When more than 30% of the pairs are `saturated` or `not_converged` the run says that the
    sequences are too divergent for dS and points to "Filtering before interpretation" (S. cerevisiae /
    C. glabrata: 3,775 of 4,388, 86%).
  - **Retest of 2026-11-14 (R01-R04).** **The genetic-code hint names the code of a whole set read with the
    wrong `-c`** (R01). *Mycoplasma genitalium* (524 CDS, code 4) with `-c 1` or `-c 11` got no hint at all
    although 376 genes have internal TGA: the rule that keeps only molecules with at least half of their records
    as evidence dropped its one chromosome (253 of 524). *Tetrahymena* with its records shuffled got "81 of the
    26,927 genes ... mitochondrial or plastid genes in a nuclear set? ... the rest with -c 1 and these with -c 6"
    (the first 5,000 flagged genes looked at, spread over scaffolds), i.e. advice to keep 26,846 genes under the
    wrong code. Now, when at least half of a set's records have internal stops, the hint is about the whole set
    ("most genes (376 of the 524) have internal stop codons under genetic code 11, and 253 of the 264 with two
    or more have none under NCBI genetic code 4 ... The whole set looks like code 4 ... Run again with -c 4"),
    the molecule rule does not apply, a molecule with 150 records or more is a genome (only the codes of whole
    genomes are candidates: 1, 4, 6, 10, 11, 12, 15), and of codes that tie those whose start codons do not fit
    a tenth of the genes drop out (Mycoplasma starts 34 of 264 genes with GTG / TTG / TTA, which code 10 does
    not allow). A set with more than 5,000 flagged genes is judged on 5,000 spread evenly over the file, not the
    first ones. `report.json`: new `records` and `whole_set` in `genetic_code_hints`. The organelle wording and
    the two-run recipe are kept for sets in which fewer than half of the records have internal stops. Before
    (0153e3b) and after, the warning of `coal pipeline pairs` for the query set, on every real CDS set of the
    deep test, the retest and the D02 evaluation (43 runs); 33 are unchanged (listed below the table):

    | Set | `-c` | Before | After |
    |---|---|---|---|
    | *M. genitalium* RefSeq (524) | 1 | none | whole set: code 4 (376 of 524 with internal stops) |
    | *M. genitalium* RefSeq | 11 | none | whole set: code 4 |
    | *T. thermophila* RefSeq (26,996) | 1 | "4,896 of the 26,927 genes": code 6, organelle wording | whole set: code 6 (4,966 of 5,000 looked at) |
    | *T. thermophila*, records shuffled | 1 | "81 of the 26,927 genes": code 6, organelle wording | whole set: code 6 (4,971 of 5,000) |
    | *P. tetraurelia* RefSeq (39,642) | 1 | "4,850 of the 39,534 genes": code 6, organelle wording | whole set: code 6 (4,995 of 5,000) |
    | human Ensembl (19,746) | 2 | "4,988 of the 19,311 genes": codes 1, 11, 12, "mitochondrial genes in a nuclear set?" | whole set: codes 1, 11 or 12 |
    | *S. cerevisiae* mitochondrial genes alone (19) | 1 | code 3, organelle wording | whole set: code 3 |
    | *C. glabrata* mitochondrial genes alone (11) | 1 | code 3, organelle wording | whole set: code 3 |
    | *K. lactis* mitochondrial genes alone (9) | 1 | codes 3, 4, 5, 13, organelle wording | whole set: codes 3, 4, 5 or 13 |

    Unchanged: *S. cerevisiae*, *C. glabrata*, *K. lactis*, *A. gossypii* whole genomes (`-c 1`: code 3 for
    their mitochondrial genes); human and mouse Ensembl (`-c 1`: code 2, 12 genes each); *D. melanogaster*
    Ensembl (twice) and RefSeq (`-c 1`: codes 4 or 5, 11 genes); no hint for *Salmonella* (twice), *Shigella*
    (twice), *Y. pestis* (twice), *E. coli* K-12 and Sakai (`-c 11`), *A. thaliana* (whole, 1,000, TAIR10
    RefSeq), *A. lyrata* (whole, 1,000), the human giant-gene and exon sets, the human and mouse subsets and
    all-isoform files; none with the right code (*M. genitalium* `-c 4`, *T. thermophila* and *P. tetraurelia*
    `-c 6`, *S. cerevisiae* mitochondrial genes `-c 3`), and none where the code reads fewer codons as stops than
    the genes' own, so that no gene has an internal stop to explain (*S. cerevisiae* genome `-c 3`,
    *D. melanogaster* `-c 5`). No run names a wrong single code.
    **A `--gene-map` in prefixed form no longer merges the genes of two species** (R02). With `--prefix-ids`, a
    map naming the prefixed records (`q_sameid|gene0<TAB>G0`, `s_sameid|gene0<TAB>G0`) kept the gene ID as it
    was, so `--longest-isoform` took the two species' records for transcripts of one gene `G0` and dropped one
    of each (exit 0), although docs/pipeline.md says `S_cer|GENE1` and `C_gla|GENE1` are two genes. The gene
    ID of a prefixed record now gets the record's prefix too (unless it already has it), as for a map in
    unprefixed form. The retest's files (300 + 300 CDS, IDs `gene0`..`gene299` in both, a map of 4 lines):
    before 598 genes, 2 records dropped, 350 pairs; after 600 genes, 0 dropped, 352 pairs, `dnds.tsv`
    byte-identical to the run with the map in unprefixed form.
    **An unpinned MAFFT in COAL's cache is named** (R03). A MAFFT that an earlier COAL installed unpinned
    (7.525) was kept by `coal pipeline install mafft`, which printed only its path (exit 0), and runs used it
    without a word (the report's `tools` had the version). It is still kept (replacing it is the user's
    decision), but `install` and every run that uses it print one note: "MAFFT 7.525 in COAL's cache (...) is not
    the pinned 7.526 (an earlier COAL installed it unpinned); it is used as it is. To install 7.526: coal
    pipeline install mafft --force"; in a run it is also in the report's `warnings`. The same for a cached
    Clustal Omega that is not 1.2.4. A MAFFT on PATH or in `COAL_MAFFT` is not checked against the pin.
    **The self-pair warning only where it is true** (R04). "--query and --subject are the same file ...: every
    gene's best hit is itself, so the pairs are self-pairs with dN = dS = 0" also came with `--pairs`, whose pairs
    are given: Salmonella as both sets with 5 pairs of different genes (4 ok, 1 saturated) got it. With `--pairs`
    it is now left out, and the lines that do pair a gene with itself are counted instead ("N of the M pairs of
    --pairs ... pair a gene with itself ...: their rows have dN = dS = 0"); without `--pairs` it is as before. The
    same `--cds` file given twice (by the same name or through a link) is refused up front ("--cds names X
    twice ... give it once", exit 1), instead of every record failing with "the ID occurs in two input files, X
    and X".
    No number changes (R01-R04): Salmonella / Shigella and S. cerevisiae / C. glabrata pairs and the four-yeast
    clusters run give byte-identical `dnds.tsv`, `orthologs.tsv`, `genes.tsv`, `failed.tsv`, clusters and
    alignments to the deep test's runs.
  - **Second retest of 2026-11-14 (Q01-Q05).** **The genetic-code hint never names a single wrong code on the
    test sets** (Q01, Q03). *C. elegans* mitochondrial genes (code 5) got "code 3" ("Run again with -c 3" for the
    12 genes alone): they avoid CGA / CGC as yeast mitochondrial genes do, and the start-codon check ran only when
    several codes tied, but 9 of the 12 start with ATT / TTG, which code 3 does not allow. A code that wins alone
    now passes the same start-codon test (at least two and a quarter of the genes that another candidate's starts
    fit); if not, it and the codes the codons do not argue against are the candidates and the start codons choose
    among them (4 and 5), and when only it would be left there is no hint. The one-in-ten start-codon bar of
    whole sets (R01) is used only for genome-size sets (a molecule of 150+ records, or 50+ evidence genes): two odd
    start codons of *Ciona*'s 9 mitochondrial genes had dropped its code 13. **`pipeline clusters` judges the
    code per `--cds` file** (Q02): over all files together a *Mycoplasma* with *E. coli* K-12 (`-c 11`) or
    *Tetrahymena* with human and mouse (`-c 1`) got no hint, and yeast + human mitochondrial genes a list without
    3 and 2; now one hint per file, naming it (`cds FILE: ...`; `file` in `report.json`). The molecule rule
    counts a molecule's records with one or more internal stops (not only those with two or more). **A byte copy
    of `--query` as `--subject` gets the same-file warning** (Q04; the files are compared by size, then SHA-256,
    a `.gz` by its content). **The same `.gz` twice in `--cds` is refused before it is decompressed**, and a
    run refused after decompressing (IDs in two `--cds` files) removes the `work/inputs` it made (Q05).
    The warning of `coal pipeline pairs` for the query set before (a007063) and after, on 108 runs: the 46 of
    the R01 evaluation and the retest's 62 (19 more NCBI genomes and their organelle subsets, BRAKER-style IDs,
    *Mycoplasma* + K-12 mixtures); 104 are unchanged:

    | Set | `-c` | Right code | Before | After |
    |---|---|---|---|---|
    | *C. elegans* RefSeq (30,560) | 1 | 5 (mito) | "10 of the 1,559 genes": code 3 | codes 4 or 5 |
    | *C. elegans* mitochondrial genes alone (12) | 1 | 5 | whole set: code 3, "Run again with -c 3" | whole set: codes 4 or 5 |
    | *Ciona intestinalis* mitochondrial genes alone (13) | 1 | 13 | whole set: codes 4 or 5 | whole set: codes 4, 5 or 13 |
    | *M. genitalium* + 229 K-12 records in one file | 11 | 4 | none | "253 of the 376 genes": codes 3, 4, 5, 9, 10, 13 or 21 |

    | Outcome (108 runs) | Before | After |
    |---|---|---|
    | a single code, the right one | 28 | 28 |
    | candidates that include the right code | 21 | 25 |
    | candidates without the right code | 1 | 0 |
    | a single wrong code | 2 | 0 |
    | no hint | 56 | 55 |

    With the right `-c` all 22 runs that have one are silent, before and after (*M. genitalium*, *M. pneumoniae*,
    *Mesoplasma florum*, *Ureaplasma parvum*, *Spiroplasma citri* `-c 4`; *T. thermophila*, *P. tetraurelia*,
    *Ichthyophthirius*, *Pseudocohnilembus* `-c 6`; *Blepharisma* `-c 4`; *Stentor* `-c 1`; the mitochondrial
    genes of *S. cerevisiae* and *D. hansenii* `-c 3`, *S. pombe* `-c 4`, *C. elegans*, *Apis* and
    *D. melanogaster* `-c 5`, *Ciona* `-c 13`, *S. purpuratus* and *S. mansoni* `-c 9`, human `-c 2`;
    *A. thaliana* plastid `-c 11`). `pipeline clusters`, before / after: *M. genitalium* + K-12 (+ *Salmonella*)
    `-c 11`: none / code 4 for the *Mycoplasma* file; *Tetrahymena* + human + mouse Ensembl `-c 1`: none / code 6,
    code 2, code 2; *S. cerevisiae* + human Ensembl: codes 4, 5, 9, 10, 13, 21 / code 3 and code 2;
    *S. cerevisiae* + *D. melanogaster* RefSeq: the same list / code 3 and codes 4 or 5; *S. cerevisiae* + human
    and *C. elegans* mitochondrial genes: code 3 (39 genes, the human and nematode ones among them) / code 3,
    code 2, codes 4 or 5; *Tetrahymena* + *S. cerevisiae*: whole set code 6 / code 6 and code 3.
    No number changes (Q01-Q05): Salmonella / Shigella and S. cerevisiae / C. glabrata pairs and the four-yeast
    clusters run give byte-identical `dnds.tsv`, `orthologs.tsv`, `genes.tsv`, `failed.tsv`, clusters and
    alignments to the deep test's runs.
  - **Final retest of 2026-11-14 (P02, P04, P05).** **Codes COAL does not support** (P02): *Ca.* Gracilibacteria
    and *Ca.* Absconditicoccus (code 25, TGA = Gly) read with `-c 11` got "codes 4 and 10 ... Run again with the
    -c of its organism", a list without their code, which `-c` does not take. A whole set's hint whose
    candidates all depart from the standard stops now also names the NCBI codes 24-33 that fit as well ("; or
    NCBI code 25 (Candidate Division SR1 and Gracilibacteria: TGA = Gly), which COAL does not support (-c takes
    1-6, 9-16 and 21-23)"), says that no `-c` gives right rows for such a code, and `report.json` lists them as
    `unsupported_candidates`; docs/pipeline.md lists the codes COAL does not support. **A copy of `--query` that
    is not a byte copy** (P04): a re-wrapped copy, or one with a record more, ran into 300 of 300 self-pairs
    without a word; without `--pairs`, 90% or more of the smaller file's records with the same ID and sequence in
    both (white space and case aside) now give the copy warning ("has the same genes as --query: 300 of the 300
    records ..."). Cheap: the first 200 IDs of `--subject` are looked up among those of `--query` (one header scan
    through mmap, 0.04 s for two 5 MB genomes) before anything is hashed (0.2 s for a 5 MB copy). **A missing
    input is refused before any `.gz` input is decompressed** (P05: the copies were made first and stayed in
    `work/inputs`); the error is as before.
    No number changes (P01-P05): Salmonella / Shigella and S. cerevisiae / C. glabrata pairs and the four-yeast
    clusters run give byte-identical `dnds.tsv`, `orthologs.tsv`, `unpaired.tsv`, `genes.tsv`, `failed.tsv`,
    clusters and alignments to the deep test's runs.
  - **Integration pass with orthologr of 2026-11-14 (C01-C03, C06, C08).** **`--preset pwalign` now returns
    pwalign's protein alignment, not only one with its score** (C01). Most pairs have several alignments with the
    best score; COAL's traceback put the gaps leftmost, pwalign (which fills its matrices over the reversed
    sequences) does otherwise, so orthologr's `engine = "coal"` and its native engine aligned about half of the
    pairs differently and their dN/dS differed by more than 1% in 5-27% of the pairs. New option `--tie-break
    leftmost|pwalign` (`coal pipeline`, `coal codon-dnds`, `coal-cpp-align`; `AlignParams.tie_break` in C++ and
    Python, `tie_break` in coalr): `leftmost` is COAL's choice and stays the default; `pwalign`, the default of
    `--preset pwalign`, is the mirror image of the same kernel run on the reversed pair with the sequences
    exchanged (no measurable cost: 4,164 yeast pairs 2.63 / 2.57 s, 945 human / mouse pairs 6.78 / 6.51 s on one
    thread). The preset also uses pwalign's copy of BLOSUM62 (new matrix `BLOSUM62-PWALIGN`: X scores -1
    against every residue, other B and Z values), so that scores agree for proteins with X too. `report.json`
    records `parameters.tie_break`. orthologr's native protein alignments (pwalign 1.6.0) that COAL reproduces
    row for row, before / after, and orthologr native vs COAL `--preset pwalign` (YN dS of the pairs where both
    give one, equal to the 6 digits orthologr prints):

    | Genome pair | Pairs | Identical protein alignments | YN dS equal to print precision | YN dS differing by > 1% |
    |---|---|---|---|---|
    | *A. thaliana* / *A. lyrata* | 995 | 505 → 995 / 995 | 677 / 994 → 994 / 994 | 267 → 0 |
    | *S. cerevisiae* / *C. glabrata* | 4,164 | 998 → 4,164 / 4,164 | 322 / 565 → 587 / 587 | 206 → 0 |
    | *E. coli* | 3,663 | 3,441 → 3,663 / 3,663 | 3,585 / 3,640 → 3,639 / 3,639 | 47 → 0 |
    | *Salmonella* / *Shigella* | 2,820 | 2,269 → 2,820 / 2,820 | 2,414 / 2,676 → 2,673 / 2,673 | 214 → 0 |
    | human / mouse | 946 | 397 → 945 / 945 | 484 / 924 → 926 / 926 | 256 → 0 |

    LWL agrees in the same way, Comeron to the 4 digits orthologr prints; what still differs is not the alignment:
    pairs whose YN iteration does not converge (COAL: `not_converged`, no estimate; KaKs_Calculator: its last
    round) and orthologr's NA rule for Comeron (docs/kaks.md). The pwalign fixture (tests/fixtures/pwalign, 85
    synthetic pairs and pwalign's alignments of them): 85 identical, 48 with `leftmost`. **A pair is computed in
    input order in every route** (C02): `coal pipeline clusters` (and `coalr::coal_dnds_clusters()`) aligned and
    estimated a pair aligned on its own with the representative first, the pairs route and orthologr with the
    gene that comes first in the input (orthologr's `dNdS_orthogroups()`: the first species' gene). Which gene
    is first decides among co-optimal alignments, and Comeron's estimate (libsequence counts the first sequence's
    sites). Now the gene of the earlier `--cds` file (then the earlier record) is first; the rows still list the
    representative first. 914 two-member clusters of four yeasts against the pairs route with the same pairs:
    differing rows YN 5 → 0, LWL 5 → 0, Comeron 809 → 0 (`--preset pwalign --comeron-flavour gestimator`); with
    the defaults 27 YN, 27 LWL and 210 Comeron rows of the 11,091 change. Clusters aligned by Clustal Omega / MAFFT
    are not affected (their alignment depends on the input order in another way; docs/pipeline.md, "Which gene
    comes first"). **Several
    `-m`** (C03): the summary line, the saturation notice and `report.json` described the first method of the
    rows only, unnamed (yeast with `-m YN -m Comeron -m LWL`: "46% saturated", LWL's share, where YN's is 86%), and
    failed.tsv listed only its error rows. Now the summary line and the notice name every method, `counts.pairs`
    names its method and `counts.pairs_by_method` has all of them, and failed.tsv has the error rows of every
    method (once per pair when all methods failed alike). **A cluster whose pairs have no estimate is
    `saturated`, `no_data` or `not_converged` in cluster_summary.tsv**, not `failed` (C06; `failed` when a pair
    failed), with a note that counts its pairs: orthologr's rule (four yeasts, three methods: 2,728 cluster rows were
    `failed` that orthologr calls `saturated` or `not_converged`).
    **`report.json` of a run from Python or R** (C08) has the equivalent `coal pipeline` command line in
    `command` (it was null) and `entry` (`cli`, `python-api`, or `coalr <version> coal_pipeline_pairs()`).
  - **pwalign's copies of every matrix: `--preset pwalign` with BLOSUM45, BLOSUM80 or PAM250 gives pwalign's
    alignments too.** pwalign's matrices (`data(<NAME>, package = "pwalign")`, pwalign 1.6.0) are not NCBI's:
    its BLOSUM45 scores B, Z, X and J otherwise (X -1 against every residue), its BLOSUM80 is the matrix on another
    scale (A-A 5, C-C 9, `*` -6 where NCBI's has 7, 13 and -8), and its BLOSUM62 has a J row (I or L) where
    `BLOSUM62-PWALIGN` read J as X. New matrices `BLOSUM45-PWALIGN` and `BLOSUM80-PWALIGN`, and J in
    `BLOSUM62-PWALIGN`, in C++, Python and coalr (the matrices now have 25 letters: NCBI's 24, then J, which
    NCBI's matrices read as X, as before). pwalign's PAM250 is NCBI's, value for value (and has no J): there is no
    PAM250 copy. **`--preset pwalign` now aligns with pwalign's copy of the matrix named** (`--matrix BLOSUM45`:
    `BLOSUM45-PWALIGN`; `--matrix BLOSUM62`: `BLOSUM62-PWALIGN`, where it was NCBI's BLOSUM62; PAM250 stays
    PAM250) in `coal pipeline`, `AlignParams.pwalign(matrix)` (C++ and Python, new argument) and `coal-cpp-align`
    (`-m` after `--preset pwalign`); `report.json` records the matrix used. A *-PWALIGN name given explicitly is
    used as given, with or without the preset; without the preset a plain name is NCBI's matrix, as before. New
    `pwalign_matrix_name()` (C++) / `coal.align.pwalign_matrix()` (Python): the name the preset uses. 104 real
    protein pairs (yeast, *Arabidopsis*, human / mouse, *E. coli*, *Salmonella* / *Shigella*, yeast mitochondria;
    some with B, Z, X, J, U, O, `*`), tests/fixtures/pwalign_matrices: identical alignments and scores to
    `pwalign::pairwiseAlignment(type = "global", substitutionMatrix = M)` for BLOSUM62, BLOSUM45, BLOSUM80 and
    PAM250, 104 / 104 each (before: 91, 83, 0 and 104; pwalign refuses U and O, and J with PAM250, so for the
    pairs with those letters the comparison is with pwalign's alignment of the pair with them as X, which is how
    COAL reads them). Other presets and explicit NCBI matrices give the same output as before, byte for byte.
  - **Pipeline messages and the order of MSA pairs (tutorial writers and testers, 2026-11-14).** **The pairs of a
    cluster aligned by Clustal Omega / MAFFT are computed in input order** (the gene of the earlier `--cds` file,
    then the earlier record, first), as the pairwise-aligned pairs already were; the rows still list the clusters
    file's order. New `O <rank> ...` line of `coal-cpp codon-dnds` (C++ `Group::rank`): every pair of a group of
    any size is computed with the member of the lower rank first. Only estimates that depend on which gene is
    first change: Comeron with libsequence's two-position branches (`--comeron-twosubs libsequence`, the default
    of `--comeron-flavour gestimator`). Four yeasts, 1,500 CDS each, DeepClust's clusters (443 MSA clusters,
    2,783 pairs, 1,548 of them listed later gene first), `-m YN NG LWL Li Comeron --comeron-flavour gestimator
    --li-flavour seqinr`: 1,459 Comeron rows change with Clustal Omega (1,455 with MAFFT; dN by 0.2% in the
    median, 1.3% at the 90th percentile; 12 change status) and 139 (143) Comeron rows of cluster_summary.tsv;
    Comeron rows equal to the same clusters aligned pairwise 174 → 240 (191 → 272), as many as Li's; no YN, NG,
    LWL or Li row changes. With the default options `dnds.tsv`, `cluster_summary.tsv`, `genes.tsv` and
    `failed.tsv` are byte for byte the same (both aligners). coalr's `coal_dnds_clusters()` does the same for
    clusters with `msa` rows or an `aligner` (the earlier gene of `cds` first). **`--resume` reuses the stages of a renamed or moved
    output folder** (or a copy that keeps the files' times): the signatures named every file by its full path,
    also those inside the output folder, so a renamed folder started again. Files inside the folder are now named
    relative to it; inputs keep their full path (the first `--resume` of a folder made before this change redoes
    its stages once). **The warning for several CDS records per gene names alternative start sites**: it said
    "several transcripts per gene" also for bacterial genes that NCBI annotates with two or three CDS of one locus
    tag (*E. coli* mrcB, infB), and now says "the input has several CDS records per gene (isoforms or
    alternative start sites; ...)", mentions the reciprocal best hits in pair mode only and counts "1 gene" in
    the singular. **docs/pipeline.md** says that a cluster without pairs (singleton, or fewer than two usable
    members) is `ok` in cluster_summary.tsv with `n_pairs` 0 (nothing failed; no estimate: filter on `n_ok` > 0;
    the value is unchanged), and links tutorials 2 and 3. **coalr: `preset = "pwalign"` takes pwalign's copy of
    the matrix named** (`matrix = "BLOSUM45"` → `BLOSUM45-PWALIGN`), as `coal pipeline` does.
  - **Review of 2026-11-15.** **A file in `-o DIR` that COAL did not write is not replaced without
    `--overwrite`.** When a file the run would write (`genes.tsv`, `dnds.tsv`, `failed.tsv`, `report.json`,
    `orthologs.tsv`, `unpaired.tsv`, `clusters.tsv`, `cluster_summary.tsv`, `realign.tsv`, and the optional
    outputs asked for) exists in the folder and no COAL run recorded it there (`work/outputs.json`; without it the
    outputs of `DIR/report.json`, or of the `--report` file when that is a COAL report of this folder), the run
    stops with exit 1, like the other refusals before a run starts, before anything changes and names the files;
    before, a folder with the user's own `genes.tsv`, `dnds.tsv` and `report.json` was overwritten with exit 0 and
    no word, and `-o examples/dnds` replaced the repository's `orthologs.tsv`. New option `--overwrite`
    (`overwrite=True` from Python, `overwrite = TRUE` in coalr) replaces them as before; it does not cover a
    foreign `codon_alignments/` / `protein_alignments/` folder, which stays refused. Runs into a folder of an
    earlier COAL run, `--resume`, `--report DIR/report.json` (the default spelled out), and
    `--pairs DIR/orthologs.tsv` / `--clusters DIR/clusters.tsv` of an earlier run (now listed under
    `read_in_place` in `work/outputs.json`, kept by the runs that neither read nor write the file, so a later run
    may write it again) behave as before. The check is made when the run starts. If `work/` is deleted after such
    a run, COAL goes by `DIR/report.json`, whose outputs list the `orthologs.tsv` it read, so a later DIAMOND run
    into `DIR` replaces a pairs file of your own there (docs/pipeline.md; keep your own file outside `DIR`).
    **`--diamond-args '-c 1 -b 4'` works**: a value of `--diamond-args`, `--deepclust-args` or `--realign-args`
    that starts with one of COAL's own short options (`-c`, `-o`, ...) failed with "expected one argument"; the
    word after the option (or an abbreviation argparse accepts, `--diamond-arg`) is now its value. Usage errors
    (exit 2) that show the `--diamond-args="-c 1 -b 4"` form: no value, a value that is one of COAL's options or
    an abbreviation of one (`--resum`), an unquoted string (`--diamond-args -b 4`), and the option given twice
    (argparse kept the last string without a word). Double quotes that the shell left in the words (cmd.exe
    `\"-c 1\"` gives `"-c` and `1"`) are not doubled in the quoted form shown: `--diamond-args="-c 1"`, not
    `--diamond-args=""-c 1""` (follow-up audit 2026-11-15). report.json's `command` is the command line as typed.
    `$COAL_<TOOL> is set to '...'` shows the path as set (Windows: no doubled backslashes).
- **Protein alignment for the codon step (docs/alignment.md).** `<coal/align.hpp>` and the `coal-cpp-align` tool: a
  global Needleman-Wunsch/Gotoh aligner (BLOSUM62/45/80, PAM250, affine gaps; EMBOSS needle defaults, a pwalign
  preset), full length (overhangs gapped, never cut off), AVX2 row loops picked at run time, multithreaded; optimal
  scores equal Biopython's and EMBOSS needle's. The expansion of DIAMOND `qseq_gapped`/`sseq_gapped` hits to
  full-length alignments without realigning. Python: `coal.align` (the same, with a pure-Python fallback) and
  `coal.msa` (Clustal Omega / MAFFT for clusters: find, install into COAL's cache, run, batch; stop codons handled).
  Clustal Omega (GPL-2.0-or-later) is installed on demand, not bundled (bundling would be licence-compatible now
  that COAL is GPL-3.0-only; it is not done). Benchmark on orthologr's 1,000-gene set:
  `benchmarks/results/align_layer_linux_20261111/`. `coal::core` now links `Threads::Threads`.

  Settled before the release (states that no release ever had are not listed as fixes):

  - **The platform wheels ship `coal-cpp-align` next to `coal-cpp`** (`coal/bin/` inside the package), so
    `coal.align.align_pairs` runs the native aligner after `pip install coal-align` on Linux, macOS and Windows
    instead of the pure-Python kernel (about 1,000 times slower). `packaging/wheels/build_platform_wheel.py` takes
    `--align-binary`; `build_linux_binary.sh` and the release workflow build and check both binaries;
    `packaging/wheels/check_wheel_binaries.py` runs in the environment the wheel was installed into (both lookups
    return the wheel's binaries, one alignment, one `coal pipeline pairs` run with an empty `PATH`).
    `coal-cpp-align --version` is new. Verified on Linux x86-64 with the repository's build script outside the
    container; the macOS and Windows steps of the workflow are changed the same way and have not run yet.
  - **An interrupt can no longer leave a Clustal Omega / MAFFT process behind.** `coal.msa` started an aligner
    and registered it a moment later; an interrupt in between left it running while the run reported that nothing
    was. Starting and registering now happen under the lock that `kill_running()` takes, and after
    `kill_running()` the interrupted batch starts no further aligner (its clusters get the error `interrupted`);
    the next batch runs normally.
  - **Very long pairs no longer take 600 MB or more** (deep test of 2026-11-14). The traceback is held at most
    256 MiB at a time (`AlignParams::traceback_bytes`); a larger pair is aligned with checkpoints, recomputing
    one block of rows at a time with the same code, so its alignment is identical. Human / mouse Ensembl pipeline
    (17,118 pairs with titin): peak 711 MB before, 227 MB after, 49.9 s → 53.2 s on 2 cores; `dnds.tsv` and
    the codon alignments byte-identical, and identical to the deep test's run. Forced on every pair (test hook,
    blocks of 1 and more rows): Salmonella / Shigella and S. cerevisiae / C. glabrata `dnds.tsv`, codon and
    protein alignments byte-identical.
- **R package `coalr` (packaging/r/coalr, docs/r.md).** COAL's C++ core compiled into an R package with Rcpp:
  `coal_dnds_pairs()` (CDS of ortholog pairs to dN, dS, dN/dS: the columns of `dnds.tsv`), `coal_dnds_clusters()`
  (every pair of each cluster and the per-cluster summary of `cluster_summary.tsv`; protein alignments from given
  rows, an aligner function such as `coal_msa_aligner("clustalo")`, or the built-in pairwise aligner),
  `coal_kaks()`, `coal_codon_align()` (PAL2NAL's job), `coal_align_pairs()` and `coal_translate()`. Everything runs
  inside the R process: no program is started, no file is written, worker threads never touch R, a user interrupt
  stops the run within about 0.1 s, and bad data is a row with a status, not an error. `coal_pipeline_pairs()` /
  `coal_pipeline_clusters()` run the whole `coal pipeline` through reticulate when Python with `coal-align` is
  there. The package holds a copy of `cpp/src` (`scripts/sync_r_sources.py`; `tests/test_r_bindings.py` checks it
  is current), so `remotes::install_github("Plumat-COAL/COAL", subdir = "packaging/r/coalr")` works. On
  orthologr's 1,000-gene example its 999 pairs equal `coal pipeline pairs` in every column and take 1.5 s on one
  thread. `R CMD check --as-cran`: no error or warning. Built and tested on Linux only so far.

  Settled before the release (states that no release ever had are not listed as fixes):

  - Version 0.6.0 (0.6.0.9000 in the development tree); `min_codons` / `max_ds` in `coal_dnds_pairs()`,
    `coal_dnds_clusters()` and `coal_kaks()`; its messages name COAL 0.6 and the release wheel (not PyPI).
  - **coalr: a path-like string that names no file is the error "file not found"**, not a sequence (docs/r.md).
  - **coalr: warnings where a call silently did the wrong thing** (new-user audit 2026-11-14, U02, U18, U22).
    `coal_dnds_pairs()` without `pairs` pairs the i-th query with the i-th subject, as before, and now warns when
    the IDs of such pairs differ, with how many of them have a protein identity below 20% (`pairs = "position"`
    asks for that pairing without the warning). When every pair fails, one warning names the common cause (the
    columns of `pairs` are not the IDs or are the other way round, protein given as CDS). A mistyped upper-case
    file name (`"CDS.FASTA"`), which is read as a sequence, gets an error that ends "it does not exist (looked in
    <working folder>)". `R CMD check --as-cran`: 0 errors, 0 warnings; the tests use at most two threads unless
    `NOT_CRAN=true`. `coal_pipeline_pairs()` / `coal_pipeline_clusters()` without COAL's Python pipeline name,
    in a development coalr (0.6.0.9000), the build from a COAL checkout instead of the latest release.
  - **coalr: `coal_info()$features`** names behaviours that a development build's version number does not tell,
    so that orthologr can recognise an older coalr of the same version (0.6.0.9000) and warn:
    `"comeron_twosubs"`, `"comeron_min_sites"` and `"comeron_min_sites_positive"` (the D9 refinement of #76: the
    min-sites rule only where it gives a distance above 0). Names are only added; a coalr without `features`
    predates all three. No number changes.
- **`coal --batch --batch-parallel N`: N jobs at a time.** Worker processes each import COAL once and take jobs
  from a queue (`0` = one per CPU core; the default is still one job after the other), every job's output
  byte-identical to running it alone; a job's messages are printed together when it ends, so the lines of two
  jobs never interleave, and the batch report keeps job order. Ctrl-C and SIGTERM end every worker and remove
  half-written outputs; a worker killed from outside fails its job and the batch goes on. coal-cpp refuses the
  option like the other batch options. 64 jobs of 2,000 × 300 (16-core WSL 2): 4.15 s one after the other, 1.20 s
  with 4 workers, 0.71 s with 16; tiny jobs gain nothing (docs/usage.md).

  Settled before the release (states that no release ever had are not listed as fixes):

  - **`--batch-parallel`: a worker killed while it sends a result no longer hangs the batch.** Jobs and results
    went through two queues shared by all workers; a worker killed in the middle of writing a result died holding
    the queue's lock, with half a message in the pipe, and the batch waited for ever. Every worker now has a pipe
    of its own, the batch process hands out the jobs and watches the workers themselves, and a killed worker fails
    its job with the documented one line (`error: the job's worker process died (out of memory?)`) while the
    batch goes on; the dead job's lock file is removed. Nothing else about `--batch` changes.
- **What to cite for what you ran.** `coal cite` prints a table (codon alignment: PAL2NAL; each Ka/Ks method: its
  estimator and the program COAL's code is derived from, KaKs_Calculator 2.0, PAML, libsequence or seqinr; the
  pipeline: DIAMOND, DeepClust, Clustal Omega, MAFFT) and the references; `coal cite YN` those of a method;
  `coal cite out/report.json` those of one pipeline run. The same table is in README.md, docs/kaks.md,
  `CITATION.cff` (references) and on the website. Every cited work has its DOI (second review, 2026-11-13: each
  reference was checked against the publisher's or the institution's page): `coal cite` prints it as a
  `https://doi.org/` link, `CITATION.cff` has `doi:` fields, the documents and the website link it. COAL itself
  still has no DOI. **DIAMOND DeepClust is cited by its paper**, Buchfink, Barbé, Ashkenazy, Reuter, Kennedy &
  Drost (2026), Nature Methods 23:724-727 (key `buchfink2026`); the key of its preprint, `buchfink2023`, is no
  longer listed but still resolves (`coal cite buchfink2023`, and a report that names it), followed by the paper.
  The seqinr reference has the book's full title.
  On the website's Cite page the table wraps and, on a phone, stacks the two cells of each row (it showed only
  its first column at 390 px).
- **A documented Python API for Ka/Ks and the pipeline.** docs/python-api.md has a new section "dN/dS and the
  pipeline": `coal.dnds_pipeline.run_pairs` / `run_clusters` (signatures, `PipelineResult`, the exit codes, and
  that a failed run returns exit code 1 with `report["error"]` instead of raising; `TypeError` / `ValueError` for
  a wrong call), `codon_dnds`, the column constants `PAIR_COLUMNS` / `SUMMARY_COLUMNS` / `FAILED_COLUMNS`,
  `METHODS`, `SCHEMA`, `find_diamond` / `install_diamond`, `coal.citations`, and the entry points of `coal.align`
  and `coal.msa`, with an example on `examples/dnds` that the tests run. The package ships `py.typed` and the
  public functions have accurate type hints. The versioning promise now covers the column names and order of
  `dnds.tsv`, `cluster_summary.tsv` and `failed.tsv` and the report schema id `coal.pipeline/1`. `import coal`
  exports what it did (no facade); its docstring no longer calls the rest of the package internal.
- **coalr: a vignette, and one name for the methods.** `vignette("from-two-cds-files-to-dnds")` goes from two
  CDS files to a filtered dN/dS table on an example shipped with the package (in process, in seconds, no
  network; built by `R CMD check`). `method` and `methods` are accepted by every function (a `method =` given
  to `coal_pipeline_pairs()` used to be ignored), `min_identity` is the R form of `--min-identity`, and `?coalr`
  and docs/r.md have one table of the command line's options and the R arguments, with the thread defaults.
- **`coal --help` lists every command, and `coal help [COMMAND]` works.** The help of `coal` (and of `coal-cpp`)
  ends with one line per command: `kaks`, `dnds`, `pipeline`, `translate`, `codon-dnds`, `cite`, `help`,
  `--cpp-path`. `coal kaks --help`, `coal dnds --help`, `coal translate --help` and `coal codon-dnds --help` print
  their help in an install without `coal-cpp` too. The "needs coal-cpp" message of the Ka/Ks commands and the
  pipeline says what the install is (the Python engine: codon alignment only) and names the platform wheel.
- **`examples/dnds`:** two small synthetic species with 12 known ortholog pairs (a range of dN/dS, all `ok`), a
  ready codon alignment, and the true pairs, for `coal pipeline`, `coal dnds` and `coal kaks`
  (`examples/make_dnds_example.py`; the older examples give only saturated pairs in `coal dnds`). README,
  docs/usage.md and the website use it.
- **GUI: a "Colour residues by class" setting.** Off, the alignment view and its downloadable page show the
  letters only (the mismatch outlines, frameshift hatching and the key's mark entries stay): for colour-blind
  readers who prefer no tints, and for printing. On by default; saved with the other settings.
- **App: Ka/Ks and dN/dS pipeline tabs, and the MCP tools `dnds_pairs` and `pipeline_pairs`.** The web GUI
  (`coal-gui`, `python -m coal_app`, the Windows app) now has four tabs: Align (unchanged), **Ka/Ks**, **dN/dS
  pipeline** and Settings.
  - **Ka/Ks** runs `coal dnds` on an uploaded or pasted aligned codon FASTA (every pair, or consecutive records)
    or a KaKs_Calculator AXT file: methods YN (default), NG, LWL, MLWL, LPB, MLPB, MYN, Comeron and Li (Comeron and
    Li: genetic code 1 only), any of COAL's genetic codes, and the optional flags "fewer codons than" /
    "dS above" (`--min-codons` / `--max-ds`: status `warning`, nothing dropped). It shows the counts per status,
    the table (first 300 rows; the TSV download has all), what each status in the run means, the equivalent
    command and what to cite; an example button loads `examples/dnds/codon.fasta`.
  - **dN/dS pipeline** runs `coal pipeline pairs` on two CDS FASTA files, one per species (gzip accepted):
    methods, genetic code, one transcript per gene (the longest), ortholog pairs by DIAMOND reciprocal best hits
    (`rbh`) or best hits (`bh`), the same flags, and optionally the codon alignments. It shows the four stages
    and the pipeline's progress lines while it runs, has a Stop button, and offers `dnds.tsv`, `report.json` and a
    zip of all result files. It runs only when `coal-cpp` and DIAMOND are found; otherwise the tab says how to
    install DIAMOND (`coal pipeline install diamond`) and "Check again" re-checks.
  - **MCP:** the stdio server (`python -m coal_app --mcp`, mcp 2.x and 1.x) and the GUI's MCP endpoint now have
    four tools: `codon_align`, `explain_report`, **`dnds_pairs`** (Ka/Ks for pairs of an aligned codon FASTA or an
    AXT file; arguments `pairs`, `methods`, `codontable`, `min_codons`, `max_ds`, `save_as`) and
    **`pipeline_pairs`** (two species' coding sequences -> ortholog pairs -> dN/dS per pair; `pairs_path` takes a
    table of given pairs, so DIAMOND is not needed, and `out_dir` a new or empty folder). Both return
    `completed`, the status counts, the first 20 rows, the rows that are not `ok` with their note, what each status
    means, `next_step`, the equivalent `coal` command, what to cite and the limits that applied; a problem the
    caller can fix before the run (missing file, unknown method or code, too large an input, no `coal-cpp` or
    DIAMOND) is a tool error. The server's instructions say to call `codon_align` first and then `dnds_pairs` on
    its `output_path` for dN/dS from a protein alignment and DNA. docs/mcp.md describes all four tools.
  - **Needs `coal-cpp`:** Ka/Ks and the pipeline (tabs and tools) need COAL's C++ engine, which the platform
    wheels and the Windows app include; an install with the Python engine only makes codon alignments, and the
    tabs and tools say so (in a development version such as 0.6.0.dev0 they name the build from the checkout,
    not the latest release, as `coal` does).
  - **Limits** (checked before the engine starts): Ka/Ks up to 100 MB and 500,000 results (pairs x methods) and
    one hour, in the GUI and over MCP; the pipeline up to 1,000 MB per file in the GUI (one run at a time, Stop
    button) and 500 MB per file and one hour over MCP. On a shared server (`COAL_HOSTED=1`, the Hugging Face
    Space): Ka/Ks up to 5 MB and 20,000 results, only the bundled examples, the run folder or uploads, no
    `save_as`; the pipeline is refused. The GUI's new tabs read only uploads, pasted text and the bundled examples
    (no typed paths). Results go to the run folder that is deleted after 24 hours (or `out_dir` / `save_as`); the
    pipeline's `work` folder is removed after every run, also after Stop. Larger runs belong on the command line
    (`coal dnds`, `coal pipeline pairs`); the pipeline's cluster mode and `coal kaks` stay command-line only.

- **First-use messages (new-user audit of 2026-11-14).** `coal pipeline`: a `--pairs` file none of whose pairs
  names a record of `--query` and of `--subject` (other columns, swapped columns, another ID type) stops with exit
  1 before any alignment, naming the IDs it looked for and the records' style, or saying that the columns look
  swapped (it ran to a table of `error` rows, exit 2); a `--clusters` file without one known member likewise; a
  first line of `--pairs` whose two fields are records of neither file, while later lines are, is read as a header
  (a notice says so); a CSV `--pairs` file gets "expected two tab-separated columns ... is this a CSV file?". A
  genome, chromosome or contigs given as coding sequences gets a warning (NCBI's `*_cds_from_genomic.fna` is what
  it needs). `--matrix` is checked when the command line is read (a typo came after the DIAMOND search, as a
  failure of `coal-cpp codon-dnds`). The first unusable record and its reason are on stderr, not only in
  failed.tsv. After Ctrl-C the line names `--resume`. `-o` onto a folder with an earlier run's results says that
  they are replaced. A Windows path on Linux / macOS names its WSL form (in the pipeline's messages and in
  those of `coal` itself, with either engine, and of `coal-cpp` run directly: "cannot open 'C:\Users\me\p.fa': no
  such file. This is a Windows path; under WSL (Windows Subsystem for Linux) the same file is
  /mnt/c/Users/me/p.fa"; the same words in both engines, compared by tests/test_engine_ux.py; one `/` between the
  parts and none at the end, so `C:\` is `/mnt/c` and `C:\\data\\x.fa` with doubled backslashes `/mnt/c/data/x.fa`,
  in all three: the pipeline's messages had a copy of the hint of their own that still gave `/mnt/c/` and
  `/mnt/c//data//x.fa` and now use the alignment's, follow-up audit 2026-11-15), and paths are shown as typed (no doubled backslashes). `-codontable N` is accepted as `-c N`; orthologr's and DIAMOND's option names
  (`--sensitivity-mode`, `--comp-cores`, `--threads`) and close typos of an option get a suggestion. `coal
  translate` warns when most genes have internal stop codons, with the first one and, as the pipeline does, the
  genetic code under which they have none. `coal dnds` / `coal kaks` on an input without any pair (an empty file,
  blank or comment lines only) exit 2 with one line, as docs/kaks.md's exit rule says (fuzz finding K05; it was 0).
  `coal-gui --check` reports coal-cpp and DIAMOND. In a development version (0.6.0.dev0) the coal-cpp-missing
  message names the build from the checkout, not the latest release (0.5.0 makes codon alignments only); README,
  docs/install.md and docs/tutorials/README.md say so at the top, and install.md's "from source" builds coal-cpp
  and says how to keep `COAL_CPP`.

### Changed

- **MAFFT and Clustal Omega are installed as one exact package per platform** (decision D6b of 2026-11-14). A
  version pin left the build open (conda-forge has three linux-64 builds of mafft 7.526, Bioconda twelve of
  clustalo 1.2.4). `coal pipeline install mafft|clustalo`, `--install-tools` and `coal.msa.install_*()` now ask
  for `channel::name=version=build` (`conda-forge::mafft=7.526=h4bc722e_0`, `bioconda::clustalo=1.2.4=h503566f_10`
  on linux-64; four platforms each) and keep the install only if conda's record of the package has the pinned
  SHA-256. The install prints the package and its SHA-256; the manifest and `report.json` (`tools.<tool>.package`)
  record it. A build the channel no longer serves falls back to the version pin with a warning. docs/alignment.md,
  "Exact builds". (MAFFT now comes from conda-forge, where it came from before too: conda-forge is the first
  channel and Bioconda has no 7.526.)
- **Comeron: Comeron's two pathways for codons that differ at two positions, by default** (decision D8 of
  2026-11-14). libsequence's `TwoSubs` (and so orthologr 0.4.3's `gestimator()`, and COAL until now) counts the
  branches codon 1 - intermediate 1, intermediate 1 - intermediate 2, codon 1 - codon 2 and codon 2 -
  intermediate 2 instead of the pathways codon 1 - intermediate 1 - codon 2 and codon 1 - intermediate 2 - codon
  2 (an argument-order slip; libsequence's repository is archived), which makes the estimate depend on which
  sequence is first. Swap test, 300 random 150-codon pairs: Ka changes in 279 and Ks in 190 pairs with
  libsequence's branches, in **0** with Comeron's pathways. New option `--comeron-twosubs comeron|libsequence`
  (`coal kaks`, `coal dnds`, `coal codon-dnds`, `coal pipeline`; coalr `comeron_twosubs`; C++
  `Options::comeron_twosubs`): `comeron` is the default of the default flavour, `libsequence` that of
  `--comeron-flavour gestimator`, which therefore still gives gestimator's numbers. `report.json` records the
  order used (`parameters.comeron_twosubs`). Effect on the default flavour (`coal pipeline pairs --preset pwalign -m
  Comeron`, new default against `--comeron-twosubs libsequence --comeron-min-sites 0`, which restore the previous
  output byte for byte):

  | Genome pair | Pairs | dS changed | dS median / p95 / max relative change | mean dS ratio new / old | dN/dS median / p95 change | status changes |
  |---|---:|---:|---|---:|---|---:|
  | *A. thaliana* / *A. lyrata* | 995 | 474 of 994 | 0% / 2.9% / 12% | 1.006 | 0% / 5.4% | 0 |
  | human / mouse | 946 | 809 of 932 | 1.0% / 4.8% / 50% | 1.015 | 2.1% / 8.2% | 6 |
  | *S. cerevisiae* / *C. glabrata* | 4,164 | 807 of 847 | 5.9% / 39% / 166% | 1.108 | 9.9% / 36% | 314 |

  (yeast: 308 pairs go from `ok` to `saturated`, 6 the other way; these are distant species, most of whose pairs
  are saturated either way.) In a neutral simulation (300 codons) the median error of dN against the truth is
  2.0% (libsequence's branches 2.3%), of dS -2.1% (-3.0%). docs/kaks.md, "Two-position codon differences".
- **Comeron: a site class with fewer than 3 sites cannot make dN or dS `NA`** (decision D9 of 2026-11-14). In
  short alignments the 2-fold transversional class has one or two sites, and a single change there made the whole
  estimate `NA`. Such a class's undefined logarithm is now set to 0 (gestimator's treatment), with status
  `warning` and the note `Ka: a site class with fewer than 3 sites has an undefined correction (its logarithm set
  to 0)`; no value that was not `NA` changes. New option `--comeron-min-sites N` (coalr `comeron_min_sites`;
  `0` restores the previous rule). Simulated neutral pairs, 600 per length: dN `NA` at 30 codons 22 (3.7%) to 0,
  the recovered values' mean absolute error 8.8% (all 30-codon pairs: 11.2%); 50 codons 2 to 0; 300 codons none
  before or after. Real data: *A. thaliana* 0 rows, human / mouse 3, yeast 14 get a value. In the gestimator
  flavour (whose numbers are gestimator's either way) only the note changes, so that orthologr can tell such a row
  from a saturated one. docs/kaks.md, "Small site classes".

  **Refined: the rule applies only where it gives a positive distance.** In very short alignments the small
  class's logarithm set to 0 could leave dS at exactly 0 although a synonymous change was seen (`ATGGCTAAA` /
  `ATGGCCAAA`: dS 0), which reads as "no synonymous change" (orthologr labels it `dS_zero`). Where a small class
  was the reason and the distance is not above 0, the row is now what `--comeron-min-sites 0` gives: `NA` with
  status `saturated` (default flavour); gestimator's 0 with the note `Ks: correction undefined, gestimator's
  value (undefined logarithms set to 0)` (gestimator flavour, the one orthologr's COAL engine uses, so that
  orthologr gives `NA`, not a plain 0). 7,000 simulated neutral pairs, 1,000 per length (an independent re-check
  of D9): dS was 0 this way in 53 pairs at 5 codons, 4 at 10, 1 at 15 and none at 20, 30, 50 or 100; 57 of the
  58 had true synonymous substitutions. All 58 are now `saturated` (gestimator flavour: noted); dN never became 0
  (none of 168 small-class values); no finite value changes, and the `ok` rows are the same (6,738). COAL's
  hand-made test cases: 17 rows (15 dS, 4 dN values), all with at most 4 codons compared. `coal pipeline pairs`
  on *A. thaliana* / *A. lyrata* (1,000 genes) and on yeast mitochondria writes byte-identical output files,
  also with `-m Comeron` in both flavours; in the real data of D9 none of the 17 small-class rows is 0 (re-check).
  The help of `--comeron-min-sites` (`coal dnds`, `coal kaks`, both engines) says so: a small class "cannot make
  dN or dS NA unless it is then 0".

- **Licence: GPL-3.0-only from the next release.** COAL (the Python and C++ engines, the `pal2nal` drop-in, the app
  and everything else in this repository) is now under the GNU General Public License, version 3 only (`LICENSE`,
  `NOTICE`; package metadata, `CITATION.cff`, the Bioconda recipe, the website and the Windows app say so). **0.5.0
  and every earlier release remain under the MIT License**: nothing already published changes. Why: the Ka/Ks
  (dN/dS) engine that the next release adds is derived from KaKs_Calculator 2.0 and from PAML's yn00, which are
  both GPL version 3 (without an "or later" clause), so COAL as a whole has to be distributed under GPL-3. All
  contributors agreed to the change. What it means for users: running COAL (the commands, the app, a pipeline step)
  does not affect the licence of your data, results, scripts or workflows. A program that includes COAL's code or
  links its library (the C++ `coal::core`, `import coal` in Python) and is distributed must be under a
  GPL-3-compatible licence; programs that only start `coal` / `coal-cpp` as a separate process are not affected.
  The third-party notices already in the repository (the website's font, the real-data test sets) are unchanged.
- **The development tree has the version `0.6.0.dev0`** (`coal --version`, `coal-cpp --version`, reports), so a
  build of it is never taken for the MIT release 0.5.0; coalr is 0.6.0.9000. The release sets 0.6.0
  (`scripts/release_version.py`). The run report's schema accepts such a version.
- **`coal --help` is longer:** after the options it lists the other commands (see "Added"); a script that
  compares the help text byte for byte sees a difference. The options and their help lines are unchanged.
- **Two runs writing the same `-o` file at the same time are refused** instead of sharing one `<out>.coal-partial`
  file, which could end a run with exit 0 and an empty or another run's output. A run with `-o` keeps
  `<out>.coal-lock` while it runs; a lock left by a killed run is replaced once its process is gone (a lock
  file that is still empty is one being created by another run, not a stale one: the second iteration of the
  stress test caught two runs reading each other's empty lock and both going ahead).

  Also part of the lock as released:

  - **The output lock (`<out>.coal-lock`; `coal -o` and `coal pipeline`) holds under simultaneous starts, through
    links and across machines** (the Python side). Several runs that started together and found a stale or empty
    lock could all take it; a stale lock is now replaced by one run at a time (`<out>.coal-lock.takeover`, held
    for an instant) and only if it still is the lock that was judged stale. The lock is next to the real output,
    so a symbolic link to a pipeline folder no longer has a lock of its own. The lock line is now
    `<pid> <host> <process start> <token>` (the pid stays first, which is all `coal-cpp` and earlier versions
    read): a lock from another machine counts as held and the message names the machine, a process id that is
    alive but started at another time (Linux) is not the holder, a second call into the same output from the same
    process is refused, and releasing never removes a lock that is not the run's own (docs/pipeline.md).
  - **A refused `--resume` (or a signal) no longer leaves `<out>.coal-lock` behind, both engines.** The lock was taken
    before the checkpoint set-up and only released by the block after it, so an error raised in between (the inputs
    changed since the checkpoint, a checkpoint of the other engine) or a signal there left the lock, and the next run
    reported a stale lock of a dead process. Found by the new `--batch-parallel` SIGTERM test on a macOS CI runner;
    a test refuses a resume in both engines and checks the lock is gone.
  - **A SIGTERM, SIGHUP or Ctrl-C at any instant leaves no lock and no `<lock>.<pid>.<hex>.tmp` (final retest of
    2026-11-14, P01).** The cause of the flaky `test_batch` SIGTERM test (2 of 84 CI jobs on PR #62): Python's
    handler raises SystemExit wherever the main thread is, and a signal right after the lock was linked, while its
    temporary file was removed, after `take_output_lock` returned and before the caller's `try`, or during the
    final release left the lock (and in one window its temporary file, which nothing ever removed). Now the
    creation and the release of a lock, and the clean-up at the end of a run, run with SIGINT / SIGTERM / SIGHUP
    blocked on POSIX (delivered when the step is done); the bookkeeping is ordered so that what a signal can
    still cut short (Windows, other threads) is cleaned up: `take_output_lock` removes the lock it made when
    anything ends it, a lock is forgotten only after its file is removed, and every lock and temporary file this
    process still holds is removed at exit and by a batch worker after each job. A batch worker also closes the
    job's report (its entries file `<out>.coal-entries.jsonl` stayed after a SIGTERM), and a checkpoint is removed
    after its files. coal-cpp registered the lock for its signal handler only after reading it back: a SIGTERM in
    between left it; the steps from creating the lock to registering it, and the release, now run with those
    signals blocked (POSIX). Numbers: the retest's deterministic windows (`coal pep nuc -o out.fa` sending itself
    SIGTERM at each instant) left the lock in 4 of 5, now in none; the retest's storm (the CI scenario, SIGTERM
    after a random 0-60 ms) left a lock, `.tmp` or entries file in 4 of 600 batches (retest) and 2 of 300 under
    load (here), now in none of 600 (300 unloaded, 300 under 3 CPU hogs); coal-cpp with SIGTERM right after the
    link: lock left, now nothing. `tests/test_batch.py` passes 100 of 100 times under 3 CPU hogs. The test
    itself now runs without checkpoints (`COAL_CHECKPOINT_EVERY`): a job that has run for 2 s when the signal
    comes keeps its checkpoint for `--resume` by design, which a loaded machine reached in 3 of 300 storm batches.
  - **`coal-cpp -o`: the output lock of the C++ engine follows the Python engine's rules.** `<out>.coal-lock`
    holds `<pid> <host> <start> <token>` in both engines; it is created with its content in one step (never seen
    empty), a lock written on another machine is held and never replaced, a process id that was reused (another
    start time; Linux) is not the holder, and a stale lock is replaced by one run at a time (`<lock>.takeover`),
    which removes it only if it is still the one that was judged stale. Before, several `coal-cpp -o` runs that
    found a stale lock at the same moment could each remove "the stale lock", the later ones removing the lock an
    earlier one had just taken, and both wrote the output; and a lock from another machine counted as stale.
    Tests: the race step by step (the unguarded code fails it), 200 trials of six simultaneous starts onto a stale
    lock (one holder each time), Python and C++ starts mixed, each engine refused by the other's lock. The
    Windows code path (process handle, `CreateHardLink`) is compiled by CI and not exercised by these tests there.
    **On Windows** the first CI run of these tests found the guard `<lock>.takeover` left behind after simultaneous
    starts (one holder each time, as required): Windows does not remove a file that another process has open, the
    runs that lost the race were reading the guard at that instant, and the removal was tried once. Both engines
    now try every removal of a lock, guard or temporary file again for about 2 s, remove the guard on every path,
    and clear a guard that is left all the same when its process is gone or it is older than 15 s (before: 30 s,
    the time a run waits). The host name in the lock has one documented source and spelling in both engines and is
    compared without regard to case on Windows (`docs/pipeline.md`). Tests on every system: a hook makes removals
    fail as they do on Windows.
- **Installers: what is in COAL's cache is checked before it is run.** Each installed binary gets a
  `<binary>.coal-manifest.json` with its SHA-256, and a cache binary is used only while it has that hash.
  **Tools installed into the cache by COAL 0.5.0 or earlier have no manifest and are no longer used
  automatically** (that includes a Clustal Omega fetched without verification): COAL says so once and
  `coal pipeline install <tool> --force` installs them again. Tools on `PATH` or given by `COAL_CLUSTALO`,
  `COAL_MAFFT`, `COAL_DIAMOND`, `--diamond` are unaffected. Also: downloads go to a new randomly named file
  instead of the fixed `<name>.part` (two installs collided there, and a link at that name was written
  through), the hash is taken from the file on disk, one install per tool runs at a time and the others wait, an
  archive that was cut short (possible with `--allow-unverified`) is a clear error instead of an `EOFError`,
  and a zip member name with `:`, an absolute path, `..` or a reserved device name is refused
  (docs/alignment.md).

  Also:

  - **Installers verify before they run anything.** Downloads are kept only if their SHA-256 is the pinned one:
    micromamba is one fixed release (2.10.0-0) with pinned hashes, not `latest`; the DIAMOND 2.2.8 archives are
    pinned. Clustal Omega is installed from Bioconda only; see the entry on the Clustal site below
    (docs/alignment.md).
  - **The download of Clustal Omega from the Clustal site is removed.** That site's address now serves an
    unrelated "people search" site and its three pre-compiled binaries (Linux, macOS, Windows) answer "not found"
    (checked from a browser on 2026-11-12). There is nothing to pin, and anything that appears there later would not
    be from the Clustal authors, so COAL no longer holds an address of that site: `install_clustalo(method=
    "official")` is refused with this explanation, and `sha256=` / `allow_unverified=True` (`--sha256`,
    `--allow-unverified`) cannot bring it back (with Clustal Omega they are an error). Linux and macOS install it
    from Bioconda (conda / mamba / the pinned micromamba) as before. **On Windows there is now no automatic Clustal
    Omega install**: the message names what works (WSL; MSYS2's package of Clustal Omega 1.2.4, `pacman -S
    mingw-w64-ucrt-x86_64-clustal-omega`, or another `clustalo.exe` you trust, on `PATH` or in `COAL_CLUSTALO`;
    `--aligner mafft` with MAFFT's own Windows package, the "all-in-one" `mafft-7.526-win64-signed.zip`, in
    `COAL_MAFFT`). Neither Bioconda nor conda-forge has a Windows package of Clustal Omega or MAFFT (checked
    2026-11-14), and COAL adds a download only with a SHA-256 computed from the file itself, which could not be
    done for Clustal Omega's, so it is not downloaded; MAFFT's package is now pinned by its own SHA-256 and
    installed by `coal pipeline install mafft` (see "Added"), and docs/install.md "Clustal Omega and MAFFT on
    Windows" has the steps. Bundling a self-built `clustalo.exe` is
    another way to support Windows and is not done yet (docs/alignment.md, "Licences").
- **`pal2nal` without arguments is a usage error:** the usage text goes to standard error with exit 1 (as the
  drop-in's help says and as PAL2NAL's own usage goes to stderr), not to standard output with exit 0, so a pipeline
  step that forgot its files no longer gets the usage text as its "alignment". `pal2nal -h` is unchanged (stdout, 0).
- **`-q` without `--report`/`--json` keeps no per-sequence entries** and writes none to the `--resume` entries
  file (a million genes of unrelated DNA wrote 55 GB of it), so memory and disk stay flat on data where every
  codon mismatches (the summary lines, the hint and the report, which read them, are not produced under `-q`).
- **The genetic-code hint no longer re-tests genes of more than 100,000 residues** (28 s of a 30 s run on one
  1,000,000-residue gene with one wrong codon; `-q` already skipped it).
- **Messages name at most three input files** ("a.fa, b.fa, c.fa and 997 more files") instead of every path.
- **coal-cpp: the exact-fit search is linear in the DNA** (a string search over the translated frames) instead of
  trying every start codon by codon: a 200,000-residue low-complexity gene behind a near-copy of itself took 138 s
  (an hour at 1,000,000 residues), now well under a second, as in the Python engine. Same answer everywhere
  (the differential tests and the engine fuzz).
- **GUI: the alignment view's residue-class colours are readable for colour-blind users.** The seven classes keep
  their Okabe-Ito hue families but now sit on an even lightness ladder, so they stay apart under protanopia,
  deuteranopia and tritanopia and in greyscale (the pastel tints of equal lightness merged there; verification
  2026-10-24, V-17); the text stays at 5:1 or better on every tint, in the light and the dark theme. A test holds
  the simulated distances. The letters, outlines and legend carried the meaning before and still do.

### Fixed

- **`coal`: a genetic code that NCBI merged into another (7 into 4, 8 into 1) names the code to use** instead of
  only "invalid choice: '7'"; and the summary of a run that stopped no longer adds "(kept in the alignment)" to
  "N of M sequences have codons that do not match their amino acid" next to "no output was written" (both
  engines; new-user audit 2026-11-14, U20, U21).
- **The summary of a `coal` run that stopped no longer says that a probably-wrong gene "was written to the
  output"** (both engines; found by the tutorial writers, docs/tutorials/01-codon-alignment.md step 5). With
  sequences paired by order, one gene flagged as probably wrong and another that stops the run (exit 1), the
  summary said both "1 sequence of 3 was written to the output but is probably wrong (... COAL still wrote its
  row ...)" and "no output was written"; no file was written. When no output was written the line is now "1
  sequence of 3 is probably wrong (first: <id>: <the first sentence of its message>) Fix or remove it before
  running again, or use --strict to refuse such genes."; with `--keep-going`, or a run that completes, it is as
  before, word for word. `report.json` `summary_text` follows.

In commands that 0.5.0 has (both engines unless said):

- **A helper process killed before it had started no longer ends the run with exit 141 and no message.** With
  the spawn and forkserver start methods (macOS, Windows, Python 3.14 on Linux) the first child of a run with
  `-j` is Python's resource tracker; if it was killed at once, the pool's first write to it failed with a broken
  pipe, which COAL took for a closed standard output. It is now the usual
  `coal: error: a worker process died (out of memory?)` (exit 1). (This was so before the review; `--batch-parallel`
  reports a worker that cannot be started the same way.)
- **Python engine: a worker killed while it sends a result no longer hangs the run.** The pool's management
  thread then waits inside a read that never completes and never reports the dead worker, so `coal -j N` waited
  for ever (seen as an intermittent CI failure on macOS with Python 3.10; reproduced on Linux with the spawn
  start method, about 1 run in 16). The main process now also watches the workers itself while it waits and ends
  with "a worker process died" (exit 1). Ending the broken pool no longer closes the call queue under the
  management thread, which on Python 3.10 / 3.11 printed a `ValueError: Queue ... is closed` traceback before
  that message (about 1 run in 5 there).
- **Thousands of DNA files:** at most about as many files as the open-file limit allows stay open (least recently
  used closed first, the limit raised towards its hard limit when needed) and their read buffers together take
  about 64 MB at most. 1,000 files failed with "too many open files" in the Python engine with workers, 60,000 open
  files took 8 GB (129 KB each on Python 3.14), and more than 65,535 files overflowed a 16-bit file number: a
  traceback in the Python engine, a silent wrap in coal-cpp.
- **A lone `\r` ends a line anywhere in a DNA file** (universal newlines, as the protein reader): the style used to
  be decided from the first 64 kB, so a file whose line ends change part-way lost every record after the first lone
  `\r` beyond that and failed with a misleading count error. PAL2NAL v14 reads such files.
- **`nohup coal ...` survives the terminal closing:** an already ignored SIGHUP stays ignored (the handler used to
  replace nohup's setting; the run died with exit 129).
- **Out of memory in the main process is one `coal: error: out of memory ...` line** (exit 1, report code
  `out_of_memory`) instead of a traceback (a worker dying of it was already handled).
- **Memory no longer grows with the number of per-codon warnings, in either engine** (finding 9 of the stress
  test, both halves). The run report keeps, per sequence, one byte for its kind, one for a few flags and the
  codon of its first definite mismatch (plus the first flagged and the first excluded gene whole, for the
  summary lines, and where up to 200 failed genes break), instead of every per-sequence entry with every message
  until the end of the run; the entries go to the `--resume` entries file next to the output as the run goes
  (or to a temporary file when a report is wanted without `-o`; nothing with `-q` and no report), and `--report`
  / `--json` read them back in order. A resumed run rebuilds its state from that file one line at a time (the
  checkpoint check no longer loads the whole file either). The Python engine's worker processes hand back each
  sequence's mismatch columns and codons instead of its warning strings (about 2 KB against 35 KB for a
  400-codon gene of unrelated DNA, per gene of the batch in flight in every worker and in the main process), and
  the main process formats a warning only when it prints it; the strings, the report, stderr, the exit status
  and `--resume` are unchanged byte for byte (`tests/golden_memory`: recorded from the previous `main`, both
  engines). Measured on 400-codon genes of unrelated DNA, every codon a warning (peak memory as the summed PSS
  over the process tree; 16 threads): 20,000 genes in default mode, coal-cpp 1,172 MB before and 103 MB after,
  the Python engine 2,822 MB and 45 MB (`-j 1`), 3,226 MB and 338 MB (`-j 16`, 7 workers); 100,000 genes,
  coal-cpp 5,654 MB and 118 MB, the Python engine 13,761 MB and 69 MB (`-j 1`), 14,647 MB and 885 MB (`-j 16`,
  16 workers: 545 MB of that is what the same input in clean form costs with 17 processes, the rest the batches
  in flight). The default mode is also faster on such
  data (100,000 genes, `-j 1`: 213 s to 158 s in Python, 185 s to 161 s in coal-cpp) since no entry store grows;
  the one cost is the main process formatting the strings the workers formatted in parallel before, which only
  shows when stderr takes most of a run. A run with `-o` on such data writes the entries file next to the
  output while it runs (about 55 KB per 400-codon gene of unrelated DNA, removed at the end).
- **App: `python -m coal_app` no longer dies on a non-UTF-8 console** (Gradio's MCP banner emoji on a cp1252
  Windows console or pipe).
- **GUI: a run that stops before pairing says why in the banner** (the engine's own sentence, e.g. "control character:
  sequence files must be plain text (is this a binary or word-processor file?)", file names only); it used to say
  only "COAL stopped: no alignment was written" with the reason inside the collapsed raw-messages box.
- **GUI: the genetic-code tip on a large input says it tested a sample** ("50 of the 50 genes COAL tested (a sample
  of the 100,000)") instead of "50 of 50 genes".
- **GUI accessibility:** the "Methods paragraph" box is announced by its label (its copy button, which sat inside the
  label and took the label for itself, is now a separate "Copy methods paragraph" button), and Gradio's hidden
  duplicate tab bar no longer holds focusable buttons (axe `aria-hidden-focus`).

### Internal

- **Engine review (2026-11-15).** The AVX2 copy of YN / MYN's transition-probability code is compiled on x86-64
  only (GCC / Clang, and only with SSE2 floating point, `__SSE2_MATH__`): on 32-bit x86 the baseline does its
  double arithmetic on the x87 unit by default, so the two copies could round differently; a 32-bit build has
  the baseline copy only. align.cpp's AVX2 row loops (integer DP) stay on 32-bit x86. `coal-cpp codon-dnds` (the
  native step of `coal pipeline`) accepts the testing flags `--kaks-baseline` and `--kaks-isa` of `dnds` /
  `kaks`, not in its help. The Ka/Ks property fuzz's P10 judges MLWL on LWL's scale (its Sd over L2 / 3 + L4):
  MLWL's kappa-based sites make a genuine dS of 15.6 (513 CAG / CAA pairs at kappa 1, seed 4243) or 162.7 (kappa
  0.5031) possible, as KaKs_Calculator prints. `benchmarks/kaks/random_short_pairs.py` classes difference 18.
  Unused imports removed (coal/codon.py, coal/align.py). New tests: `tests/test_engine_review_fixes.py`.

- **Review of 2026-11-15, M2: the reference fuzz's classifier knows differences 16, 18 and 19.**
  `tests/fuzz_kaks_vs_reference.py` reported about 190 documented differences as "unexplained": the LWL family
  with genetic codes other than the standard one (difference 16; now explained only where a codon column of the
  pair has a step that KaKs_Calculator's hard-coded standard-code rules and COAL's rules for the code count
  differently, checked by reimplementing both), MLWL at kappa exactly 2 or 0.5 (difference 18, with
  `random_short_pairs.py`'s rule: the sites must be the two formulas) and the counts of difference 19. It also
  compares the LWL family's S- and N-Substitutions now, and an LPB S-Sites (N-Sites) `NA` is explained only when
  Ks (Ka) is saturated or `NA`. New `tests/test_fuzz_classifier.py`: one fixture per difference of docs/kaks.md
  (the list is read from the docs, so a new difference fails until it has one), and fixtures that must stay
  unexplained. Seed 77115 x 500 pairs (the generator's codes) and seed 990017 x 1,500 (`--codes 1,4,11`): 0
  unexplained against the 15-digit and the corrected KaKs_Calculator and against yn00 (weighting 1 and 0).

- **Follow-up audit (2026-11-15), the low findings.** nightly.yml's `cpp-extras` root step passes `PATH`,
  `COAL_CPP_BIN`, `COAL_REQUIRE_CPP=1` and `COAL_IMPL=cpp` through `sudo env` by name: sudo-rs (Ubuntu 26) ignores
  `--preserve-env` ("preserving the entire environment is not supported"), so `tests.test_robustness.TestDiskFull`
  ran without them; that class now fails, instead of comparing nothing with the Python engine alone, when
  `COAL_REQUIRE_CPP` is set and coal-cpp is not found (tests/test_repo_ci.py pins the step).
  `tests/cpp_wrapper.py` says what the coal-cpp stand-in does not do (escape cmd.exe's metacharacters; on
  Windows, pass a kill on to the real coal-cpp), and tests/test_kaks_dispatch.py that its clang-cl detection away
  from a build tree has not been checked on a real clang-cl build. `benchmarks/kaks/random_short_pairs.py`
  explains an MLWL difference at kappa exactly 2 or 0.5 (difference 18) only when COAL's S-Sites / N-Sites are
  the site formula of the exact kappa and KaKs_Calculator's those of the formula on the other side (it accepted
  any difference on such a row); both random sets give 0 unexplained against the 15-digit and the corrected
  KaKs_Calculator (single: none of the class; multi: 36 fields in 10 pairs, all of them the other formula). The
  property fuzz's P10 checks MLWL's LWL-scale value (Sd / (L2 / 3 + L4), Nd / (L0 + 2 L2 / 3)) whatever MLWL's own
  dN / dS (it was looked at only above 10, so an LWL-scale value just above 10 with MLWL's own below passed);
  KNOWN stays empty, and 3 seeds x 1,500 pairs (all methods) and 3 x 5,000 MLWL-only pairs (all genetic codes)
  have no finding.

- **Bioconda recipe (packaging/bioconda), second review.** The test line for `pal2nal` without arguments expected
  exit 0 and the usage on stdout and would have failed the first build: it now tests what the command does (exit
  1, usage on stderr, nothing on stdout). `python` is unpinned in `host:` / `run:` and Pythons below 3.10 are
  skipped by a selector (the form for a compiled, non-noarch package); `linux-aarch64` and `osx-arm64` are
  additional platforms; the package is the command line and the Python package only (`build.sh` no longer
  installs the static library, the headers and the CMake package); the two test commands with tab characters are
  quoted YAML strings. Rendered with `conda-render`, `build.sh` and all 23 test commands run against a virtual
  environment; a real `conda build` and `bioconda-utils lint` are still to do (packaging/bioconda/README.md).
  `tests/test_bioconda_recipe.py`.

- **Package metadata follows PEP 639:** `license = "GPL-3.0-only"` (an SPDX expression) with
  `license-files = ["LICENSE", "NOTICE"]` in both `pyproject.toml` files, no licence classifier; the build needs
  setuptools 77.0.3 or later (Python 3.9+, below the floor of 3.10; the Bioconda recipe asks for it too). The
  wheels' metadata says `License-Expression: GPL-3.0-only` and carries both files under `licenses/`; the
  deprecation warnings of every build are gone.

Documentation, tests, continuous integration and packaging:

- **Continuous integration.** A workflow of its own for the R package (`.github/workflows/r-package.yml`): `R
  CMD check --as-cran --no-manual` of `packaging/r/coalr` on Linux (R release and the previous release), macOS
  and Windows, and the check that the vendored C++ sources are `cpp/src`; it runs only when the package, the C++
  core or the sync script change, so the main matrix is not slowed. The nightly run has a ThreadSanitizer job
  (`tests/tools/tsan_threads.py`: `coal-cpp kaks`, `dnds`, `codon-dnds -j 4` with all nine methods and both
  flavours, and `coal-cpp-align -j 4`, each also compared with `-j 1`; 28 commands, no report locally). The C++
  job requires `tests.test_kaks_comeron_li`, `tests.test_short_options` and `tests.test_output_lock` (a missing
  binary fails), checks the manylinux wheel with both binaries in a fresh environment, and installs and
  requires DIAMOND on macOS too: DIAMOND's macOS archive is a universal binary with a native arm64 slice.
- **`coal-cpp`: test hooks only in test builds.** The `COAL_TEST_*` environment variables with which the tests
  pause a run or line up several runs are compiled into `coal-cpp` only with the CMake option
  `-DCOAL_TEST_HOOKS=ON` (default `OFF`): release builds, wheels and packages do not contain them
  (`docs/development.md` "Test hooks").
- **Line ends in a Windows checkout** (`.gitattributes`): the R package (an R source package has LF line ends;
  its vendored C++ sources are compared byte for byte), the Ka/Ks goldens, the real-data cases, the stored
  benchmark results and shell scripts are LF in every checkout, as the fixtures and goldens already were.
  `scripts/sync_r_sources.py --check` reads `cpp/src` with either line end and names a copy that differs only by
  its line ends as that. (With `core.autocrlf`, as on the Windows CI machines, it called all 13 vendored files
  "out of date".)
- **Docs: Linux benchmark numbers.** `docs/benchmarks.md` has the 2026-10-24 run on the same desktop as the Windows
  table but under WSL 2 (up to 1,000,000 sequences, PAL2NAL v14 against both engines, raw data in
  `benchmarks/results/linux_wsl_4159e10/`): COAL's times and memory hold, PAL2NAL's Perl runs about 3.6× faster there
  than on Windows, so the speed-ups are 123× (coal-cpp) and 29× (Python) at 1M on that setup.
- **Docs:** the desktop README's privacy line says none of *your data* leaves the computer and names the one helper
  script Gradio's page requests from a CDN (not needed; the page works with it blocked); the installer name is the
  version-free `COAL-Setup-windows-x64.exe`; `docs/usage.md` says a `kill -9` of a run writing to standard output can
  leave `.coal_<random>.tmp` in the temp folder; `tests/README.md` says which of the "17 local codon alignments" are
  in the repository (12) and which are local-only (5); SKILL.md's "after the engine update" and "not yet verified on
  Windows" sections are current (both engines verified on Windows 2026-10-24); `demo/mcp_demo.md` regenerated
  (`codon_align` keeps going past a broken gene and returns `completed: true` with `exit_code: 2`).
- **Tests:** `test_docs_explainer` finds Playwright's Chromium on every platform (it looked for a Linux-only path and
  was always skipped elsewhere); `test_cpp_vs_reference` copies instead of symlinking where symlinks need a privilege
  (a stock Windows account); `test_resume_with_workers` ignores a venv launcher's same-argv child on Windows; new
  tests for the GUI's download buttons (the one gap found by mutation testing: 12 of 15 planted bugs were caught),
  the stop reason, the sampled tip and the methods label. The PAL2NAL fuzz harness no longer reports a difference
  for `-codontable 10` genes that both programs place by position 1 or 2 nt off the generator's frame (the
  neutralised TGA then straddles two codons of the alignment grid; seeds 20261024 and 20261025 each had one such
  case, with COAL's output equal to PAL2NAL's on the original input). `test_cpp_vs_reference` is skipped on
  Windows with the reason: its shim runs coal-cpp as a child process there, so the reference tests that interrupt
  a run watch the wrong process (before, the test errored on the symlink privilege before getting that far).
- `tests/tools/stress_*.py`: the stress-test generator, harness and helpers; `tests/test_stress_findings.py`.
- **CI: the app's real-engine tests.** The GUI and MCP app jobs build coal-cpp and install DIAMOND, and run the
  app's real-engine Ka/Ks and pipeline tests (`TestRealEngine`) on Linux and Windows.
- **CI: `tests.test_engine_ux` runs in the C++ job** (and the nightly sanitizer job); it was skipped in every
  job before (#85).
- **CI: every coal-cpp test runs on every push** (#86). The C++ engine job runs every module in `tests/` against
  coal-cpp on Linux, Windows and macOS, with DIAMOND, Clustal Omega and MAFFT installed (`.github/ci_cpp_tests.py`:
  a "not built" skip fails anywhere), so the 77 coal-cpp tests of 16 modules that ran in no workflow now run, and
  `tests/test_repo_ci.py` checks that every module with a coal-cpp skip is in that job and the sanitizer job; a
  new x86_64 macOS job (`cpp-macos-x86`, AppleClang) runs the AVX2 Ka/Ks copy; the nightly sanitizer job runs the
  same set (less `test_wheels`) under ASan / UBSan, a new nightly job `cpp-extras` runs the tests that need more
  than the binary, and the TSan build works again (GCC 15's false `-Werror=free-nonheap-object` under sanitizers
  is a warning in that job only); `tests/test_kaks_dispatch.py` tells a build without the AVX2 copy by the build,
  not by the OS.
- **CI: the tutorials run as tests.** `scripts/run_tutorials.py` runs the commands of `docs/tutorials/` (01, 02,
  03, 05) on their real NCBI data, one bash session per tutorial, and fails when a documented output line, an exit
  status or one of 15 numbers in the prose no longer matches what the code prints; documented lines that are known
  to differ are listed in `tests/fixtures/tutorials/known_drift.json` (empty now) and reported as warnings
  (`--strict-known` fails on them too). The CI job `tutorials` (not a required check) runs it when the tutorials,
  `src/coal/`, `cpp/` or the runner change; `tests/test_tutorials.py` tests the runner offline. The tutorials'
  code blocks are tagged `bash` (commands) and `text` (outputs), as the runner reads them. The app's two timeout
  tests force `COAL_ENGINE=python`, so they no longer fail when coal-cpp is found.
- **CI: release dry run and smoke test.** `release.yml` can be dispatched on a branch with `publish=false` (a dry
  run of every build job, nothing published; before, it required a release tag, so main could not be built before
  tagging). `release-smoke.yml` (dispatch only, read-only token) installs the files of a release.yml run on Linux
  x86_64 and arm64, macOS and Windows with Python 3.10 and 3.13 (the platform wheel, the pure wheel and the sdist,
  each in a fresh environment), runs them (alignments with both engines, `coal dnds`, `coal pipeline pairs
  --pairs`, the app's `--check`, the Windows app's portable zip and silent install), and requires every file
  written with `-o` to be byte-identical across systems, Pythons and install kinds.
- **Tests: property-based and differential fuzz tests of the Ka/Ks engine.** `tests/test_kaks_properties.py`
  (seeded random pairs, 1 to 2,000 codons, all nine methods, all genetic codes: no crash or hang, documented
  statuses only, no nan/inf, invariance to row swap, codon-column order, `-j` and input format, `python -m coal`
  equal to `coal-cpp`; a fast subset by default, `COAL_SLOW=1` for more) and `tests/fuzz_kaks_vs_reference.py`
  (random pairs against KaKs_Calculator 2.0 and PAML yn00, every difference classified by the documented ones;
  skipped unless the references are given). K01, K02, K03, K05 and K06 are fixed and their `TestKnownBugs`
  tests pass; `tests/test_kaks_myn_kappa_inf.py` (and `tests/cpp/kaks_myn_kappa_inf.cpp` through the C++ API)
  covers the infinite kappa. ci.yml's C++ job runs `tests.test_kaks_properties` and the MYN test; the nightly run
  has a job `kaks-properties` (COAL_SLOW unit tests and a 1,500-pair campaign) and runs 300 property pairs under
  the sanitizers.
- **CI: the tutorials every night.** nightly.yml has a job `tutorials`: ci.yml's tutorials job without its paths
  filter, so drift that comes from a new tool pin or from NCBI's files is caught too (CONTRIBUTING.md "Nightly
  checks" lists all six nightly jobs).
- **Docs:** docs/mcp.md describes the four MCP tools, their arguments, results, errors and limits; docs/install.md,
  the Claude skill (`.claude/skills/coal/SKILL.md`) and the MCP stage demo (`demo/mcp_demo.md`) no longer say the
  app or its MCP server has two tools or covers alignment only. docs/install.md has a section "Clustal Omega and
  MAFFT on Windows"; docs/r.md describes coalr's positional-pairing and all-failed warnings; docs/usage.md says
  what `coal translate` warns about; docs/pipeline.md describes the first-use messages of the new-user audit (a
  `--pairs` / `--clusters` file that names no record stops with exit 1, the header line, CSV, genome-like inputs,
  the earlier-results notice, option suggestions) and the Windows routes; the release checklist
  (docs/development.md) removes the dev-version banner (#81, #83). docs/alignment.md's Windows install advice
  names the same routes as install.md (MSYS2's Clustal Omega 1.2.4, MAFFT's own Windows package, why COAL
  downloads neither). The "Tested with" lines of tutorials 1, 2, 3 and 5 no longer name an old commit: CI reruns
  those four and compares their outputs with the text. Tutorial 4 (R: coalr and orthologr) is not run by CI; its
  line names the commit it was run with (009bdd5) and says that it is not rerun.
- **Tests added with the review of 2026-11-13:** the pair-table header is the same in `dnds_cli.cpp`,
  `kaks_cli.cpp`, `dnds_pipeline.py` and coalr (`tests/test_pair_table.py`); a `coal-cpp` built without test hooks
  makes the reference-shim resume tests skip instead of fail; five fixed sleeps in the lock, MSA and pipeline
  tests are event or polling waits (test hook `COAL_TEST_LOCK_BEFORE_TAKEOVER_UNTIL`); the YN / MYN shortcut is
  compared with the full computation on every platform (`--yn-always-pmat`, a testing option).
- **Bioconda recipe** (`packaging/bioconda`): written for 0.6.0 with the version and sha256 marked as
  placeholders for release time; both binaries; GPL-3.0-only with `LICENSE` and `NOTICE`; the new scope in the
  summary; DIAMOND, Clustal Omega and MAFFT under `run_constrained`; tests of `coal dnds`, `coal kaks`,
  `coal pipeline --help` and `coal cite`.
- **Website and README** present the full scope (codon alignment, Ka/Ks, the pipeline, coalr, orthologr) with
  one runnable example each; the website's new sections are marked `data-since="0.6.0"` and stay hidden on the
  deployed site until that release. No document tells users to `pip install coal-align` as if it were on PyPI.

## 0.5.0

- **Run report: a new top-level `engine` field** (`coal.report/1`, next to `source`): `"coal-python"` or
  `"coal-cpp"`, the same strings as the `--resume` checkpoint's `engine`, so a pipeline or the app can tell which
  engine wrote a report (also when `coal` handed a large input to `coal-cpp`). The per-job reports of a
  `coal --batch` report (`coal.batch/1`) carry it too. Apart from `engine` and `runtime_seconds` the two engines'
  reports stay identical (#18).
- **Changed: `summary.not_checked` of a run that stopped at an error now counts only the sequences the run never
  reached**, the number the `coal: summary:` line "N sequences were not checked" gives; the sequences checked
  before the error are counted as what COAL found (`ok`, `warning`, `failed`), so `ok + warning + failed +
  not_checked == total`. It used to count every sequence but the failed one (e.g. `examples/one_broken_gene`
  without `--keep-going`: `ok 0, failed 1, not_checked 19`; now `ok 13, failed 1, not_checked 6`). The per-sequence
  entries are unchanged (`not_checked` with `not_written` before the error, `not_reached` after it). Both engines;
  documented in `docs/report-schema.md` (#18).
- **`pip install coal-align` brings the C++ engine** (issue #16): each release now has one wheel per platform with
  `coal-cpp` inside, at `coal/bin/coal-cpp` in the installed package: Linux x86_64 and aarch64 (`manylinux_2_28`,
  glibc 2.28+, `libstdc++` linked statically), macOS 11+ (one `universal2` wheel for Apple silicon and Intel) and
  Windows x64 (static CRT). pip picks it by itself from PyPI; from a GitHub Release, install the wheel file for your
  platform. The pure `py3-none-any` wheel and the sdist stay as the fallback everywhere else. `coal` finds the bundled
  binary with no settings: the lookup order is now `COAL_CPP`, then `coal/bin` inside the package, then next to
  Python, then `PATH` (docs/cpp-engine.md "Where coal looks for coal-cpp"). New `coal --cpp-path` prints the
  `coal-cpp` that `coal` would use (exit 1 with the reason when there is none or it is another version), for
  pipelines that call the binary directly. The wheels are made by `packaging/wheels/build_platform_wheel.py` (no
  cibuildwheel, no compiled extension) in `release.yml`'s new `cpp-wheels` job and checked there (auditwheel, lipo, a
  fresh-venv install whose output equals the Python engine's); `tests/test_wheels.py` does the same with a local build.
- **Nightly fuzz workflow** (`.github/workflows/nightly.yml`, 03:17 UTC, or by hand with a seed): the fuzz against
  the real PAL2NAL v14 with a fresh seed each night, once per engine (`COAL_FUZZ_ENGINE=cpp` in
  `tests/fuzz_vs_pal2nal.py` runs `coal-cpp`), the new `tests/tools/engine_fuzz.py` (coal-cpp vs the Python engine
  on malformed inputs and every option; Linux and macOS) and the C++ tests plus the generator under
  AddressSanitizer and UndefinedBehaviorSanitizer (`COAL_SKIP_MEMORY_TESTS=1` skips the peak-memory test there).
  Differing cases are uploaded as artifacts; CONTRIBUTING.md "Nightly checks" says how to replay them.
- **coal-cpp: no read past the end of a short gene's DNA** when the expected position of an anchor (after a
  frameshift digit or a mismatch) lies beyond the sequence end; found by the first sanitizer run. Output and exit
  codes are unchanged (the out-of-bounds read could never have produced a match).
- **MCP server on mcp 1.x (the Windows app's) is quiet and names its version:** `python -m coal_app --mcp` /
  `COAL.exe -m coal_app --mcp` with the `mcp` 1.x that Gradio's extra installs reported the mcp library's version
  ("1.10.1") as `serverInfo.version` and logged about 37 KB of INFO lines to standard error per session (into
  Claude Desktop's MCP log). It now reports COAL's version, as on mcp 2.x, and logs warnings and errors only. The
  message for a missing `mcp` package no longer claims that version 2 is required (1.x works too).
- **Packaging:** the source distributions of `coal-align` and `coal-align-app` no longer ship `tests/` (the tests
  need the repository: helpers, fixtures, examples, coal-cpp), and the `coal-align-app` wheel and sdist now carry
  the `LICENSE` (MIT at the time; the app had none).
- **App: a settings file that is not text no longer stops COAL from starting.** `settings.json` in another encoding
  (or binary garbage) gave a traceback at start-up; it now counts as a damaged file: the defaults are used and the
  Settings tab says so, as for a half-written file.
- **App: an unusable results folder is one plain sentence, not a traceback.** With `COAL_RUNS_DIR` pointing where no
  folder can be made (a path under a missing directory, an existing file, no permission), the GUI shows "Could not
  create the results folder <folder>: <reason>" in the banner, the MCP tool `codon_align` raises a tool error with
  that sentence, and `coal-gui` / `python -m coal_app` print it and exit 1 instead of failing in `os.makedirs`.

- **Device inputs are refused** (robustness review): `coal /dev/zero nuc.fa` or `coal - nuc.fa </dev/zero` used to
  copy the device to a temporary file without end, filling the disk silently. Both engines now stop before reading
  with `coal: error: cannot read /dev/zero: it is a device, not a file` (exit 1; `standard input` for `-`). A
  terminal on standard input and `/dev/null` behave as before; pipes, FIFOs and `<(...)` are unchanged.
- **`--report PATH` no longer replaces what is at PATH** (robustness review): it used to rename a temp file over
  anything there, so `--report nuc.fa` destroyed an input, `--report out.fa` the output, and a symlink, FIFO or
  device node (`/dev/null` became a regular file!) was replaced. Both engines now give the report path the rules of
  `-o`: an input or the output file is a usage error before anything runs (`coal: error: the report path r.json is
  an input file` / `is the output file`, exit 2; also for a job line's `--report` and the batch report), symlinks
  are written through, and FIFOs, devices and hard-linked files are written in place. A report into a folder that
  does not exist is a usage error too (`the report path results/r.json: folder /data/results does not exist`, exit
  2) instead of a full run that ends by losing its report, and the report's temp file is removed on SIGTERM /
  SIGHUP / Ctrl-C while it is written, in both engines.
- **A killed worker no longer hangs the Python engine** (robustness review): when the out-of-memory killer (or
  `kill -9`) took one worker of a `-j` run, the main process could block forever in the pool's shutdown. It now
  ends the pool without waiting and stops with `coal: error: a worker process died (out of memory?)` (exit 1, report
  code `worker_died`, no partial output, no temp files). Workers also let SIGTERM end them (the inherited handler
  used to be swallowed inside the work) and **die with the main process** (`PR_SET_PDEATHSIG` on Linux, a parent
  watchdog everywhere), so a `kill -9` of `coal` no longer leaves workers holding its stdout open and a pipeline
  reading it waiting for EOF forever.
- **Python 3.14 on Linux: every multi-worker run failed** with `a worker process died (out of memory?)` after
  half a second. 3.14 made `forkserver` the default start method there, so the workers' parent is the fork server,
  not `coal`, and the watchdog above took that for a dead parent. The workers now wait for EOF on a pipe whose only
  write end `coal` holds, whatever the start method; that also fixes **Windows, where a killed `coal` left its worker
  processes running** (there `getppid()` keeps returning the dead parent's pid). `COAL_MIN_NT_PER_JOB` (default 4,000,000)
  can now be set in the environment, as `COAL_BATCH_NT`, so tests can force worker processes on a small input; CI
  runs the suite on 3.14 (Linux).
- **Worker processes end as soon as the main process fails** (Python engine): after an error on the writing side
  of a `-j` run (a full disk, a file-size limit) the pool was only ended when the cyclic garbage collector got to it
  (the suspended result generator and an exception held in a local kept it alive), so terminated workers and pool
  threads could outlive the error by seconds. The generator is now closed on the way out and the pool ended at once.
- **Checkpoints are checked before `--resume` trusts them** (robustness review): an edited or corrupted
  `<out>.coal-checkpoint.json` (`completed` beyond the input or negative, a string or float where an integer belongs,
  a list at the top level, a missing key, a deleted entries file, a partial file shorter than recorded, a truncated
  or deeply nested file, a folder at the checkpoint's path) made both engines write a truncated output with exit 0,
  or crash with a traceback / raw C++ text. Both engines now check every field's type and range and the partial
  files' sizes, treat any problem as the existing "unreadable checkpoint" notice and start a normal run; a
  checkpoint that cannot be written (`cannot write the checkpoint ...: is a directory`) is a notice and the run
  goes on without one; and `coal-cpp`'s JSON parser has a nesting limit (64) instead of overflowing the stack.
  The `drop` and `bad` column lists are bounded by the run's protein and codon columns too: a crafted
  `"bad": [4000000000000]` made `coal-cpp` fail with a raw `std::bad_alloc`; it is an unreadable checkpoint now, in
  both engines.
- **`--resume` of a `--keep-going` run in fasta or paml format** (Python engine): a checkpoint written after a
  skipped sequence was rejected as "unreadable" (the skipped sequence's index was recorded for clustal / codon
  only, and the checkpoint check requires one index per skipped sequence), so the run started over. The index is
  now recorded in every format, as `coal-cpp` does, and the resumed output is byte-identical to an uninterrupted run.
- **A checkpoint written by the other engine** (robustness review): `coal --resume` now picks the engine that wrote
  the checkpoint (so a `coal-cpp` installed or removed between the run and the resume does not get in the way),
  and when an engine is asked to resume the other one's run it says so instead of blaming the inputs: `out.fa was
  interrupted by the Python engine (coal-python), which must resume it: run again with COAL_ENGINE=python, or start
  a new run without --resume`.
- **`-o /dev/fd/1` and `-o /dev/stdout` work when stdout is a pipe** (robustness review): both engines resolved
  the path to `pipe:[N]` and failed with `cannot write /dev/fd/1: No such file or directory`. An existing target
  that is not a regular file is now written through before any path resolution.
- **`--report /dev/fd/1` and `--report /dev/stdout` work when stdout is a pipe**, in both engines, for the same
  reason (the report failed with `could not write the report /dev/fd/1: No such file or directory` after the
  whole run).
- **`--json | head` is a silent exit 141** (robustness review) in both engines (the Python engine printed a
  traceback and exited 1, `coal-cpp` exited 0), also for `--batch ... --json`. `coal-cpp` now streams the `--json`
  report to standard output one sequence at a time instead of building it whole in memory first (the bytes are
  the same).
- **`coal-cpp` falls back when `TMPDIR` is unusable** (robustness review): a `TMPDIR` naming a missing folder or a
  file made it fail with `filesystem error: temp_directory_path`; it now searches `TMPDIR`, `TEMP`, `TMP`, `/tmp`,
  `/var/tmp`, `/usr/tmp` and the current folder as Python's `tempfile` does, so both engines behave the same.
- **Disk-full errors name the file** (robustness review): the Python engine printed `[Errno 28] No space left on
  device` (and could hang at exit after terminating its workers), `coal-cpp` `could not write the output (disk
  full?)` (and, with the rows going through a full temp folder, exit 0 with an empty output). Both now say
  `cannot write out.fa: no space left on device`, `cannot write a temporary file in /tmp: ...` or `cannot copy
  standard input to a temporary file in /tmp: ...`, exit 1, and leave no temp file behind.
- **Python API checks its arguments** (robustness review): `coal.codon_align` raises `CoalError` with
  `code="duplicate_id"` for a repeated protein or DNA ID (it used to pair them silently; the command line refuses
  the file), and both `codon_align` and `coal.align_files` raise `TypeError` with a message naming the argument for
  a dict, generator, bare string, `None`, a non-`str` ID or a pair of the wrong length, instead of a `KeyError` or
  `AttributeError` from inside the engine. `align_files(pep, "nuc.fa")` with one DNA path as a string or
  `os.PathLike` now means that one file (the string used to be iterated character by character). `nuc_paths` is
  one path or a `list` / `tuple` of paths: a generator (or any other iterator) of paths, which used to work, is a
  `TypeError` now; pass `list(...)`.
- **Smaller robustness fixes** (both engines, identical messages): a NUL or other control byte (except tab, LF,
  VT, FF, CR) anywhere in a protein or DNA file is `<file>:<line>: control character: sequence files must be plain
  text (is this a binary or word-processor file?)` (exit 1, code `input_format`) instead of a confusing row-length
  error, a mangled ID or silently dropped bytes; IDs longer than 200 characters are shortened on stderr to their
  first 100 and last 20 characters (the output and the report keep them whole); `-o -` means standard output; the
  Python engine reports a DNA file that shrank between indexing and reading as `could not read the DNA of record N
  (the file changed after it was read, or cannot be read)`, as `coal-cpp` does, instead of reading another
  record's bytes; an empty DNA file among several is documented as accepted.
- **Mismatching long genes no longer take quadratic time.** When a protein does not fit its DNA exactly, both engines
  locate PAL2NAL's 10-residue chunks in the DNA; every chunk was a scan from the start of the DNA (coal-cpp; the
  Python engine scanned for every chunk that was not where the previous ones put it), so a 35,000-codon gene with one
  wrong codon cost about 10 s per gene with coal-cpp (and the genetic-code hint repeated that under 17 other codes).
  The DNA is now translated once per frame and its 10-residue windows indexed, so a chunk of ordinary residues is a
  lookup; chunks with a frameshift digit, an X or the first Met keep the old scan. Same chunks found, same output and
  messages: 10 genes × 35,000 codons with one wrong codon each, 104 s → 3.0 s (coal-cpp), 27 s → 16 s (Python);
  exactly fitting genes are unchanged. The Python engine's fallback is also about 1.5× faster per short gene.
- **coal-cpp uses far less memory on runs with many warnings or with `--report`.** It kept a JSON tree for every
  sequence with a warning (about 4 KB each: 1 GB at 300,000 mismatching sequences) and, for `--report`/`--json`, built
  every sequence's entry before writing (3 KB each: 3.2 GB at 1,000,000 sequences). It now keeps each non-ok entry as
  its compact JSON line (about 0.6 KB; the Python engine's dict entries take 1.4 KB) and streams the report one
  sequence at a time:
  300,000 mismatching sequences 994 → 266 MB; 1,000,000 clean sequences with `--report` 3.2 GB → 204 MB. The
  Python engine streams its report the same way and no longer builds a view of every sequence for the genetic-code
  hint (1,000,000 sequences with `--report`, one process: 635 → 367 MB). The report JSON is unchanged.
- **Docs:** `docs/usage.md` gives the measured memory per sequence of both engines (about 350 bytes with the Python
  engine after the change above, 190 with coal-cpp; the old figures, 300 bytes and 122/139/145 MB, matched coal-cpp
  only), where each engine puts its temporary files and how much disk the clustal/codon formats need, that
  `-q --progress` still prints progress, and that the checkpoint cadence is best-effort under I/O contention.

- **Less memory with many workers:** COAL's own processes now limit numpy's BLAS thread pool to one thread
  (`OPENBLAS_NUM_THREADS`, `OMP_NUM_THREADS`, `MKL_NUM_THREADS`, unless already set). COAL uses no BLAS, but every
  worker that imported numpy started a pool sized to the machine. On a 20-thread Windows machine each process that
  imported numpy committed about 630 MB (9-11 GB at the default `-j`); the BLAS pool is the likely cause, to be
  confirmed by re-measuring. Programs that use COAL as a library are not affected.
- **Faster start-up** for runs on small files (one COAL call per gene family in a pipeline): about 40 ms less per call
  (104 -> 64 ms on Linux), by not importing `dataclasses`, `subprocess`, `hashlib`, `json` and `shlex` unless needed.

- **`-` reads standard input** (`zcat nuc.fa.gz | coal pep.aln -`), for the protein alignment or one DNA file, in
  `coal`, `coal-cpp` and the `pal2nal` drop-in, also on Windows (no `/dev/stdin` needed). The data is read as bytes
  into a temporary file, as for other pipes, so the output is byte-identical to giving the file's path; messages name
  it `standard input`; a second `-` is a usage error. The run report records `"path": "-"` with the sha256 of the data
  read (the report of a run reading `/dev/stdin` or another pipe now hashes the data too, in both engines, instead of
  the Python engine's hash of an empty re-read).
- **C++ library for embedding** (e.g. in DIAMOND): `cmake --install` puts `coal-cpp`, `libcoal_core.a`, the headers
  (`<coal/core.hpp>`, `<coal/tables.hpp>`, `<coal/seqio.hpp>`, `<coal/version.hpp>`) and a CMake package in place;
  `find_package(coal CONFIG)` or `add_subdirectory` / FetchContent give `coal::core`. As a subproject only the library
  is built (`COAL_BUILD_CLI`, `COAL_BUILD_TOOLS`, `COAL_INSTALL` to add more). Example: `examples/cpp-embed/`; docs:
  docs/cpp-engine.md "Embedding the library".
- **orthologr stand-in takes every PAL2NAL option:** `packaging/orthologr/pal2nal.pl` sent every command line to
  `coal-cpp` when installed, which rejects PAL2NAL-only options such as `-nostderr` and prints `coal`'s help for
  `-h`. It now runs `coal-cpp` only for the options it understands (`-output`, `-codontable`, `-nogap`,
  `-nomismatch`, `-blockonly` and the files) and COAL's `pal2nal` drop-in for anything else; `COAL_PAL2NAL` still
  overrides.
- **orthologr:** `packaging/orthologr/pal2nal.pl` is a stand-in for the `pal2nal.pl` bundled in the orthologr R
  package; copied over it, orthologr's `codon_aln()` / `dNdS()` run COAL (coal-cpp when installed) with their
  unchanged command line and get the same output. Notes on what to expect and on a proper `tool = "coal"` option
  for orthologr: `packaging/orthologr/README.md`.

- **Python API: every error is a `CoalError`** (issue #18). `coal.codon_align` and `coal.align_files` with an
  unsupported `table` now raise `CoalError` with `code="unknown_table"` (the message lists the supported NCBI codes)
  instead of a bare `ValueError`.
- **Python API: `coal.seqio.InputFormatError` is a `CoalError`** with `code="input_format"` (issue #18), so one
  `except CoalError` catches a file that is not FASTA or CLUSTAL too; it stays a `ValueError` for code that catches
  that. The `coal` command's messages, exit code and report are unchanged.
- **`pal2nal | head` no longer prints a traceback** (issue #18): the drop-in's usage text into a closed pipe is a
  silent exit 141, as for `coal` and for the alignment output.
- **Batch mode (`coal --batch JOBS`, #15):** thousands of small conversions in one process, one `coal` command line
  per line of JOBS (`-` = stdin), the command line's options as the defaults of every job, `--batch-out-dir DIR` for
  jobs without `-o`, a `coal.batch/1` report with `--report`/`--json`, per-job stderr prefixes
  (`coal: [k/n PEP] warning: ...`) and a `coal: batch:` summary line. A failing job does not stop the batch (exit 1
  at the end; 2 if jobs only left sequences out; mistakes in JOBS are reported with their line number before
  anything runs). Every output is byte-identical to the single run; about 40× less wall time on tiny jobs. Python
  only: `coal-cpp --batch` is a clear error and batches are never dispatched to it. `--resume` and cross-job
  parallelism are not supported in a batch yet.
- **The app no longer waits forever for a stopped engine:** when an engine run hits the app's time limit (or the
  demo pipeline's), it is sent SIGTERM and, after the grace period, SIGKILL; the kill now takes the engine's whole
  process group, so its worker processes die with it. Before, a run that outlived the grace period on a slow machine
  was killed alone, its orphaned workers kept its output pipes open and the app waited for them indefinitely (the
  CI app suite hung until the job timeout). Windows already killed the tree. Test commands that start the engine
  now give it no stdin and always carry a timeout.
- **`coal-cpp` mirrors Python 3.13's argparse** (it mirrored 3.10's): `-qx` is "unrecognized arguments: -x" rather
  than an ignored explicit argument, `-hx` prints help, `-cod=2` is accepted, invalid choices are quoted
  (`invalid choice: '7' (choose from '1', '2', ...)`), and only the `--` marker is dropped from the positionals (a
  second `--` is a file name). `tests/test_cpp_cli.py` compares the two command lines live on an interpreter
  whose argparse behaves that way (probed, not keyed on the version) and skips elsewhere with the differences named.
- **Documented for 1.0:** `docs/python-api.md` (the public Python API: `coal.codon_align`, the new re-export
  `coal.align_files` = the `coal` command as a function, `coal.CoalError` with its `code`/`seq_id`/`index`, and the
  semantic-versioning promise for command-line flags, output bytes, exit codes, the report schema and these
  functions), `docs/report-schema.md` with a JSON Schema (`docs/report.schema.json`, draft 2020-12) for the run report
  `coal.report/1` (every field, status, pairing mode, message and hint code; additive changes only within `/1`),
  `CONTRIBUTING.md`, GitHub issue forms and a pull-request template. New tests run the documented examples
  (`tests/test_python_api.py`) and validate freshly generated reports from both engines against the schema
  (`tests/test_report_schema.py`); `tests/report_example.json` is now a real `examples/one_broken_gene --keep-going`
  report. No change to COAL's behaviour or output.
- **orthologr:** `packaging/orthologr/pal2nal.pl` is a stand-in for the `pal2nal.pl` bundled in the orthologr R
  package; copied over it, orthologr's `codon_aln()` / `dNdS()` run COAL (coal-cpp when installed) with their
  unchanged command line and get the same output. Notes on what to expect and on a proper `tool = "coal"` option
  for orthologr: `packaging/orthologr/README.md`.
- **Repository tidied for outside users:** a short README with a section on using COAL in pipelines; the detailed
  documentation is in `docs/`. Hackathon logs and pitch material are removed (kept on the `hackathon-archive` tag),
  review-round markers are gone from comments, and engine test modules are named after what they test. No change to
  COAL's behaviour or output.
- **`coal-cpp` reports the same error as `coal` when a protein file has two defects.** A FASTA record is now
  handed back before the `>` line that ends it is checked, so a short row is reported before a later duplicate or
  empty header; and a CLUSTAL file's rows are parsed before its `#` mark line is read, so a malformed CLUSTAL line is
  reported before a later non-ASCII character. Found by differential fuzzing; `tests/test_cpp_error_precedence.py`.
- **Package renamed to `coal-align`** (`coal` on PyPI is another project). The `coal` command, the `coal` Python
  module and `pal2nal` are unchanged; release files are now `coal_align-<version>-py3-none-any.whl` and
  `coal_align-<version>.tar.gz` (`coal_src.tar.gz` stays).
- **The app is its own package, `coal-align-app`:** the web GUI, the MCP server, the viewer, the Windows app and the
  Hugging Face Space moved to `apps/coal-app/` (package `coal_app` unchanged; `coal-gui` and `python -m coal_app` work
  as before). `coal-align` is now the engine alone (`coal`, `pal2nal`, the examples) with no dependencies and no app
  code, so pipelines (OrthologR, DIAMOND workflows) install only that. The `gui` and `mcp` extras are on the app:
  `pip install "coal-align @ .../coal_src.tar.gz" "coal-align-app[gui] @ .../coal_app_src.tar.gz"`, or from a checkout
  `pip install -e . -e "apps/coal-app[gui]"`. Both packages always have the same version, and the app depends on that
  exact `coal-align`. Releases add `coal_align_app-<version>-py3-none-any.whl`, `coal_align_app-<version>.tar.gz`
  and the version-free `coal_app_src.tar.gz`. The Windows launchers are now `apps\coal-app\desktop\run_from_source.bat`
  and `apps\coal-app\desktop\setup_mcp.bat`. The app's tests (`apps/coal-app/tests/`) are grouped by feature; the
  checks of the docs, the website and CI moved to `tests/` (`test_docs_*.py`, `test_repo_ci.py`).

- **Bioconda recipe prepared** (`packaging/bioconda/`): `meta.yaml`, `build.sh` (coal-cpp via CMake into `bin/`, then the
  Python package), a fixture pair and a test that `coal` finds the installed `coal-cpp` and hands a 1 MB input to it.
  Not submitted yet; `packaging/bioconda/README.md` explains the sha256 step and the pull request to `bioconda-recipes`.

## 0.4.2

- **COAL in Claude Desktop, from the installer:** the Windows installer's last page has **Connect COAL to Claude
  Desktop** (ticked). Claude Desktop then starts COAL itself (`COAL.exe -m coal_app --mcp`) and can run it on your
  files; uninstalling removes the entry. `codon_align` can save the alignment where you ask (`save_as`), e.g. next to
  the input files.
- **C++ engine included:** the Windows app ships `coal-cpp.exe`, and `coal` hands inputs of 1 MB or more to it when it
  is found and is the same version (same output, faster on large inputs); `COAL_ENGINE=python|cpp` chooses.
- **Quieter, clearer engine output:** a wrong genetic code is reported once, before at most 20 per-codon warnings
  (`--max-warnings`, 0 = all; the `pal2nal` drop-in still prints them all); `--help` is grouped and short.
- **MCP results tell the truth:** `codon_align` returns `completed` (COAL finished), `all_ok` / `needs_attention`
  (whether every sequence matched its protein), `counts`, `suggested_codontable` and `next_step`, and a short summary
  of COAL's messages instead of raw stderr. The old `ok` field (which meant "finished") is gone. A missing file, an
  unknown genetic code or an unsupported format is now a tool error (`isError`), and the arguments have descriptions
  and allowed values. The server reports its version.
- **Windows, from source:** `desktop\run_from_source.bat` (the window) and `desktop\setup_mcp.bat` (the MCP server
  for Claude, printing ready-to-paste registrations) use separate environments, pick Python 3.12 when available,
  check the Python version and repair half-finished installs.
- **Windows app:** the MCP endpoint is off by default (it was on a random port no client could be registered with);
  `COAL.exe --mcp` turns it on.

## 0.4.1

The first release with a working Windows installer (the 0.4.0 build failed, so 0.4.0 has no downloads).

- **Clearer app:** a short header, section headings (Files, Genetic code, Output, Results), the genetic code first,
  a "Start alignment" button, and no empty results area before a run. When the genes fit another genetic code, the
  result says "Wrong genetic code?" with a one-click re-run instead of a page of warnings. Long lists of genes with the
  same problem are folded. One settings place (Gradio's own settings link is hidden). Amber accent colour.
- **Large inputs:** the alignment view is an honest preview at scale (problems first, overview strips, centre on any
  column), and a Settings tab remembers your defaults.
- **Better messages:** errors say what to do next (e.g. "Add --keep-going to write the other sequences"), swapped or
  unaligned input files are recognised, missing files are named plainly, and frameshift locations are only given as
  "at column N" when the evidence confirms it.
- **Windows app fixes:** the packaged app no longer crashes on start-up (UTF-8 output), and file paths with any
  characters work.
- **Examples** set their own genetic code, so one example can't leave the wrong code behind for the next.
- **Tested on real data:** 17 published codon alignments are rebuilt exactly; 12 ship as tests.

## 0.4.0

COAL turns a protein alignment plus the unaligned coding DNA into a codon alignment, like PAL2NAL v14, with the same
output apart from a short list of documented, deliberate differences (see the README).

- **Windows app:** an installer (or portable zip) with a local GUI: a plain-language health check for every gene, an
  alignment view that marks where a gene breaks, one-click re-run with a suggested genetic code, and a ready-to-paste
  methods sentence. Your data stays on your computer.
- **All four PAL2NAL output formats** (`-output fasta|clustal|paml|codon`), `-blockonly`, and a drop-in `pal2nal`
  command that accepts an unchanged PAL2NAL command line.
- **Run report** (`--report FILE.json` / `--json`) with each gene's status, the genetic-code suggestion and the methods
  sentence; `--keep-going` writes every other gene when one can't be aligned (exit 2); `--strict` refuses
  probably-wrong data; `--max-warnings`, `--progress` and `--resume` for very large runs.
- **Safer pairing:** IDs that match only after removing a version suffix, case or `sp|…|` fields are paired by ID,
  and COAL says so.
- **AI agents:** an MCP server with two tools, `codon_align` and `explain_report`.
- **C++ engine** (`coal-cpp`, built from source with CMake): same options and byte-identical output, 3-6x faster
  than the Python engine in our measurements.
- **Tested:** 285 stored checks against PAL2NAL v14's own output and 20,000 random comparisons with the real PAL2NAL
  v14, with 0 unexplained differences.

## 0.3.0

- **Large inputs:** a streaming pipeline with constant memory, automatic multi-process `-j` (only when each worker
  gets enough work) and an optional numpy fast path (`".[fast]"`); 500,000 sequences: 2.5 GB -> 143 MB of memory.
- All four PAL2NAL output formats and `-blockonly`, `--strict`, the drop-in `pal2nal` command, the GUI, MCP tools
  and the Windows app (see README.md).

## 0.2.0

- NCBI genetic codes 1-23 (`-codontable`), extra sequence around the gene (UTRs, stop codon) skipped, PAL2NAL's
  frameshift digits, `-nogap` and `-nomismatch`, and PAL2NAL's handling of codons that don't match.
- About twice as fast; on Windows 29x faster than PAL2NAL on 10,000 sequences and 119x on 10,000-codon genes,
  with identical output.

## 0.1.0

- First version: a pure-Python, dependency-free `coal` command that builds a codon alignment from a protein
  alignment and coding DNA, checked against stored PAL2NAL v14 outputs.
