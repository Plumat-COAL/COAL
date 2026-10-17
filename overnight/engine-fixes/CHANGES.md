# Engine fixes (branch keanu_overnight/engine-fixes)

## Summary

(written at the end of the run)

## Changes

Source of the findings: independent engine review, items a–i. No earlier run had fixed any of them
(branch started clean from main, `overnight/` did not exist). Tests: `tests/test_review_fixes.py`
(class names `TestEng<ID>_…`). Run with `PYTHONPATH=src python -m unittest discover -s tests`.

| ID | Where (file:function) | Before | After | Why (finding / decision, with evidence) | Test proving it | Merge note |
|---|---|---|---|---|---|---|
| ENG-a | src/coal/pipeline.py:run | `-o PATH` always wrote a temp file next to PATH and `os.replace`d it: a symlink was replaced by a regular file (target never written); FIFOs/devices were replaced or failed with "Permission denied: '/dev/.coal_….tmp'"; new files had mode 0600. | PATH is resolved with `os.path.realpath`. If it is a regular file or does not exist: temp file in the real target's directory, `chmod 0o666 & ~umask`, then rename onto the real target (atomic, no partial output, symlink kept). Otherwise (FIFO, device, …): output is staged in the system temp dir and then copied into PATH opened for writing. | Review (a); reproduced: symlink replaced, mode 0600, FIFO replaced. | `TestEngA_OutputPath` (symlink, umask 022 → 0644, FIFO in temp dir) | Teammate's new output formats and the C++ port must use the same write rule. If the other branch also changed output writing, keep the realpath/regular-file check and the chmod. |
| ENG-b | src/coal/pipeline.py:run (new thin wrapper `run` + `_spool`; old body renamed `_run`, unchanged) | Inputs were opened by path several times (protein 3×, DNA re-opened and seeked): `cat p.fa \| coal /dev/stdin n.fa` and `coal <(cat p.fa) n.fa` gave "number of sequences differs (protein: 0 …)"; piped DNA gave "File or stream is not seekable". | Any input path that is not a regular file is first copied to a temp file `.coal_in_*` in the system temp dir; all temp copies are deleted in a `finally`. Regular files are used in place (no extra I/O). Error messages for a spooled input name the temp copy (suffix = original basename). | Review (b); reproduced with both commands. | `TestEngB_StreamInputs` (protein via pipe, DNA via pipe + no leftover copy, FIFOs in temp dir) | Signature of `run` unchanged (`*args, **kw` forwarded). If the teammate edited `run`'s body, apply their edits to `_run`. C++ port: needs the same spooling (or a single-pass reader). |
