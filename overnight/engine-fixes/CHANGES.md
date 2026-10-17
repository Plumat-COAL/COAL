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
