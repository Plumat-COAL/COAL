# Integration changes (keanu_overnight/integration = app + windows-app + engine-fixes)

Changes needed only because the branches were combined. Same marking rules as the other CHANGES.md files.

| ID | Where | Before | After | Why | Test proving it | Merge note |
|---|---|---|---|---|---|---|
| INT-1 | src/coal_app/gcode.py `_RAW[3]` | App's copy of genetic code 3 had starts ATA/ATG only | ATA/ATG/GTG, identical to the engine | Engine fix ENG-g added GTG (current NCBI table 3); the app keeps its own copy of the tables because it must not import engine internals | tests_app/test_review.py `TestFinding19NoDrift` (failed on integration, passes now) | Any future change to src/coal/tables.py must be copied to src/coal_app/gcode.py; the drift test enforces it |
| INT-2 | desktop/test_smoke.py `test_utf8_bom_is_accepted` | Marked `@unittest.expectedFailure` (engine choked on a UTF-8 BOM) | Normal test that must pass | Engine fix ENG-h now skips a BOM | The test itself (was an unexpected success) | If merged with an engine that lacks ENG-h, this test fails: that is correct, take ENG-h |

## Round 2 rows

| ID | Where | Before | After | Why | Test proving it | Merge note |
|---|---|---|---|---|---|---|
| R2-D1a | tests/test_review_fixes.py (ENG-a/b/f tests) | Ran on Windows and failed there: umask test expects 0o644, `/dev/stdin` tests, SIGTERM->143, SIGINT via send_signal, symlink creation needs admin | Skipped on Windows (umask, /dev/stdin, SIGTERM, SIGINT); symlink test skips if a symlink cannot be created | windows-latest CI red on every branch (runs 38996271419 etc.); these behaviours do not exist on Windows | The skips themselves; Linux still runs all of them | Engine tests (teammate's engine/C++ must mirror only if it has the same tests) |
| R2-D1b | src/coal/pipeline.py `_run` (stdout copy) | Only `BrokenPipeError` on stdout was a quiet exit 141 | On Windows, `OSError(EINVAL)` while writing to stdout (how Windows often reports a closed pipe) is turned into the same quiet exit 141 | `coal ... | head -1` on Windows would print a traceback / "Invalid argument" | tests/test_review_fixes.py `TestR2D1_WindowsClosedStdout` (simulates the Windows error with a fake stdout; fails without the fix) | ENGINE CODE: teammate's engine/C++ must mirror. UNVERIFIED on real Windows: test_broken_pipe_is_quiet runs there unchanged and expects 141 + empty stderr |
| R2-D1c | tests_app/test_packaging.py `test_examples_found_after_wheel_install` | venv interpreter path `Scripts\python` (no .exe) on Windows -> "Could not locate Python interpreter" | `Scripts\python.exe` on Windows | windows-latest App tests failed (regression from APP-C12) | The test itself on windows-latest (UNVERIFIED here: Linux only) | none |
