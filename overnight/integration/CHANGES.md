# Integration changes (keanu_overnight/integration = app + windows-app + engine-fixes)

Changes needed only because the branches were combined. Same marking rules as the other CHANGES.md files.

| ID | Where | Before | After | Why | Test proving it | Merge note |
|---|---|---|---|---|---|---|
| INT-1 | src/coal_app/gcode.py `_RAW[3]` | App's copy of genetic code 3 had starts ATA/ATG only | ATA/ATG/GTG, identical to the engine | Engine fix ENG-g added GTG (current NCBI table 3); the app keeps its own copy of the tables because it must not import engine internals | tests_app/test_review.py `TestFinding19NoDrift` (failed on integration, passes now) | Any future change to src/coal/tables.py must be copied to src/coal_app/gcode.py; the drift test enforces it |
| INT-2 | desktop/test_smoke.py `test_utf8_bom_is_accepted` | Marked `@unittest.expectedFailure` (engine choked on a UTF-8 BOM) | Normal test that must pass | Engine fix ENG-h now skips a BOM | The test itself (was an unexpected success) | If merged with an engine that lacks ENG-h, this test fails: that is correct, take ENG-h |
