## What this changes

<!-- One paragraph: the behaviour before and after, and why. Link the issue if there is one. -->

## Checklist

- [ ] Both engines: the change is made in `src/coal` **and** `cpp/src` (or it does not touch behaviour), and
      `COAL_CPP_BIN=$PWD/build/coal-cpp python -m unittest discover -s tests` passes with the binary built.
- [ ] App tests pass if `apps/coal-app` is touched: `python -m unittest discover -s apps/coal-app/tests`.
- [ ] Command-line change: `python tests/test_cpp_cli.py --write-help` was run and `cpp/src/help_text.hpp` committed;
      `docs/usage.md` and `README.md` updated.
- [ ] Output bytes, exit codes, the run report or the Python API change: called out below and in `CHANGELOG.md`
      (report fields are additive only within `coal.report/1`; see `docs/python-api.md#versioning-promise`).
- [ ] New deliberate difference from PAL2NAL v14: entry in `tests/test_differential.py` `DELIBERATE` (and the fuzz's
      `DELIBERATE` if it can generate it), bullet in `docs/differences.md`, row and counts in `tests/README.md`
      (`tests/test_docs_numbers.py` checks them).
- [ ] Tests added or updated for the change.

## Compatibility

<!-- Delete what does not apply. -->
- Output bytes: unchanged / changed (which inputs, which options)
- Exit codes: unchanged / changed
- Run report (`coal.report/1`): unchanged / fields added (which) / breaking (needs `/2`)
- Python API: unchanged / added / changed
