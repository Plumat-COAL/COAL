"""Run the tests of the C++ engine jobs (ci.yml `cpp`, nightly.yml `sanitizers`) and fail on the skips that mean
coal-cpp was not exercised.

    python .github/ci_cpp_tests.py [--no-skip MOD,MOD...] [--allow REGEX] [--exclude MOD,MOD...] <unittest args>

e.g. `... discover -s tests` (every module in tests/) or `... tests.test_kaks tests.test_align`. On top of
unittest's own result (failures and errors fail, each also printed as a GitHub Actions annotation, as
.github/ci_unittest.py does):

* a skip whose reason says a binary is missing ("not built": `coal-cpp not built`, `build/coal-cpp (COAL_CPP_BIN)
  is not built`, `coal-cpp-align not built`, ...; or "coal-cpp-align not available") fails, in every module: the
  job built the binaries, so such a skip means a test lost track of them;
* in the --no-skip modules (the core C++ comparisons) every skip fails unless its reason matches --allow (the
  platform skips);
* each --no-skip module must have run at least one test (a renamed module cannot drop out unnoticed), and a run in
  which no test ran at all fails;
* --exclude drops whole modules from the run (named in the workflow, with the reason there).

The skips are read from the unittest result (test id and reason), not from the log, so a test that prints while it
runs cannot hide one. The skip list and the totals go to $GITHUB_STEP_SUMMARY when it is set.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import unittest

# Like `python -m unittest`: make modules in the current folder importable (sys.path[0] is .github/ otherwise).
sys.path.insert(0, os.getcwd())

NOT_BUILT = re.compile(r"not built|coal-cpp-align not available", re.I)


def _esc(text: str) -> str:
    return text.replace("%", "%25").replace("\r", "").replace("\n", "%0A")


def _names(value: str) -> list[str]:
    return [m.strip().removeprefix("tests.") for m in value.split(",") if m.strip()]


def module_of(test_id: str) -> str:
    """The test module of a test id: `test_x.Class.test_y`, `tests.test_x.Class.test_y`, the class / module
    fixtures' `setUpClass (test_x.Class)`, and discover's `unittest.loader.ModuleSkipped.test_x`."""
    m = re.search(r"\(([\w.]+)\)\s*$", test_id)
    dotted = m.group(1) if m else test_id
    if dotted.startswith("unittest.loader."):
        return dotted.rsplit(".", 1)[-1]
    parts = [p for p in dotted.split(".") if p != "tests"]
    return parts[0] if parts else dotted


def _tests(suite):
    for t in suite:
        if isinstance(t, unittest.TestSuite):
            yield from _tests(t)
        else:
            yield t


def _drop(suite: unittest.TestSuite, exclude: set[str]) -> unittest.TestSuite:
    return unittest.TestSuite(t for t in _tests(suite) if module_of(t.id()) not in exclude)


def check(result: unittest.TestResult, ran: dict[str, int], no_skip: list[str], allow: str | None) -> list[str]:
    """The rule violations of a finished run (empty: none)."""
    problems = []
    allowed = re.compile(allow, re.I) if allow else None
    for test, reason in result.skipped:
        tid, mod = test.id(), module_of(test.id())
        if NOT_BUILT.search(reason):
            problems.append(f"skipped for want of a binary: {tid}: {reason}")
        elif mod in no_skip and not (allowed and allowed.search(reason)):
            problems.append(f"skipped in a module that must not skip ({mod}): {tid}: {reason}")
    for mod in no_skip:
        if not ran.get(mod):
            problems.append(f"{mod}: no test of this module ran (is it still in tests/?)")
    if not result.testsRun:
        problems.append("no test ran at all")
    return problems


def main(argv: list[str]) -> int:
    for stream in (sys.stdout, sys.stderr):   # a test id or reason the console's code page cannot show (Windows)
        try:
            stream.reconfigure(errors="backslashreplace")
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(prog="ci_cpp_tests", allow_abbrev=False)
    ap.add_argument("--no-skip", default="", help="comma-separated modules in which every skip fails (but --allow)")
    ap.add_argument("--allow", default=None, help="regex: skip reasons allowed in the --no-skip modules")
    ap.add_argument("--exclude", default="", help="comma-separated modules not to run")
    opts, rest = ap.parse_known_args(argv)
    no_skip, exclude = _names(opts.no_skip), set(_names(opts.exclude))

    prog = unittest.main(module=None, argv=["ci_cpp_tests"] + rest, exit=False, verbosity=2,
                         testRunner=_Runner(exclude))
    res = prog.result
    ran: dict[str, int] = {}   # tests collected per module (after --exclude)
    for t in _Runner.ran:
        ran[module_of(t.id())] = ran.get(module_of(t.id()), 0) + 1

    for kind, items in (("FAIL", res.failures), ("ERROR", res.errors),
                        ("UNEXPECTED SUCCESS", [(t, "") for t in getattr(res, "unexpectedSuccesses", [])])):
        for test, tb in items:
            tail = "\n".join(tb.strip().splitlines()[-12:])
            print(f"::error title={kind}: {_esc(test.id())}::{_esc(tail)}", flush=True)
    problems = check(res, ran, no_skip, opts.allow)
    for p in problems:
        print(f"::error title=C++ tests not fully exercised::{_esc(p)}", flush=True)

    lines = [f"C++ tests: {res.testsRun} run in {len(ran)} modules, {len(res.skipped)} skipped, "
             f"{len(res.failures)} failures, {len(res.errors)} errors"
             + (f" (not run: {', '.join(sorted(exclude))})" if exclude else "")]
    lines += [f"- skipped: {t.id()}: {r}" for t, r in res.skipped]
    print("\n".join(lines), file=sys.stderr, flush=True)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as fh:
            fh.write("\n".join(lines) + "\n")
    return 0 if res.wasSuccessful() and not problems else 1


class _Runner(unittest.TextTestRunner):
    """unittest's text runner, without the --exclude modules, remembering which tests it ran."""
    ran: list = []

    def __init__(self, exclude: set[str], **kw):
        super().__init__(verbosity=2, **kw)
        self._exclude = exclude

    def run(self, test):
        suite = _drop(test, self._exclude) if self._exclude else test
        _Runner.ran = list(_tests(suite))
        return super().run(suite)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
