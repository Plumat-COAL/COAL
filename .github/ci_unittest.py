"""Run unittest and report every failure/error as a GitHub Actions annotation.

Raw CI logs can't always be downloaded, but annotations can be read through the API, so each failing test
is printed as `::error title=<test>::<traceback tail>`. Exit status is the same as unittest's.
Usage: python .github/ci_unittest.py <unittest args...>   (e.g. `discover -s tests` or `test_wrapper test_smoke`)
"""
import os
import sys
import unittest

# Like `python -m unittest`: make modules in the current folder importable (sys.path[0] is .github/ otherwise).
sys.path.insert(0, os.getcwd())


def _esc(text):
    return text.replace("%", "%25").replace("\r", "").replace("\n", "%0A")


def main(argv):
    prog = unittest.main(module=None, argv=["ci_unittest"] + argv, exit=False, verbosity=2)
    res = prog.result
    for kind, items in (("FAIL", res.failures), ("ERROR", res.errors),
                        ("UNEXPECTED SUCCESS", [(t, "") for t in getattr(res, "unexpectedSuccesses", [])])):
        for test, tb in items:
            tail = "\n".join(tb.strip().splitlines()[-12:])
            print(f"::error title={kind}: {_esc(test.id())}::{_esc(tail)}", flush=True)
    return 0 if res.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
