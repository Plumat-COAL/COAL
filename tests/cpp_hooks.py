"""Whether the coal-cpp under test has the test hooks (docs/development.md "Test hooks").

The COAL_TEST_* hooks of the C++ engine are compiled only into builds configured with -DCOAL_TEST_HOOKS=ON; a
release build, a wheel's binary or a package has none. A test that needs them calls need_test_hooks(self, binary):
without the hooks it is skipped with the reason, or fails when COAL_REQUIRE_CPP is set (CI: the binary under test
must be a test build), unless COAL_ALLOW_NO_TEST_HOOKS=1 says that the binary is a release build on purpose (the
tests of a built wheel).
"""
from __future__ import annotations

import os
import subprocess

_cache: dict = {}
REASON = ("this coal-cpp has no test hooks: configure the build with -DCOAL_TEST_HOOKS=ON "
          "(docs/development.md \"Test hooks\")")


def has_test_hooks(binary) -> bool:
    """The binary answers the probe COAL_TEST_HOOKS=probe (cpp/src/main.cpp); one without hooks ignores it."""
    key = str(binary)
    if key not in _cache:
        try:
            r = subprocess.run([key, "--version"], capture_output=True, text=True, timeout=60,
                               env={**os.environ, "COAL_TEST_HOOKS": "probe"})
            _cache[key] = r.returncode == 0 and r.stdout.strip() == "coal-cpp: test hooks compiled in"
        except (OSError, subprocess.SubprocessError):
            _cache[key] = False
    return _cache[key]


def need_test_hooks(test, binary) -> None:
    if binary is None or has_test_hooks(binary):
        return
    if os.environ.get("COAL_REQUIRE_CPP") and os.environ.get("COAL_ALLOW_NO_TEST_HOOKS") != "1":
        test.fail(f"{binary}: {REASON} (COAL_REQUIRE_CPP is set)")
    test.skipTest(f"{binary}: {REASON}")
