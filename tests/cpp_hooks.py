"""Whether the coal-cpp under test has the test hooks (docs/development.md "Test hooks").

The COAL_TEST_* hooks of the C++ engine are compiled only into builds configured with -DCOAL_TEST_HOOKS=ON; a
release build, a wheel's binary or a package has none. A test that needs them calls need_test_hooks(self, binary), or
check_command(cmd, env) before it starts a command with a hook variable in its environment (the command may be
coal-cpp itself, or `python -m coal`, which hands the run to coal-cpp unless COAL_ENGINE=python): without the hooks
it is skipped with the reason, or fails when COAL_REQUIRE_CPP is set (CI: the binary under test
must be a test build), unless COAL_ALLOW_NO_TEST_HOOKS=1 says that the binary is a release build on purpose (the
tests of a built wheel).
"""
from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path

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


def engine_binary(cmd, env) -> str | None:
    """The coal-cpp that the command `cmd` runs in the environment `env`: the command itself, or the binary the
    Python command hands the run to (coal.dispatch: unless COAL_ENGINE=python); None for the Python engine."""
    if Path(str(cmd[0])).name.lower().startswith("coal-cpp"):
        return str(cmd[0])
    if (env.get("COAL_ENGINE") or "auto").strip().lower() == "python":
        return None
    if env.get("COAL_CPP"):
        return env["COAL_CPP"]
    try:
        from coal import dispatch
        return dispatch.find_cpp()
    except Exception:  # noqa: BLE001  (no coal on the path of this process: nothing to ask)
        return None


def check_command(cmd, env) -> None:
    """Before a command is started with COAL_TEST_* hook variables: skip the test (unittest.SkipTest), or fail it
    under COAL_REQUIRE_CPP, when the coal-cpp it would run has no test hooks."""
    binary = engine_binary(cmd, env)
    if binary is None or not os.path.isfile(binary) or has_test_hooks(binary):
        return
    if os.environ.get("COAL_REQUIRE_CPP") and os.environ.get("COAL_ALLOW_NO_TEST_HOOKS") != "1":
        raise AssertionError(f"{binary}: {REASON} (COAL_REQUIRE_CPP is set)")
    raise unittest.SkipTest(f"{binary}: {REASON}")
