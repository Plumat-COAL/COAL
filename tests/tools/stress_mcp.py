#!/usr/bin/env python3
"""Drive COAL's stdio MCP server (python -m coal_app --mcp) for the stress test: list the tools, call codon_align.

    python tests/tools/stress_mcp.py --python .venv-mcp/bin/python --list
    python tests/tools/stress_mcp.py --python .venv-mcp/bin/python --call codon_align '{"protein_alignment": "...", ...}'

Speaks JSON-RPC 2.0 over stdio (newline-delimited, as the MCP stdio transport does): initialize, notifications/
initialized, tools/list, tools/call. Prints the response, the wall time and the server's peak memory (psutil).
Standard library (+ psutil for memory).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import threading
import time

try:
    import psutil
except ImportError:
    psutil = None


class Client:
    def __init__(self, cmd):
        self.p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                  bufsize=1, encoding="utf-8")
        self.n = 0
        self.peak = 0
        self.err_lines = []
        threading.Thread(target=self._drain, daemon=True).start()
        if psutil:
            threading.Thread(target=self._mem, daemon=True).start()

    def _drain(self):
        for line in self.p.stderr:
            self.err_lines.append(line.rstrip()[:300])

    def _mem(self):
        proc = psutil.Process(self.p.pid)
        while self.p.poll() is None:
            try:
                rss = proc.memory_info().rss + sum(c.memory_info().rss for c in proc.children(recursive=True))
                self.peak = max(self.peak, rss)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
            time.sleep(0.05)

    def send(self, method, params=None, notify=False):
        msg = {"jsonrpc": "2.0", "method": method}
        if params is not None:
            msg["params"] = params
        if not notify:
            self.n += 1
            msg["id"] = self.n
        self.p.stdin.write(json.dumps(msg) + "\n")
        self.p.stdin.flush()
        if notify:
            return None
        while True:
            line = self.p.stdout.readline()
            if not line:
                raise RuntimeError(f"server closed stdout; stderr: {self.err_lines[-5:]}")
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if r.get("id") == self.n:
                return r

    def close(self):
        try:
            self.p.stdin.close()
            self.p.wait(timeout=10)
        except Exception:
            self.p.kill()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--call", nargs=2, metavar=("TOOL", "JSON"))
    ap.add_argument("--timeout", type=float, default=600)
    a = ap.parse_args(argv)
    c = Client([a.python, "-m", "coal_app", "--mcp"])
    t0 = time.monotonic()
    r = c.send("initialize", {"protocolVersion": "2024-11-05", "capabilities": {},
                              "clientInfo": {"name": "coal-stress", "version": "0"}})
    print("initialize:", json.dumps(r.get("result", r))[:300])
    c.send("notifications/initialized", {}, notify=True)
    if a.list:
        r = c.send("tools/list", {})
        for t in r.get("result", {}).get("tools", []):
            print("TOOL", t["name"], json.dumps(t.get("inputSchema", {}))[:1500])
    if a.call:
        tool, args = a.call[0], json.loads(a.call[1])
        t1 = time.monotonic()
        r = c.send("tools/call", {"name": tool, "arguments": args})
        dt = time.monotonic() - t1
        out = json.dumps(r.get("result", r), indent=None)
        print(f"call {tool}: {dt:.2f}s peak {c.peak / 1e6:.0f} MB")
        print(out[:3000])
    c.close()
    print(f"total {time.monotonic() - t0:.2f}s, server stderr lines: {len(c.err_lines)}")
    for line in c.err_lines[-5:]:
        print("  stderr:", line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
