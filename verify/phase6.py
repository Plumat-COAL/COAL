#!/usr/bin/env python3
"""Phase 6: the stdio MCP server, driven by a minimal JSON-RPC client (no mcp package needed on the client side).

    python verify/phase6.py OUTDIR            # uses .venv-mcp/bin/python -m coal_app --mcp
"""
from __future__ import annotations

import json
import os
import random
import select
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adversarial import CODE1, fasta, gene, genes, write  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
MCP_PY = ROOT / ".venv-mcp" / "bin" / "python"
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append({"name": name, "ok": bool(ok), "detail": str(detail)[:800]})
    print(("OK   " if ok else "FINDING-CANDIDATE ") + name + ("" if ok else f": {str(detail)[:400]}"), flush=True)


class Client:
    """A stdio MCP client: newline-delimited JSON-RPC, with timeouts on every read."""

    def __init__(self, cmd, cwd, env=None, timeout=60):
        self.p = subprocess.Popen(cmd, cwd=cwd, env=env, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.timeout = timeout
        self.n = 0
        self.stray = []  # stdout bytes that are not JSON-RPC messages

    def send(self, method, params=None, notify=False):
        msg = {"jsonrpc": "2.0", "method": method}
        if params is not None:
            msg["params"] = params
        if not notify:
            self.n += 1
            msg["id"] = self.n
        self.p.stdin.write((json.dumps(msg) + "\n").encode())
        self.p.stdin.flush()
        return None if notify else self.recv(self.n)

    def recv(self, want_id):
        deadline = time.monotonic() + self.timeout
        buf = b""
        while time.monotonic() < deadline:
            r, _, _ = select.select([self.p.stdout], [], [], 0.2)
            if not r:
                if self.p.poll() is not None:
                    raise RuntimeError(f"server exited {self.p.returncode}: {self.p.stderr.read()[-500:]!r}")
                continue
            chunk = os.read(self.p.stdout.fileno(), 65536)
            if not chunk:
                raise RuntimeError("server closed stdout")
            buf += chunk
            while b"\n" in buf:
                line, buf = buf.split(b"\n", 1)
                if not line.strip():
                    continue
                try:
                    m = json.loads(line)
                except ValueError:
                    self.stray.append(line[:200])
                    continue
                if m.get("id") == want_id:
                    return m
                # notifications / other ids are ignored
        raise TimeoutError(f"no reply to {want_id} within {self.timeout} s")

    def call(self, name, args, timeout=None):
        if timeout:
            self.timeout = timeout
        return self.send("tools/call", {"name": name, "arguments": args})

    def close(self):
        try:
            self.p.stdin.close()
        except OSError:
            pass
        try:
            self.p.wait(timeout=10)
        except subprocess.TimeoutExpired:
            self.p.kill()
            self.p.wait(timeout=10)
        return self.p.stderr.read()


def result_of(resp):
    """The structured result of a tools/call reply: structuredContent, else the first text content parsed."""
    r = resp.get("result") or {}
    if r.get("isError"):
        return {"isError": True, "text": " ".join(c.get("text", "") for c in r.get("content", []))}
    if "structuredContent" in r:
        sc = r["structuredContent"]
        if isinstance(sc, dict) and set(sc) == {"result"}:  # a tool returning a plain string (FastMCP wraps it)
            return {"text": sc["result"]}
        return sc
    for c in r.get("content", []):
        if c.get("type") == "text":
            try:
                return json.loads(c["text"])
            except ValueError:
                return {"text": c["text"]}
    return r


def main(outdir):
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(66)
    p, n = genes(rng, 6, 40)
    write(out, "pep.fasta", fasta(p)); write(out, "nuc.fasta", fasta(n))
    j, s = n[2]; n[2] = (j, s[:30] + CODE1["W" if p[2][1][10] != "W" else "A"] + s[33:])
    write(out, "nuc_warn.fasta", fasta(n))
    nb = list(n); j, s = nb[4]; nb[4] = (j, s[:60] + "G" + s[60:])
    write(out, "nuc_broken.fasta", fasta(nb))
    aas2 = dict(CODE1, W="TGA", M="ATA")
    p2, n2 = [], []
    for k in range(5):
        pp, ss = gene(rng, 40, aas2); p2.append((f"m{k}", pp)); n2.append((f"m{k}", ss))
    write(out, "pep_mito.fasta", fasta(p2)); write(out, "nuc_mito.fasta", fasta(n2))
    sp = out / "with space"; sp.mkdir(exist_ok=True)
    write(sp, "pep.fasta", fasta(p)); write(sp, "nuc.fasta", (out / "nuc.fasta").read_text())
    big_p, big_n = genes(rng, 12_000, 40)   # about 1.5 MB of DNA: the coal-cpp hand-over
    write(out, "big_p.fasta", fasta(big_p)); write(out, "big_n.fasta", fasta(big_n))
    env = dict(os.environ, COAL_RUNS_DIR=str(out / "runs"), COAL_CPP=str(ROOT / "build" / "coal-cpp"))
    c = Client([str(MCP_PY), "-m", "coal_app", "--mcp"], cwd=str(out), env=env)
    try:
        init = c.send("initialize", {"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "verify", "version": "1"}})
        check("MCP initialize answers with serverInfo", "result" in init and "serverInfo" in init["result"], json.dumps(init)[:300])
        info = init.get("result", {}).get("serverInfo", {})
        check("MCP serverInfo.version is COAL's version (0.5.0)", info.get("version") == "0.5.0", info)
        c.send("notifications/initialized", notify=True)
        tl = c.send("tools/list")
        tools = {t["name"]: t for t in tl.get("result", {}).get("tools", [])}
        check("MCP tools/list: exactly codon_align and explain_report", set(tools) == {"codon_align", "explain_report"}, sorted(tools))
        (out / "tools_list.json").write_text(json.dumps(tl, indent=1), encoding="utf-8")
        ca = tools.get("codon_align", {})
        props = ca.get("inputSchema", {}).get("properties", {})
        check("codon_align schema: the nine parameters with descriptions",
              set(props) == {"protein_alignment_path", "dna_fasta_path", "output_format", "codontable", "nogap", "nomismatch", "strict", "blockonly", "save_as"}
              and all(pr.get("description") for pr in props.values()), sorted(props))
        check("codon_align: codontable limited to the 17 codes", sorted(props.get("codontable", {}).get("enum", [])) == [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23], props.get("codontable"))
        check("codon_align: output_format enum of the 4 formats", sorted(props.get("output_format", {}).get("enum", [])) == ["clustal", "codon", "fasta", "paml"], props.get("output_format"))
        check("codon_align has a description (docstring) mentioning all_ok / needs_attention", "all_ok" in ca.get("description", "") and "needs_attention" in ca.get("description", ""), ca.get("description", "")[:200])
        er = tools.get("explain_report", {})
        check("explain_report schema: report_path only, described", set(er.get("inputSchema", {}).get("properties", {})) == {"report_path"}, er.get("inputSchema"))

        def align(**kw):
            base = {"protein_alignment_path": str(out / "pep.fasta"), "dna_fasta_path": str(out / "nuc.fasta")}
            base.update(kw)
            return result_of(c.call("codon_align", base, timeout=600))
        r = align()
        check("codon_align clean: completed, all_ok, counts ok=6", r.get("completed") is True and r.get("all_ok") is True and r.get("needs_attention") is False and r.get("counts", {}).get("ok") == 6, r)
        check("codon_align clean: output_path and report_path exist", Path(r.get("output_path", "/nonexistent")).is_file() and Path(r.get("report_path", "/nonexistent")).is_file(), (r.get("output_path"), r.get("report_path")))
        cli = subprocess.run([str(ROOT / ".venv" / "bin" / "python"), "-m", "coal", str(out / "pep.fasta"), str(out / "nuc.fasta"), "--keep-going"], capture_output=True, timeout=120)
        check("codon_align clean: alignment byte-identical to the CLI", Path(r["output_path"]).read_bytes() == cli.stdout)
        clean_report = r["report_path"]
        r = align(dna_fasta_path=str(out / "nuc_warn.fasta"))
        check("codon_align with one mismatch: needs_attention, counts warning=1, problems named", r.get("needs_attention") is True and r.get("counts", {}).get("warning") == 1 and r.get("problems"), r)
        r = align(dna_fasta_path=str(out / "nuc_broken.fasta"))
        check("codon_align with one broken gene: completed (keep-going), failed=1, needs_attention, next_step says what to do", r.get("completed") is True and r.get("counts", {}).get("failed") == 1 and r.get("needs_attention") is True and r.get("next_step"), r)
        broken_report = r.get("report_path")
        r = align(protein_alignment_path=str(out / "pep_mito.fasta"), dna_fasta_path=str(out / "nuc_mito.fasta"))
        check("codon_align mito under code 1: suggested_codontable == 2", r.get("suggested_codontable") == 2 and r.get("needs_attention") is True, r)
        mito_report = r.get("report_path")
        r = align(protein_alignment_path=str(out / "pep_mito.fasta"), dna_fasta_path=str(out / "nuc_mito.fasta"), codontable=2)
        check("codon_align mito with code 2: all_ok", r.get("all_ok") is True, r)
        for fmt in ("clustal", "paml", "codon"):
            r = align(output_format=fmt)
            cli = subprocess.run([str(ROOT / ".venv" / "bin" / "python"), "-m", "coal", str(out / "pep.fasta"), str(out / "nuc.fasta"), "--keep-going", "-output", fmt], capture_output=True, timeout=120)
            check(f"codon_align output_format={fmt}: identical to the CLI", r.get("all_ok") and Path(r["output_path"]).read_bytes() == cli.stdout, r.get("output_path"))
        for flag in ("nogap", "nomismatch", "strict", "blockonly"):
            r = align(**{flag: True}, dna_fasta_path=str(out / "nuc_warn.fasta"))
            cli = subprocess.run([str(ROOT / ".venv" / "bin" / "python"), "-m", "coal", str(out / "pep.fasta"), str(out / "nuc_warn.fasta"), "--keep-going", "-" + flag if flag != "strict" else "--strict"], capture_output=True, timeout=120)
            same = Path(r.get("output_path", "/nonexistent")).exists() and Path(r["output_path"]).read_bytes() == cli.stdout
            check(f"codon_align {flag}=true: same alignment as the CLI with that option", same, (r.get("completed"), r.get("counts"), r.get("output_path")))
        # bad inputs are tool errors
        r = align(protein_alignment_path=str(out / "missing.fasta"))
        check("codon_align missing file: tool error (isError)", r.get("isError") and "missing.fasta" in r.get("text", ""), r)
        r = align(protein_alignment_path="pep.fasta")
        check("codon_align relative path: resolved against the server's folder or a clear error (docs say absolute)", r.get("isError") or r.get("all_ok") is True, r)
        r = align(protein_alignment_path=str(sp / "pep.fasta"), dna_fasta_path=str(sp / "nuc.fasta"))
        check("codon_align paths with spaces: ok", r.get("all_ok") is True, r)
        r = align(protein_alignment_path=str(out / "nuc.fasta"), dna_fasta_path=str(out / "pep.fasta"))
        check("codon_align files swapped: a clear error or a failed result, never all_ok", not r.get("all_ok"), r)
        r = result_of(c.call("codon_align", {"protein_alignment_path": str(out / "pep.fasta"), "dna_fasta_path": str(out / "nuc.fasta"), "codontable": 7}))
        check("codon_align codontable 7 (unknown): tool / validation error", r.get("isError") or "error" in json.dumps(r).lower(), r)
        r = result_of(c.call("codon_align", {"protein_alignment_path": str(out / "pep.fasta"), "dna_fasta_path": str(out / "nuc.fasta"), "output_format": "nexus"}))
        check("codon_align output_format nexus: tool / validation error", r.get("isError") or "error" in json.dumps(r).lower(), r)
        # save_as
        target = out / "saved.fasta"
        r = align(save_as=str(target))
        check("codon_align save_as new file: written there", target.is_file() and r.get("output_path") == str(target), r.get("output_path"))
        before = target.read_bytes()
        write(out, "saved.fasta", "DO NOT OVERWRITE\n")
        r = align(save_as=str(target))
        check("codon_align save_as existing file: not overwritten (documented), result says where it went or errors",
              target.read_text() == "DO NOT OVERWRITE\n", (r.get("output_path"), r.get("isError"), r.get("text", "")[:200]))
        folder = out / "savefolder"; folder.mkdir(exist_ok=True)
        r = align(save_as=str(folder))
        check("codon_align save_as existing folder: saved inside it", r.get("output_path", "").startswith(str(folder) + "/") and Path(r["output_path"]).is_file(), r.get("output_path"))
        r = align(save_as=str(out / "no_such_dir" / "x.fasta"))
        check("codon_align save_as into a missing folder: tool error or clear message, no crash", r.get("isError") or r.get("completed") is not None, r)
        # large input -> coal-cpp hand-over through the app
        t0 = time.monotonic()
        r = align(protein_alignment_path=str(out / "big_p.fasta"), dna_fasta_path=str(out / "big_n.fasta"))
        big_time = time.monotonic() - t0
        rep = json.loads(Path(r["report_path"]).read_text()) if r.get("report_path") else {}
        check("codon_align large input (>= 1 MB): completed, all_ok, report says which engine ran", r.get("all_ok") is True and rep.get("engine") in ("coal-cpp", "coal-python"), (r.get("counts"), rep.get("engine"), f"{big_time:.1f}s"))
        check("codon_align large input: engine is coal-cpp (hand-over through the app layer)", rep.get("engine") == "coal-cpp", rep.get("engine"))
        # explain_report
        for name, path in (("clean", clean_report), ("broken", broken_report), ("mito", mito_report)):
            r = result_of(c.call("explain_report", {"report_path": path}))
            txt = r.get("text") if isinstance(r, dict) else str(r)
            check(f"explain_report {name}: plain text, mentions counts", isinstance(txt, str) and len(txt) > 40 and not r.get("isError"), str(txt)[:300])
            (out / f"explain_{name}.txt").write_text(str(txt), encoding="utf-8")
        r = result_of(c.call("explain_report", {"report_path": str(out / "pep.fasta")}))
        check("explain_report on a non-report file: tool error", r.get("isError"), r)
        r = result_of(c.call("explain_report", {"report_path": str(out / "nope.json")}))
        check("explain_report missing file: tool error", r.get("isError"), r)
    finally:
        err = c.close()
    check("MCP server wrote nothing to stdout outside the protocol", not c.stray, c.stray[:3])
    check("MCP server stderr holds no tracebacks", b"Traceback" not in err, err[-500:])
    (out / "server_stderr.txt").write_bytes(err)
    (out / "results.json").write_text(json.dumps(CHECKS, indent=1), encoding="utf-8")
    bad = [x for x in CHECKS if not x["ok"]]
    print(f"\n{len(CHECKS)} checks, {len(bad)} finding candidates")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "phase6_out")
