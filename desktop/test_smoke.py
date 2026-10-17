"""Automated smoke test for the COAL desktop app: start the wrapper's local app WITHOUT the
window, run every example dataset through it, check the results, shut down. All with timeouts.

    python -m unittest desktop/test_smoke.py -v                 # source checkout
    set COAL_SMOKE_EXE=desktop\\dist\\COAL\\COAL.exe              # (PowerShell: $env:COAL_SMOKE_EXE=...)
    python -m unittest desktop/test_smoke.py -v                 # against the PyInstaller build

Layers (each skips cleanly when its requirements are missing):

1. ExamplesViaEngineDispatch: the GUI's own logic (coal_app.tools.codon_align) on all four
   examples, with the engine started the way a frozen build starts it (`<wrapper> -m coal`, or
   `COAL.exe -m coal`). Needs only the standard library. Runs anywhere.
2. WindowsFileQuirks: CRLF files in a folder with spaces (like OneDrive), and a UTF-8 BOM.
3. StandInLocalApp: the same route with a stand-in for Gradio's HTTP layer (runs anywhere).
   RealLocalApp: the real local app (Gradio) started by LocalServer exactly as the window
   starts it, each example sent through its HTTP API (/gradio_api/call/codon_align), then a
   clean shutdown. Needs Gradio (pip install -e .[gui]); skipped otherwise.
4. WindowApi: pywebview exposes the documented names coal_desktop.run_window uses. Needs
   pywebview; skipped otherwise. (The window itself is checked by hand: SMOKE_TEST.md.)
"""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import coal_desktop as cd  # noqa: E402
from test_wrapper import assert_loopback_only, assert_port_free  # noqa: E402

cd.setup_source_paths()  # this repo's src/ first, for us and for child processes

from coal_app import core, tools  # noqa: E402

SMOKE_EXE = os.environ.get("COAL_SMOKE_EXE")
HAS_GRADIO = importlib.util.find_spec("gradio") is not None
HAS_WEBVIEW = importlib.util.find_spec("webview") is not None
EXAMPLES = ("clean", "mito", "messy_ids", "one_broken_gene")
RUN_TIMEOUT = 120  # seconds per engine run


def engine_cmd() -> list[str]:
    """How a frozen build runs the engine: the app's own executable with `-m coal`."""
    if SMOKE_EXE:
        return [SMOKE_EXE, "-m", "coal"]
    return [sys.executable, os.path.join(HERE, "coal_desktop.py"), "-m", "coal"]


def server_cmd(port: int) -> list[str]:
    if SMOKE_EXE:
        return [SMOKE_EXE, "--serve", "--port", str(port), "--parent-pid", str(os.getpid())]
    return cd.default_server_command(port)


def read_fasta(path):
    seqs, name = {}, None
    with open(path, encoding="ascii") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif line and name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


class Checks:
    """Expected outcome per example; valid for today's engine AND for the engine branch
    (--keep-going, ID normalisation), so the test doesn't break when that merges."""

    def check_codon_alignment(self, path, n_expected=None):
        seqs = read_fasta(path)
        self.assertTrue(seqs, path)
        lengths = {len(s) for s in seqs.values()}
        self.assertEqual(len(lengths), 1, "rows of a codon alignment must have equal length")
        self.assertEqual(lengths.pop() % 3, 0, "codon alignment length must be a multiple of 3")
        if n_expected is not None:
            self.assertEqual(len(seqs), n_expected)

    def check_result(self, name, code, r):
        s = r.get("summary", {})
        hints = r.get("hints") or []
        problems = {p["id"]: p["status"] for p in r.get("problems", [])}
        if name == "clean":
            self.assertTrue(r["completed"], r.get("engine_messages"))
            self.assertEqual((s.get("ok"), s.get("warning"), s.get("failed")), (8, 0, 0))
            self.check_codon_alignment(r["output_path"], 8)
            self.assertIn("COAL", r["methods_text"])
            self.assertTrue(r["coal_command"].startswith("coal "))
        elif name == "mito" and code == 1:
            self.assertGreater(s.get("warning", 0), 0)
            self.assertIn(2, [h.get("suggested_codontable") for h in hints], hints)
        elif name == "mito" and code == 2:
            self.assertTrue(r["completed"], r.get("engine_messages"))
            self.assertEqual((s.get("ok"), s.get("warning"), s.get("failed")), (24, 0, 0))
            self.check_codon_alignment(r["output_path"], 24)
            self.assertIn("-codontable 2", r["coal_command"])
        elif name == "messy_ids":
            if r["completed"]:  # engine with ID normalisation
                self.assertEqual(s.get("ok"), 10)
            else:        # today: stops, but explains the ID mismatch
                self.assertIn("id_mismatch", [h.get("code") for h in hints], hints)
        elif name == "one_broken_gene":
            self.assertEqual(problems.get("OG0001_sp14"), "failed", problems)
            self.assertEqual(len([p for p in problems.values() if p == "failed"]), 1)
            if r["completed"]:  # engine with --keep-going: the other 19 are written
                self.assertEqual(s.get("ok"), 19)
                self.check_codon_alignment(r["output_path"], 19)
        self.assertTrue(os.path.isfile(r["report_path"]))

    RUNS = (("clean", 1), ("mito", 1), ("mito", 2), ("messy_ids", 1), ("one_broken_gene", 1))


@unittest.skipUnless(all(os.path.isfile(p) for n in EXAMPLES for p in core.example_paths(n)),
                     "examples/ not found")
class ExamplesViaEngineDispatch(unittest.TestCase, Checks):
    def setUp(self):
        self._cmd, self._timeout = core._coal_cmd, core.ENGINE_TIMEOUT
        core._coal_cmd = engine_cmd
        core.ENGINE_TIMEOUT = RUN_TIMEOUT
        core.detect_features.cache_clear()

    def tearDown(self):
        core._coal_cmd, core.ENGINE_TIMEOUT = self._cmd, self._timeout
        core.detect_features.cache_clear()

    def test_engine_is_reachable(self):
        f = core.detect_features()
        self.assertIn("fasta", f.output_formats)
        self.assertTrue(f.version and f.version != "unknown", f)

    def test_every_example(self):
        for name, code in self.RUNS:
            with self.subTest(example=name, genetic_code=code):
                p, d = core.example_paths(name)
                self.check_result(name, code, tools.codon_align(p, d, "fasta", code))


class WindowsFileQuirks(unittest.TestCase):
    """Files as Windows users have them: CRLF line endings, a path with spaces, a BOM."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="coal smoke ")
        self.dir = os.path.join(self.tmp, "OneDrive - Uni of X", "Thesis data")
        os.makedirs(self.dir)
        self.env = cd.child_env()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def copy(self, src, name, crlf=True, bom=False):
        with open(src, "rb") as f:
            data = f.read()
        if crlf:
            data = data.replace(b"\r\n", b"\n").replace(b"\n", b"\r\n")
        dst = os.path.join(self.dir, name)
        with open(dst, "wb") as f:
            f.write((b"\xef\xbb\xbf" if bom else b"") + data)
        return dst

    def run_engine(self, p, d, out):
        return subprocess.run(engine_cmd() + [p, d, "-o", out], capture_output=True, text=True,
                              env=self.env, timeout=RUN_TIMEOUT)

    def test_crlf_and_spaces_give_identical_output(self):
        p, d = core.example_paths("clean")
        ref = os.path.join(self.tmp, "ref.fasta")
        self.assertEqual(self.run_engine(p, d, ref).returncode, 0)
        out = os.path.join(self.dir, "my result.fasta")
        r = self.run_engine(self.copy(p, "my protein.aln.fasta"), self.copy(d, "my cds.fasta"), out)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(ref, "rb") as a, open(out, "rb") as b:
            self.assertEqual(a.read(), b.read())

    # CHANGED(INT-2): was @unittest.expectedFailure; the engine now skips a UTF-8 BOM (ENG-h), so this must pass.
    def test_utf8_bom_is_accepted(self):
        p, d = core.example_paths("clean")
        r = self.run_engine(self.copy(p, "bom.aln.fasta", bom=True), d, os.path.join(self.dir, "o.fasta"))
        self.assertEqual(r.returncode, 0, r.stderr)


def http(method, url, body=None, timeout=60):
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))  # never via a proxy
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(url, data=data, method=method,
                                 headers={"Content-Type": "application/json"} if data else {})
    with opener.open(req, timeout=timeout) as r:
        return r.status, r.read().decode("utf-8", "replace")


def gradio_call(base, api_name, args, timeout=RUN_TIMEOUT):
    """Gradio 5 HTTP API: POST /gradio_api/call/<name> -> event_id; GET .../<event_id> -> SSE.
    (From Gradio's 'Querying Gradio Apps with Curl' guide; untested here, no Gradio.)"""
    _, txt = http("POST", f"{base}gradio_api/call/{api_name}", {"data": args}, timeout=30)
    event_id = json.loads(txt)["event_id"]
    _, sse = http("GET", f"{base}gradio_api/call/{api_name}/{event_id}", timeout=timeout)
    event = None
    for line in sse.splitlines():
        if line.startswith("event:"):
            event = line.split(":", 1)[1].strip()
        elif line.startswith("data:") and event in ("complete", "error"):
            payload = line.split(":", 1)[1].strip()
            if event == "error":
                raise AssertionError(f"Gradio API error for {api_name}: {payload}")
            return json.loads(payload)[0]
    raise AssertionError(f"no result in Gradio API response: {sse[:500]}")


# WINCI-2: the stand-in now starts exactly like the product's child (`coal_desktop.py --serve --port N --parent-pid P`
# -> main() -> dispatch() -> serve()), with only launch_gradio replaced; and like Gradio's queue, POST returns the
# event id at once and the run happens in a worker thread (GET waits for it).
STAND_IN = r"""
import json, sys, threading, uuid, http.server
sys.path.insert(0, %(desktop)r)
import coal_desktop as cd
cd.setup_source_paths()
from coal_app import tools

JOBS = {}

class H(http.server.BaseHTTPRequestHandler):
    # Stand-in for Gradio's HTTP API (same URLs and SSE shape), calling the same tool function (ui.py: gr.api).
    def _send(self, code, body, ctype="application/json"):
        b = body.encode()
        self.send_response(code); self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(b))); self.end_headers(); self.wfile.write(b)
    def do_GET(self):
        if self.path == "/":
            return self._send(200, "<title>COAL codon alignment</title>", "text/html")
        eid = self.path.rsplit("/", 1)[-1]
        if self.path.startswith("/gradio_api/call/codon_align/") and eid in JOBS:
            done, box = JOBS.pop(eid)
            done.wait()
            event, data = box[0]
            return self._send(200, "event: " + event + "\ndata: " + json.dumps(data) + "\n\n", "text/event-stream")
        self._send(404, "{}")
    def do_POST(self):
        if self.path != "/gradio_api/call/codon_align":
            return self._send(404, "{}")
        args = json.loads(self.rfile.read(int(self.headers["Content-Length"])))["data"]
        eid, done, box = uuid.uuid4().hex, threading.Event(), []
        def job():
            try:
                box.append(("complete", [tools.codon_align(*args)]))
            except Exception as e:  # reported to the test like Gradio's "event: error"
                box.append(("error", repr(e)))
            finally:
                done.set()
        JOBS[eid] = (done, box)
        threading.Thread(target=job, daemon=True).start()
        self._send(200, json.dumps({"event_id": eid}))
    def log_message(self, *a): pass

def launcher(port, mcp):
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", port), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return "http://127.0.0.1:%%d/" %% port, lambda: (srv.shutdown(), srv.server_close())

cd.launch_gradio = launcher  # the only difference from the product's --serve child
sys.exit(cd.main(sys.argv[1:]))
"""


class LocalAppTests(Checks):
    """Start the local app with the wrapper's own LocalServer, use it, stop it."""

    @classmethod
    def command(cls, port):
        return server_cmd(port)

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="coal_smoke_")
        cls.server = cd.LocalServer(command=cls.command, start_timeout=120,
                                    log_path=os.path.join(cls.tmp, "local_app.log"))
        try:
            cls.url = cls.server.start()
        except Exception:
            cls.server.stop()
            raise

    @classmethod
    def tearDownClass(cls):
        code = cls.server.stop()
        shutil.rmtree(cls.tmp, ignore_errors=True)
        if code not in (0, None):
            raise AssertionError(f"local app exit code {code} on shutdown (expected 0)")

    def test_1_listens_only_on_loopback(self):
        self.assertTrue(self.url.startswith("http://127.0.0.1:"), self.url)
        assert_loopback_only(self, self.server.port)  # WIN-W5: the socket itself, not the URL
        status, page = http("GET", self.url, timeout=30)
        self.assertEqual(status, 200)
        self.assertIn("COAL", page)

    def test_2_every_example_through_the_server(self):
        for name, code in self.RUNS:
            with self.subTest(example=name, genetic_code=code):
                p, d = core.example_paths(name)
                # CHANGED(SHIP-10): was 6 values; codon_align has 8 parameters since FIN-A6 added strict and blockonly, and
                # the real Gradio API rejects a short positional list (the first CI build of the installer failed here)
                r = gradio_call(self.url, "codon_align", [p, d, "fasta", code, False, False, False, False, ""])
                self.check_result(name, code, r)

    def test_3_shutdown_frees_the_port(self):
        port = self.server.port
        self.assertEqual(self.server.stop(), 0)
        self.assertFalse(self.server.running())
        assert_port_free(self, port)  # WIN-W5


@unittest.skipUnless(HAS_GRADIO or SMOKE_EXE, "Gradio not installed (pip install -e .[gui])")
class RealLocalApp(LocalAppTests, unittest.TestCase):
    """The real local app (Gradio, or COAL.exe --serve)."""


class StandInLocalApp(LocalAppTests, unittest.TestCase):
    """Same process route (LocalServer -> serve() child -> coal_app.tools -> engine run), with a
    stdlib stand-in for Gradio's HTTP layer. Runs anywhere; checks everything except Gradio."""

    @classmethod
    def command(cls, port):
        path = os.path.join(cls.tmp, "stand_in_app.py")
        if not os.path.exists(path):
            with open(path, "w", encoding="utf-8") as f:
                f.write(STAND_IN % {"desktop": HERE})
        return [sys.executable, path] + cd.default_server_command(port)[2:]  # WINCI-2: the product's arguments


@unittest.skipUnless(HAS_WEBVIEW, "pywebview not installed")
class WindowApi(unittest.TestCase):
    """WIN-W5: the pywebview names AND parameters run_window relies on (pywebview 6.x docs).
    This checks signatures only; whether the window behaves is checked by hand (SMOKE_TEST.md)."""

    @staticmethod
    def params(f):
        import inspect
        return list(inspect.signature(f).parameters)

    def test_documented_signatures(self):
        import inspect
        import webview
        from webview.dom import DOMEventHandler
        from webview.window import Window
        cw = self.params(webview.create_window)
        for n in ("title", "html", "width", "height", "min_size", "text_select", "zoomable", "background_color"):
            self.assertIn(n, cw, f"create_window({n}=)")
        st = self.params(webview.start)
        self.assertEqual(st[:2], ["func", "args"])
        for n in ("gui", "private_mode", "storage_path", "debug"):
            self.assertIn(n, st, f"start({n}=)")
        self.assertIn("ALLOW_DOWNLOADS", webview.settings)
        dh = [p for p in self.params(DOMEventHandler.__init__) if p != "self"]
        self.assertEqual(dh[:3], ["callback", "prevent_default", "stop_propagation"])  # passed positionally
        self.assertIn("debounce", dh)
        for name in ("load_url", "get_current_url", "restore", "show", "destroy"):
            self.assertTrue(callable(getattr(Window, name, None)), name)
        self.assertIn("url", self.params(Window.load_url))
        self.assertIsInstance(inspect.getattr_static(Window, "on_top"), property)


if __name__ == "__main__":
    unittest.main()
