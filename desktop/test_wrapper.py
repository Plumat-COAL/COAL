"""Tests for the window-independent parts of desktop/coal_desktop.py (standard library only).

Run from the repo root (any OS, no pywebview or Gradio needed):

    python -m unittest desktop/test_wrapper.py -v

Every wait has a timeout; every process started here is stopped in tearDown/finally.
The local app is replaced by a tiny stand-in HTTP server that speaks the same protocol
(``coal_desktop.serve`` with a stdlib launcher), so the real port-picking, readiness,
retry, lifeline and shutdown code is exercised.
"""
from __future__ import annotations

import io
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import coal_desktop as cd  # noqa: E402

STUB = textwrap.dedent(r'''
    import http.server, os, socket, sys, threading, time
    sys.path.insert(0, {desktop!r})
    import coal_desktop as cd

    mode = os.environ.get("STUB_MODE", "normal")
    watch = os.environ.get("STUB_WATCH_PORT")
    if watch:  # the test sees this connection close when this process dies
        _keep = socket.create_connection(("127.0.0.1", int(watch)), timeout=10)

    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            body = b"stub coal app"
            self.send_response(200); self.send_header("Content-Length", str(len(body)))
            self.end_headers(); self.wfile.write(body)
        def log_message(self, *a): pass

    def launcher(port, mcp):
        if mode == "hang":
            time.sleep(600)
        if mode == "crash":
            raise RuntimeError("stub crashed on purpose")
        if mode == "nogui":
            raise ImportError("no module named gradio", name="gradio")
        srv = http.server.ThreadingHTTPServer(("127.0.0.1", port), H)
        threading.Thread(target=srv.serve_forever, daemon=True).start()
        def close():
            srv.shutdown(); srv.server_close()
            if os.environ.get("STUB_CLOSED_FILE"):
                open(os.environ["STUB_CLOSED_FILE"], "w").close()
        return f"http://127.0.0.1:{{port}}/", close

    port = int(sys.argv[1])
    sys.exit(cd.serve(port, launcher=launcher))
''').format(desktop=HERE)

PARENT = textwrap.dedent(r'''
    import sys, time
    sys.path.insert(0, {desktop!r})
    import coal_desktop as cd
    s = cd.LocalServer(command=lambda p: [sys.executable, {stub!r}, str(p)], start_timeout=30)
    s.start()
    print(s.proc.pid, flush=True)
    time.sleep(600)  # killed by the test
''')


def http_get(url: str, timeout: float = 10) -> tuple[int, bytes]:
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))  # never via a proxy
    with opener.open(url, timeout=timeout) as r:
        return r.status, r.read()


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="coal_desktop_test_")
        self.stub = os.path.join(self.tmp, "stub_server.py")
        with open(self.stub, "w", encoding="utf-8") as f:
            f.write(STUB)
        self.env = cd.child_env()
        self.servers = []

    def tearDown(self):
        for s in self.servers:
            s.stop()
        shutil.rmtree(self.tmp, ignore_errors=True)

    def server(self, mode="normal", **kw) -> cd.LocalServer:
        env = dict(self.env, STUB_MODE=mode, **{k: str(v) for k, v in kw.pop("env", {}).items()})
        s = cd.LocalServer(command=kw.pop("command", lambda p: [sys.executable, self.stub, str(p)]),
                           env=env, log_path=os.path.join(self.tmp, "local_app.log"), **kw)
        self.servers.append(s)
        return s


class TestPorts(unittest.TestCase):
    def test_free_port_is_bindable(self):
        p = cd.find_free_port()
        self.assertTrue(1024 <= p <= 65535)
        with socket.socket() as s:
            s.bind(("127.0.0.1", p))

    def test_port_busy_detection(self):
        import errno
        self.assertTrue(cd.is_port_busy_error(OSError(errno.EADDRINUSE, "x")))
        self.assertTrue(cd.is_port_busy_error(OSError("Cannot find empty port in range: 5-5.")))
        self.assertTrue(cd.is_port_busy_error(OSError(10048, "Only one usage of each socket address "
                                                             "(protocol/network address/port) is normally permitted")))
        self.assertFalse(cd.is_port_busy_error(OSError(errno.ENOENT, "No such file")))

    def test_real_bind_error_is_detected(self):
        with socket.socket() as a:
            a.bind(("127.0.0.1", 0))
            a.listen(1)
            with socket.socket() as b:
                try:
                    b.bind(("127.0.0.1", a.getsockname()[1]))
                except OSError as e:
                    self.assertTrue(cd.is_port_busy_error(e), e)
                else:
                    self.skipTest("OS allowed a second bind")


class TestLocalServer(Base):
    def test_start_serve_stop(self):
        closed = os.path.join(self.tmp, "closed")
        s = self.server(env={"STUB_CLOSED_FILE": closed}, start_timeout=30)
        url = s.start()
        self.assertTrue(url.startswith("http://127.0.0.1:"), url)
        self.assertEqual(url, f"http://127.0.0.1:{s.port}/")
        status, body = http_get(url)
        self.assertEqual((status, body), (200, b"stub coal app"))
        t0 = time.monotonic()
        code = s.stop()
        self.assertEqual(code, 0)  # clean exit: stdin closed -> serve() closed the server
        self.assertLess(time.monotonic() - t0, cd.SERVER_STOP_TIMEOUT)
        self.assertTrue(os.path.exists(closed), "close() of the app was not called")
        self.assertFalse(s.running())
        with socket.socket() as sock:  # listener gone (REUSEADDR: ignore TIME_WAIT, not LISTEN)
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            sock.bind(("127.0.0.1", s.port))
        self.assertEqual(s.stop(), 0)  # second stop is harmless
        with open(os.path.join(self.tmp, "local_app.log"), encoding="utf-8") as f:
            self.assertIn(cd.READY_PREFIX, f.read())

    def test_retries_when_port_is_taken(self):
        busy = socket.socket()
        busy.bind(("127.0.0.1", 0))
        busy.listen(1)
        calls = []

        def command(port):
            calls.append(port)
            use = busy.getsockname()[1] if len(calls) == 1 else port  # first try: a taken port
            return [sys.executable, self.stub, str(use)]
        try:
            s = self.server(command=command, start_timeout=30)
            url = s.start()
            self.assertEqual(len(calls), 2)
            self.assertEqual(http_get(url)[0], 200)
        finally:
            busy.close()

    def test_gives_up_after_three_busy_ports(self):
        busy = socket.socket()
        busy.bind(("127.0.0.1", 0))
        busy.listen(1)
        calls = []

        def command(port):
            calls.append(port)
            return [sys.executable, self.stub, str(busy.getsockname()[1])]
        try:
            s = self.server(command=command, start_timeout=30)
            with self.assertRaises(cd.ServerError):
                s.start()
            self.assertEqual(len(calls), cd.PORT_ATTEMPTS)
        finally:
            busy.close()

    def test_start_timeout_kills_the_child(self):
        s = self.server(mode="hang", start_timeout=2)
        t0 = time.monotonic()
        with self.assertRaises(cd.ServerError) as cm:
            s.start()
        self.assertLess(time.monotonic() - t0, 15)
        self.assertIn("did not start within 2 seconds", str(cm.exception))
        self.assertFalse(s.running())

    def test_child_crash_is_reported(self):
        s = self.server(mode="crash", start_timeout=30)
        with self.assertRaises(cd.ServerError) as cm:
            s.start()
        self.assertIn("stopped while starting", str(cm.exception))
        self.assertIn("stub crashed on purpose", str(cm.exception))

    def test_missing_gui_is_reported_once(self):
        calls = []

        def command(port):
            calls.append(port)
            return [sys.executable, self.stub, str(port)]
        s = self.server(mode="nogui", command=command, start_timeout=30)
        with self.assertRaises(cd.ServerError) as cm:
            s.start()
        self.assertIn("GUI package is missing", str(cm.exception))
        self.assertEqual(len(calls), 1)  # no pointless retries

    def test_child_exits_when_parent_dies(self):
        """Lifeline: kill the window process hard; the local app must exit by itself."""
        watch = socket.socket()
        watch.bind(("127.0.0.1", 0))
        watch.listen(1)
        watch.settimeout(30)
        parent_py = os.path.join(self.tmp, "parent.py")
        with open(parent_py, "w", encoding="utf-8") as f:
            f.write(PARENT.format(desktop=HERE, stub=self.stub))
        env = dict(self.env, STUB_WATCH_PORT=str(watch.getsockname()[1]))
        parent = subprocess.Popen([sys.executable, parent_py], stdout=subprocess.PIPE, env=env)
        conn = None
        try:
            conn, _ = watch.accept()          # the stub child connected
            child_pid = int(parent.stdout.readline())  # parent says the child is ready
            self.assertGreater(child_pid, 0)
            parent.kill()                      # like Task Manager "End task"
            parent.wait(timeout=10)
            conn.settimeout(15)
            t0 = time.monotonic()
            self.assertEqual(conn.recv(1), b"")  # blocks until the child's socket closes = it exited
            self.assertLess(time.monotonic() - t0, 15)
        finally:
            if parent.poll() is None:
                parent.kill()
                parent.wait(timeout=10)
            parent.stdout.close()
            if conn:
                conn.close()
            watch.close()


class TestServeProtocol(unittest.TestCase):
    def test_ready_line_and_clean_exit_on_eof(self):
        out, closed = io.StringIO(), []
        code = cd.serve(5555, launcher=lambda p, m: (f"http://127.0.0.1:{p}/", lambda: closed.append(1)),
                        stdin=io.BytesIO(b""), stdout=out)
        self.assertEqual(code, 0)
        self.assertEqual(out.getvalue().strip(), "COAL-READY http://127.0.0.1:5555/")
        self.assertEqual(closed, [1])

    def test_port_busy_exit_code(self):
        import errno

        def busy(p, m):
            raise OSError(errno.EADDRINUSE, "Address already in use")
        self.assertEqual(cd.serve(1, launcher=busy, stdin=io.BytesIO(b""), stdout=io.StringIO()),
                         cd.EXIT_PORT_BUSY)

    def test_missing_gradio_exit_code(self):
        def nogui(p, m):
            raise ImportError("No module named 'gradio'", name="gradio")
        out = io.StringIO()
        self.assertEqual(cd.serve(1, launcher=nogui, stdin=io.BytesIO(b""), stdout=out), cd.EXIT_NO_GUI)
        self.assertIn("gradio", out.getvalue())

    def test_default_command_source_mode(self):
        cmd = cd.default_server_command(1234)
        self.assertEqual(cmd[0], sys.executable)
        self.assertTrue(cmd[1].endswith("coal_desktop.py"))
        self.assertEqual(cmd[2:5], ["--serve", "--port", "1234"])
        self.assertIn("--parent-pid", cmd)
        self.assertNotIn("--mcp", cmd)
        self.assertIn("--mcp", cd.default_server_command(1234, mcp=True))

    def test_child_env_offline_and_no_proxy(self):
        old = dict(os.environ)
        try:
            os.environ["NO_PROXY"] = "corp.example"
            os.environ.pop("no_proxy", None)
            env = cd.child_env()
        finally:
            os.environ.clear()
            os.environ.update(old)
        self.assertEqual(env["GRADIO_ANALYTICS_ENABLED"], "False")
        self.assertEqual(env["NO_PROXY"].split(","), ["corp.example", "127.0.0.1", "localhost"])


class TestSingleInstance(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="coal_instance_")

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_second_instance_hands_over_and_first_receives(self):
        import threading
        got, ev = [], threading.Event()
        first = cd.SingleInstance(self.dir, on_message=lambda m: (got.append(m), ev.set()))
        self.assertTrue(first.acquire())
        try:
            second = cd.SingleInstance(self.dir)
            self.assertFalse(second.acquire())
            self.assertTrue(second.notify_primary({"action": "focus", "files": ["C:\\a b\\x.aln"]}))
            self.assertTrue(ev.wait(10))
            self.assertEqual(got, [{"action": "focus", "files": ["C:\\a b\\x.aln"]}])  # token removed
        finally:
            first.release()
        third = cd.SingleInstance(self.dir)
        self.assertTrue(third.acquire())  # free again after release
        third.release()

    def test_wrong_token_is_refused(self):
        got = []
        first = cd.SingleInstance(self.dir, on_message=got.append)
        self.assertTrue(first.acquire())
        try:
            with open(first.info_path, encoding="utf-8") as f:
                info = json.load(f)
            with socket.create_connection(("127.0.0.1", info["port"]), timeout=5) as c:
                c.settimeout(5)
                c.sendall(json.dumps({"token": "wrong", "action": "focus"}).encode() + b"\n")
                self.assertEqual(c.recv(16), b"denied\n")
            with socket.create_connection(("127.0.0.1", info["port"]), timeout=5) as c:
                c.settimeout(5)
                c.sendall(b"not json\n")
                self.assertEqual(c.recv(16), b"")  # dropped, listener still alive:
            self.assertTrue(cd.SingleInstance(self.dir).notify_primary({"action": "focus"}))
        finally:
            first.release()

    def test_lock_is_freed_when_holder_crashes(self):
        code = textwrap.dedent(f'''
            import sys, time
            sys.path.insert(0, {HERE!r})
            import coal_desktop as cd
            i = cd.SingleInstance({self.dir!r})
            print(i.acquire(), flush=True)
            time.sleep(600)
        ''')
        p = subprocess.Popen([sys.executable, "-c", code], stdout=subprocess.PIPE, text=True)
        try:
            self.assertEqual(p.stdout.readline().strip(), "True")
            self.assertFalse(cd.SingleInstance(self.dir).acquire())
            p.kill()
            p.wait(timeout=10)
            again = cd.SingleInstance(self.dir)  # WIN-R5
            # the OS released the dead process's lock (asynchronously on Windows: allow 3 s)
            self.assertEqual(cd.claim_instance(again, {"action": "focus"}, wait=3.0), "primary")
            again.release()
        finally:
            if p.poll() is None:
                p.kill()
                p.wait(timeout=10)
            p.stdout.close()

    def test_notify_without_primary_is_false(self):
        self.assertFalse(cd.SingleInstance(self.dir).notify_primary({"action": "focus"}, timeout=2))


class TestDroppedFiles(unittest.TestCase):
    def test_examples_are_sorted(self):
        for name in ("clean", "mito", "messy_ids", "one_broken_gene"):
            p = os.path.join(REPO, "examples", name, "protein.aln.fasta")
            d = os.path.join(REPO, "examples", name, "dna.fasta")
            self.assertEqual(cd.sniff_sequence_file(p), "protein", name)
            self.assertEqual(cd.sniff_sequence_file(d), "dna", name)
            r = cd.classify_dropped_files([d, p, os.path.join(REPO, "missing.fa")])  # order doesn't matter
            self.assertEqual((r["protein"], r["dna"]), (p, d))
            self.assertEqual(r["other"], [os.path.join(REPO, "missing.fa")])

    def test_clustal_bom_crlf_and_unaligned_protein(self):
        tmp = tempfile.mkdtemp()
        try:
            def w(name, text, bom=False):
                path = os.path.join(tmp, name)
                with open(path, "wb") as f:
                    f.write((b"\xef\xbb\xbf" if bom else b"") + text.replace("\n", "\r\n").encode())
                return path
            self.assertEqual(cd.sniff_sequence_file(w("a.aln", "CLUSTAL W (1.83)\n\nsp1  MKV-LA\nsp2  MKVQLA\n       *** **\n")), "protein")
            self.assertEqual(cd.sniff_sequence_file(w("b.fa", ">x\nMKWVTFISLLFLFSSAYS\n", bom=True)), "protein")
            self.assertEqual(cd.sniff_sequence_file(w("c.fa", ">x\nATGAAAGTNCTGGCA\n>y\natgaaa\n", bom=True)), "dna")
            self.assertEqual(cd.sniff_sequence_file(w("d.txt", "")), "unknown")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class TestEngineDispatch(unittest.TestCase):
    """In a PyInstaller build the local app runs the engine as `COAL.exe -m coal ...`."""

    def test_dispatch_matches_python_m_coal(self):
        env = cd.child_env()
        env["PYTHONPATH"] = os.path.join(REPO, "src")
        tmp = tempfile.mkdtemp()
        try:
            for name, extra in (("clean", []), ("mito", ["-codontable", "2"]), ("one_broken_gene", [])):
                p = os.path.join(REPO, "examples", name, "protein.aln.fasta")
                d = os.path.join(REPO, "examples", name, "dna.fasta")
                a, b = os.path.join(tmp, name + ".a"), os.path.join(tmp, name + ".b")
                ra = subprocess.run([sys.executable, "-m", "coal", p, d, "-o", a] + extra,
                                    capture_output=True, env=env, timeout=120)
                rb = subprocess.run([sys.executable, os.path.join(HERE, "coal_desktop.py"), "-m", "coal",
                                     p, d, "-o", b] + extra, capture_output=True, env=env, timeout=120)
                self.assertEqual(ra.returncode, rb.returncode, name)
                self.assertEqual(ra.stderr, rb.stderr, name)
                if os.path.exists(a):
                    with open(a, "rb") as fa, open(b, "rb") as fb:
                        self.assertEqual(fa.read(), fb.read(), name)
                else:
                    self.assertFalse(os.path.exists(b), name)
            r = subprocess.run([sys.executable, os.path.join(HERE, "coal_desktop.py"), "-m", "coal", "--version"],
                               capture_output=True, text=True, env=env, timeout=60)
            self.assertEqual(r.returncode, 0)
            self.assertIn("0.", r.stdout + r.stderr)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_other_modules_refused(self):
        r = subprocess.run([sys.executable, os.path.join(HERE, "coal_desktop.py"), "-m", "os"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 2)
        self.assertIn("only '-m coal'", r.stderr)


# ============================================================== fake pywebview (documented 6.x API shape)

import threading  # noqa: E402
import types  # noqa: E402
from unittest import mock  # noqa: E402


class FakeEvent:
    def __init__(self):
        self.handlers = []

    def __iadd__(self, h):
        self.handlers.append(h)
        return self

    def fire(self):
        for h in list(self.handlers):
            h()


class FakeDocument:
    def __init__(self, url):
        self.url = url
        self.events = types.SimpleNamespace(dragover=FakeEvent(), drop=FakeEvent())


class FakeWindow:
    def __init__(self, title, **kw):
        self.title, self.kw = title, kw
        self.url = None                      # html= content: no http URL
        self.events = types.SimpleNamespace(closed=FakeEvent(), loaded=FakeEvent(), closing=FakeEvent())
        self.dom = types.SimpleNamespace(document=FakeDocument(None))
        self.documents = [self.dom.document]
        self.calls = []

    def get_current_url(self):
        return self.url

    def load_url(self, url):
        self.calls.append(("load_url", url))
        self.url = url
        self.dom.document = FakeDocument(url)   # a new page: handlers on the old one are gone
        self.documents.append(self.dom.document)
        self.events.loaded.fire()

    def destroy(self):
        self.calls.append(("destroy",))

    def restore(self):
        self.calls.append(("restore",))

    def show(self):
        self.calls.append(("show",))


def make_fake_webview(scenario):
    """scenario(window, func) plays what the user does while webview.start() runs."""
    wv = types.ModuleType("webview")
    dom = types.ModuleType("webview.dom")
    wv.settings = {"ALLOW_DOWNLOADS": False}
    wv.windows = []
    wv.start_kwargs = {}

    class DOMEventHandler:
        def __init__(self, callback, prevent_default=False, stop_propagation=False, debounce=0):
            self.callback, self.prevent_default = callback, prevent_default

    def create_window(title, **kw):
        if getattr(wv, "fail_create", None):
            raise wv.fail_create
        w = FakeWindow(title, **kw)
        wv.windows.append(w)
        return w

    def start(func=None, args=None, **kw):
        wv.start_kwargs = kw
        if getattr(wv, "fail_start", None):
            raise wv.fail_start
        w = wv.windows[0]
        w.events.loaded.fire()  # the splash page loaded
        scenario(w, lambda: func(args))

    dom.DOMEventHandler = DOMEventHandler
    wv.dom = dom
    wv.create_window, wv.start = create_window, start
    return wv


class WindowBase(Base):
    def setUp(self):
        super().setUp()
        self.dialogs, self.notices = [], []
        self._patches = [
            mock.patch.object(cd, "_error_dialog", side_effect=lambda t, x, offer_browser=False:
                              (self.dialogs.append((t, x)), False)[1]),
            mock.patch.object(cd, "_notice", side_effect=lambda t, x: self.notices.append((t, x)), create=True),
            mock.patch.object(cd, "data_dir", return_value=self.tmp),
        ]
        for p in self._patches:
            p.start()
        self.instance = cd.SingleInstance(os.path.join(self.tmp, "inst"))

    def tearDown(self):
        for p in self._patches:
            p.stop()
        super().tearDown()

    def run_window(self, server, scenario, files=(), **fail):
        wv = make_fake_webview(scenario)
        for k, v in fail.items():
            setattr(wv, k, v)
        with mock.patch.dict(sys.modules, {"webview": wv, "webview.dom": wv.dom}):
            code = cd.run_window(server, self.instance, list(files))
        return code, wv


def run_and_close(w, startup):
    t = threading.Thread(target=startup, daemon=True)  # pywebview runs func in a thread
    t.start()
    t.join(60)
    w.events.closed.fire()


class TestWindowDropGuard(WindowBase):
    """WIN-W1: the drop guard must be on the Gradio page, not only on the splash page."""

    def test_drop_handlers_are_on_the_app_page(self):
        s = self.server(start_timeout=30)
        code, wv = self.run_window(s, run_and_close)
        self.assertEqual(code, 0)
        w = wv.windows[0]
        app_docs = [d for d in w.documents if d.url and d.url.startswith("http://127.0.0.1:")]
        self.assertEqual(len(app_docs), 1)
        drop = app_docs[0].events.drop.handlers
        self.assertEqual(len(drop), 1, "no drop guard on the app page")
        self.assertTrue(drop[0].prevent_default)
        self.assertTrue(app_docs[0].events.dragover.handlers[0].prevent_default)
        self.assertEqual(w.documents[0].events.drop.handlers, [], "splash page should not need it")


class TestWindowClose(WindowBase):
    """WIN-W2 / WIN-R9: closing the window while COAL starts is not an error, and is quick."""

    def test_close_during_start_shows_no_error(self):
        starting = threading.Event()

        def command(port):
            starting.set()
            return [sys.executable, self.stub, str(port)]
        s = self.server(mode="hang", command=command, start_timeout=60, stop_timeout=2)
        threads, close_took = [], []

        def scenario(w, startup):
            t = threading.Thread(target=startup, daemon=True)
            threads.append(t)
            t.start()
            self.assertTrue(starting.wait(30))
            t0 = time.monotonic()
            w.events.closed.fire()  # the person closes the window during "Starting COAL..."
            close_took.append(time.monotonic() - t0)
        code, wv = self.run_window(s, scenario)
        threads[0].join(60)
        self.assertFalse(threads[0].is_alive())
        self.assertEqual(code, 0)
        self.assertEqual(self.dialogs, [], "error dialog shown after the person closed COAL")
        self.assertNotIn(("destroy",), wv.windows[0].calls)
        self.assertFalse(s.running())
        self.assertLess(close_took[0], 1.0, "closed handler blocked on stopping the server")


class TestWindowFallback(WindowBase):
    """WIN-R7: any failure creating/starting the window offers the browser fallback."""

    def test_create_window_failure_offers_browser(self):
        s = self.server(start_timeout=30)
        with mock.patch.object(cd, "_run_in_browser", return_value=7) as rb:
            cd._error_dialog.side_effect = lambda t, x, offer_browser=False: (self.dialogs.append((t, x)), True)[1]
            code, _ = self.run_window(s, run_and_close, fail_create=RuntimeError("WebView2 missing"))
        self.assertEqual(code, 7)
        rb.assert_called_once()
        self.assertIn("WebView2", self.dialogs[0][1])


class TestOpenWith(WindowBase):
    """WIN-W4: "Open with" / second launch with files shows a visible message."""

    def test_second_launch_with_file_explains_what_to_do(self):
        s = self.server(start_timeout=30)
        p = os.path.join(REPO, "examples", "clean", "protein.aln.fasta")

        def scenario(w, startup):
            t = threading.Thread(target=startup, daemon=True)
            t.start()
            t.join(60)
            self.instance.on_message({"action": "focus", "files": [p]})
            w.events.closed.fire()
        code, wv = self.run_window(s, scenario)
        self.assertEqual(code, 0)
        self.assertIn(("show",), wv.windows[0].calls)
        self.assertEqual(len(self.notices), 1, self.notices)
        text = self.notices[0][1]
        self.assertIn("protein.aln.fasta", text)
        self.assertIn("box 1", text)
        self.assertNotIn(REPO, text)

    def test_first_launch_with_file_explains_what_to_do(self):
        s = self.server(start_timeout=30)
        d = os.path.join(REPO, "examples", "clean", "dna.fasta")
        code, _ = self.run_window(s, run_and_close, files=[d])
        self.assertEqual(code, 0)
        self.assertEqual(len(self.notices), 1, self.notices)
        self.assertIn("box 2", self.notices[0][1])

    def test_browser_mode_answers_second_launch(self):
        s = self.server(start_timeout=30)
        opened = []
        with mock.patch("webbrowser.open", side_effect=opened.append):
            def dialog(t, x, offer_browser=False):
                if t == "COAL is running":  # "Press OK to quit": a second launch happens meanwhile
                    self.instance.on_message({"action": "focus", "files": []})
                return False
            cd._error_dialog.side_effect = dialog
            self.assertEqual(cd._run_in_browser(s, self.instance), 0)
        self.assertEqual(len(opened), 2, "second launch should reopen the page")


class TestLaunchLogging(unittest.TestCase):
    """WIN-W3: only the instance that owns the lock rotates logs; local_app.log is rotated too."""

    def setUp(self):
        import logging
        self.tmp = tempfile.mkdtemp(prefix="coal_logs_")
        self.logs = os.path.join(self.tmp, "logs")
        os.makedirs(self.logs)
        self.root_handlers = list(logging.getLogger().handlers)
        self.env = mock.patch.dict(os.environ, {"COAL_DESKTOP_HOME": self.tmp})
        self.env.start()

    def tearDown(self):
        import logging
        root = logging.getLogger()
        for h in list(root.handlers):
            if h not in self.root_handlers:
                root.removeHandler(h)
                h.close()
        self.env.stop()
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write(self, name, text):
        with open(os.path.join(self.logs, name), "w", encoding="utf-8") as f:
            f.write(text)

    def test_second_launch_leaves_the_running_instances_log_alone(self):
        got = []
        primary = cd.SingleInstance(self.tmp, on_message=got.append)
        self.assertTrue(primary.acquire())
        try:
            self.write("coal_desktop.log", "PRIMARY\n")
            self.assertEqual(cd.main_window(cd.parse_args([])), 0)
        finally:
            primary.release()
        with open(os.path.join(self.logs, "coal_desktop.log"), encoding="utf-8") as f:
            self.assertEqual(f.read(), "PRIMARY\n")
        self.assertFalse(os.path.exists(os.path.join(self.logs, "coal_desktop.log.1")))

    def test_primary_rotates_local_app_log(self):
        self.write("local_app.log", "OLD RUN\n")
        with mock.patch.object(cd, "_run_in_browser", return_value=0):
            self.assertEqual(cd.main_window(cd.parse_args(["--browser"])), 0)
        with open(os.path.join(self.logs, "local_app.log.1"), encoding="utf-8") as f:
            self.assertEqual(f.read(), "OLD RUN\n")
        self.assertFalse(os.path.exists(os.path.join(self.logs, "local_app.log")))


class TestClaimInstance(unittest.TestCase):
    """WIN-R5: a lock that is still being released, or a primary still writing instance.json,
    gets a few seconds before COAL says it is "already running"."""

    def setUp(self):
        self.dir = tempfile.mkdtemp(prefix="coal_claim_")

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def test_unreachable_holder_gives_up_after_the_wait(self):
        import fcntl
        if cd.IS_WINDOWS:
            self.skipTest("POSIX lock used to simulate a dead listener")
        with open(os.path.join(self.dir, "instance.lock"), "a+b") as fh:
            fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            t0 = time.monotonic()
            self.assertEqual(cd.claim_instance(cd.SingleInstance(self.dir), {"action": "focus"}, wait=1.0),
                             "unreachable")
            self.assertGreaterEqual(time.monotonic() - t0, 0.9)
            self.assertLess(time.monotonic() - t0, 10)

    def test_lock_released_during_the_wait_is_taken(self):
        import fcntl
        if cd.IS_WINDOWS:
            self.skipTest("POSIX lock used to simulate a slow release")
        fh = open(os.path.join(self.dir, "instance.lock"), "a+b")
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        threading.Timer(0.5, fh.close).start()  # released like a crashed holder, a bit later
        inst = cd.SingleInstance(self.dir)
        try:
            self.assertEqual(cd.claim_instance(inst, {"action": "focus"}, wait=5.0), "primary")
        finally:
            inst.release()
            fh.close()

    def test_running_primary_is_notified(self):
        got = []
        primary = cd.SingleInstance(self.dir, on_message=got.append)
        self.assertTrue(primary.acquire())
        try:
            self.assertEqual(cd.claim_instance(cd.SingleInstance(self.dir), {"action": "focus"}, wait=3.0),
                             "notified")
        finally:
            primary.release()


if __name__ == "__main__":
    unittest.main()
