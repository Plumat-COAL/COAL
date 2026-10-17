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
            again = cd.SingleInstance(self.dir)
            self.assertTrue(again.acquire())  # the OS released the dead process's lock
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


if __name__ == "__main__":
    unittest.main()
