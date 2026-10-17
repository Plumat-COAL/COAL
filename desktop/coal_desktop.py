"""COAL desktop wrapper for Windows: the existing local app (coal_app) in a native window.

One program, three roles, chosen by the command line:

    coal_desktop.py [FILE ...]            open the COAL window (FILEs: from "Open with")
    coal_desktop.py --serve --port N      run the local app on 127.0.0.1:N (child of the window)
    coal_desktop.py -m coal ARGS...       run the engine (the local app starts it this way,
                                          because in a PyInstaller build sys.executable is COAL.exe)

Process model (see DESIGN.md, section 5.2):

    window process  --stdin pipe-->  local app process (Gradio)  -->  engine runs (python -m coal)

The parts are kept separate on purpose:

* ``find_free_port``, ``LocalServer``, ``serve``, ``SingleInstance``, ``classify_dropped_files``
  and ``dispatch`` use only the standard library. They are tested on any OS by
  ``desktop/test_smoke.py``.
* ``run_window`` and the helpers below it import pywebview. They are written from pywebview's
  documented API (6.x) and need a Windows machine to test.

The wrapper never reads or changes sequences; results are exactly what ``python -m coal`` gives.
"""
from __future__ import annotations

import hmac
import json
import logging
import multiprocessing
import os
import queue
import re
import secrets
import socket
import subprocess
import sys
import threading
import time

APP_NAME = "COAL"
WINDOW_TITLE = "COAL codon alignment"
READY_PREFIX = "COAL-READY"
EXIT_PORT_BUSY = 75          # serve(): the port was taken; the window process tries another one
EXIT_NO_GUI = 3              # serve(): the local app (Gradio) could not be imported
SERVER_START_TIMEOUT = 60.0  # first start after install can be slow while antivirus scans files
SERVER_STOP_TIMEOUT = 5.0
PORT_ATTEMPTS = 3
APP_MUTEX = "COAL-desktop-running"  # Inno Setup's AppMutex: setup/uninstall ask to close COAL first

FROZEN = bool(getattr(sys, "frozen", False))
HERE = os.path.dirname(os.path.abspath(__file__))
IS_WINDOWS = os.name == "nt"

log = logging.getLogger("coal_desktop")


# ============================================================== paths and environment

def setup_source_paths() -> None:
    """In a source checkout, use this repo's src/ (not some other installed copy) for us and
    for every child process (the local app and the engine runs it starts)."""
    if FROZEN:
        return
    src = os.path.normpath(os.path.join(HERE, "..", "src"))
    if os.path.isdir(os.path.join(src, "coal_app")):
        if src not in sys.path:
            sys.path.insert(0, src)
        parts = [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p and p != src]
        os.environ["PYTHONPATH"] = os.pathsep.join([src] + parts)


def data_dir() -> str:
    """Per-user folder for the lock, logs and the WebView2 profile (never for results).

    Windows: %LOCALAPPDATA%\\COAL. Elsewhere: ~/.local/state/coal-desktop.
    COAL_DESKTOP_HOME overrides it (used by the tests)."""
    d = os.environ.get("COAL_DESKTOP_HOME")
    if not d:
        if IS_WINDOWS and os.environ.get("LOCALAPPDATA"):
            d = os.path.join(os.environ["LOCALAPPDATA"], APP_NAME)
        else:
            d = os.path.join(os.path.expanduser("~"), ".local", "state", "coal-desktop")
    os.makedirs(d, exist_ok=True)
    return d


def bundled_examples_dir() -> str | None:
    """The examples folder inside a PyInstaller build (the spec copies it to <bundle>/examples)."""
    base = getattr(sys, "_MEIPASS", None)
    if base:
        d = os.path.join(base, "examples")
        if os.path.isdir(d):
            return d
    return None


def child_env() -> dict:
    """Environment for the local app: offline, and never route 127.0.0.1 through a proxy."""
    env = dict(os.environ)
    env["GRADIO_ANALYTICS_ENABLED"] = "False"
    no_proxy = [p for p in env.get("NO_PROXY", env.get("no_proxy", "")).split(",") if p]
    for h in ("127.0.0.1", "localhost"):
        if h not in no_proxy:
            no_proxy.append(h)
    env["NO_PROXY"] = env["no_proxy"] = ",".join(no_proxy)
    env["PYTHONUNBUFFERED"] = "1"
    env.setdefault("PYTHONIOENCODING", "utf-8")
    return env


def setup_logging(name: str = "coal_desktop.log", keep: int = 5) -> str:
    """Log to <data_dir>/logs/<name>; keep the newest `keep` old logs. Never logs sequences."""
    d = os.path.join(data_dir(), "logs")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, name)
    try:  # rotate: name.4 <- name.3 ... name.1 <- name
        for i in range(keep - 1, 0, -1):
            a, b = f"{path}.{i}", f"{path}.{i + 1}"
            if os.path.exists(a):
                os.replace(a, b)
        if os.path.exists(path):
            os.replace(path, path + ".1")
    except OSError:
        pass
    handler = logging.FileHandler(path, encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(logging.INFO)
    return path


# ============================================================== ports

def find_free_port(host: str = "127.0.0.1") -> int:
    """A port that is free right now on `host` (the OS picks it). Another program could take it
    before the local app binds it; LocalServer then retries with a new port."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind((host, 0))
        return s.getsockname()[1]


_PORT_BUSY = re.compile(r"empty port|address already in use|only one usage of each socket address|"
                        r"errno 98|errno 48|winerror 10048", re.I)


def is_port_busy_error(e: BaseException) -> bool:
    import errno
    if isinstance(e, OSError) and e.errno in (errno.EADDRINUSE, 10048):
        return True
    return bool(_PORT_BUSY.search(str(e)))


# ============================================================== child role: serve the local app

def launch_gradio(port: int, mcp: bool = False):
    """Start the existing Gradio GUI on 127.0.0.1:port without blocking; return (url, close)."""
    os.environ["GRADIO_ANALYTICS_ENABLED"] = "False"
    from coal_app import core, ui  # ImportError -> EXIT_NO_GUI
    ex = bundled_examples_dir()
    if ex:
        # core computes EXAMPLES_DIR as <coal_app>/../../examples, which doesn't exist in a
        # bundle. Proposed proper fix (COAL_EXAMPLES_DIR) is in overnight/windows-app/LOG.md.
        core.EXAMPLES_DIR = ex
    demo = ui.launch(host="127.0.0.1", port=port, mcp=mcp, demo_mode=False, prevent_thread_lock=True)
    return f"http://127.0.0.1:{port}/", demo.close


def _wait_for_parent_windows(parent_pid: int, on_exit) -> None:
    """Second lifeline on Windows: wait (blocking, no polling) on the parent's process handle."""
    import ctypes
    from ctypes import wintypes
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    k32.OpenProcess.restype = wintypes.HANDLE
    k32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    k32.WaitForSingleObject.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    SYNCHRONIZE, INFINITE = 0x00100000, 0xFFFFFFFF
    h = k32.OpenProcess(SYNCHRONIZE, False, parent_pid)
    if not h:
        return
    k32.WaitForSingleObject(h, INFINITE)
    on_exit()


def serve(port: int, mcp: bool = False, launcher=None, parent_pid: int | None = None,
          stdin=None, stdout=None) -> int:
    """Run the local app until stdin closes (the window process exited or asked us to stop).

    Protocol with LocalServer: print one line ``COAL-READY <url>`` when listening; exit with
    EXIT_PORT_BUSY if the port was taken; exit 0 after a clean shutdown."""
    launcher = launcher or launch_gradio
    stdin = stdin if stdin is not None else (sys.stdin.buffer if sys.stdin else None)
    stdout = stdout if stdout is not None else sys.stdout
    if stdout is None:  # no pipe to report readiness on: nothing useful we can do
        return 4
    try:
        url, close = launcher(port, mcp)
    except ImportError as e:
        print(f"coal_desktop: the local app could not start (missing {e.name or e}). "
              f"In a source checkout: pip install -e '.[gui]' -r desktop/requirements.txt", file=stdout, flush=True)
        return EXIT_NO_GUI
    except OSError as e:
        if is_port_busy_error(e):
            print(f"coal_desktop: port {port} is busy: {e}", file=stdout, flush=True)
            return EXIT_PORT_BUSY
        raise
    print(f"{READY_PREFIX} {url}", file=stdout, flush=True)
    stop = threading.Event()
    if parent_pid and IS_WINDOWS:
        threading.Thread(target=_wait_for_parent_windows, args=(parent_pid, stop.set),
                         daemon=True).start()
    if stdin is not None:
        def lifeline():
            try:
                while stdin.read(65536):  # returns b"" only at EOF: the parent closed the pipe or died
                    pass
            except (OSError, ValueError):
                pass
            stop.set()
        threading.Thread(target=lifeline, daemon=True).start()
    try:
        stop.wait()
    except KeyboardInterrupt:
        pass
    finally:
        try:
            close()
        except Exception as e:  # noqa: BLE001 - shutting down anyway
            print(f"coal_desktop: error while closing the server: {e}", file=stdout, flush=True)
    return 0


# ============================================================== window role: manage the child

class ServerError(RuntimeError):
    """The local app did not start. str() is a plain-English message for a dialog."""


def default_server_command(port: int, mcp: bool = False) -> list[str]:
    base = [sys.executable] if FROZEN else [sys.executable, os.path.abspath(__file__)]
    cmd = base + ["--serve", "--port", str(port), "--parent-pid", str(os.getpid())]
    if mcp:
        cmd.append("--mcp")
    return cmd


class _WindowsJob:
    """Windows Job Object with KILL_ON_JOB_CLOSE: when the window process ends, Windows kills the
    local app and every engine run it started. Untested here (Linux). Failure is only logged."""

    def __init__(self):
        import ctypes
        from ctypes import wintypes
        self._ctypes = ctypes
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)

        class IO_COUNTERS(ctypes.Structure):
            _fields_ = [(n, ctypes.c_ulonglong) for n in (
                "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
                "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]

        class BASIC(ctypes.Structure):
            _fields_ = [("PerProcessUserTimeLimit", ctypes.c_longlong),
                        ("PerJobUserTimeLimit", ctypes.c_longlong),
                        ("LimitFlags", wintypes.DWORD),
                        ("MinimumWorkingSetSize", ctypes.c_size_t),
                        ("MaximumWorkingSetSize", ctypes.c_size_t),
                        ("ActiveProcessLimit", wintypes.DWORD),
                        ("Affinity", ctypes.c_size_t),
                        ("PriorityClass", wintypes.DWORD),
                        ("SchedulingClass", wintypes.DWORD)]

        class EXTENDED(ctypes.Structure):
            _fields_ = [("BasicLimitInformation", BASIC),
                        ("IoInfo", IO_COUNTERS),
                        ("ProcessMemoryLimit", ctypes.c_size_t),
                        ("JobMemoryLimit", ctypes.c_size_t),
                        ("PeakProcessMemoryUsed", ctypes.c_size_t),
                        ("PeakJobMemoryUsed", ctypes.c_size_t)]

        k32.CreateJobObjectW.restype = wintypes.HANDLE
        k32.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
        k32.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
        k32.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
        self._k32 = k32
        self.handle = k32.CreateJobObjectW(None, None)
        if not self.handle:
            raise OSError(ctypes.get_last_error(), "CreateJobObjectW failed")
        info = EXTENDED()
        info.BasicLimitInformation.LimitFlags = 0x2000  # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        JobObjectExtendedLimitInformation = 9
        if not k32.SetInformationJobObject(self.handle, JobObjectExtendedLimitInformation,
                                           ctypes.byref(info), ctypes.sizeof(info)):
            raise OSError(ctypes.get_last_error(), "SetInformationJobObject failed")

    def add(self, proc: subprocess.Popen) -> None:
        if not self._k32.AssignProcessToJobObject(self.handle, int(proc._handle)):  # noqa: SLF001
            raise OSError(self._ctypes.get_last_error(), "AssignProcessToJobObject failed")


class LocalServer:
    """Start the local app in a child process on a free 127.0.0.1 port, wait until it says it is
    ready (with a timeout), keep its output in a log, and stop it for certain.

    ``command`` is a function port -> argv (tests pass a stand-in server)."""

    def __init__(self, command=None, start_timeout: float = SERVER_START_TIMEOUT,
                 stop_timeout: float = SERVER_STOP_TIMEOUT, log_path: str | None = None,
                 attempts: int = PORT_ATTEMPTS, env: dict | None = None):
        self.command = command or default_server_command
        self.start_timeout = start_timeout
        self.stop_timeout = stop_timeout
        self.log_path = log_path
        self.attempts = attempts
        self.env = env if env is not None else child_env()
        self.proc: subprocess.Popen | None = None
        self.url: str | None = None
        self.port: int | None = None
        self.output: list[str] = []  # last lines of the child's output (for error messages)
        self._lines: queue.Queue = queue.Queue()
        self._log_fh = None
        self._job = None
        self._lock = threading.Lock()

    # --- start
    def start(self) -> str:
        last = ""
        for attempt in range(1, self.attempts + 1):
            port = find_free_port()
            code, last = self._start_once(port)
            if code is None:
                return self.url
            if code != EXIT_PORT_BUSY:
                break
            log.warning("port %s was taken (attempt %s of %s)", port, attempt, self.attempts)
        self.stop()
        raise ServerError(last)

    def _start_once(self, port: int):
        """Returns (None, "") when ready, else (exit code or -1, message)."""
        if self.proc is not None:  # previous attempt (already exited or killed)
            self._kill()
            self._close_pipes(self.proc)
        self._close_log()
        if self.log_path:
            self._log_fh = open(self.log_path, "a", encoding="utf-8", errors="replace")
        flags = 0
        if IS_WINDOWS:
            flags = getattr(subprocess, "CREATE_NO_WINDOW", 0x08000000)  # no console flashing up
        self._lines = queue.Queue()
        self.output = []
        self.proc = subprocess.Popen(
            self.command(port), stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, env=self.env, creationflags=flags, close_fds=True)
        if IS_WINDOWS:
            try:
                self._job = self._job or _WindowsJob()
                self._job.add(self.proc)
            except Exception as e:  # noqa: BLE001 - the stdin lifeline still works
                log.warning("could not attach the local app to a job object: %s", e)
        threading.Thread(target=self._pump, args=(self.proc, self._lines), daemon=True).start()
        deadline = time.monotonic() + self.start_timeout
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                msg = (f"COAL's local app did not start within {self.start_timeout:.0f} seconds."
                       + self._log_hint())
                self._kill()
                return -1, msg
            try:
                line = self._lines.get(timeout=remaining)  # blocks; the reader thread feeds it
            except queue.Empty:
                continue  # deadline reached: handled at the top of the loop
            if line is None:  # child's output closed: it exited
                code = self._wait_exit(5)
                tail = "\n".join(self.output[-8:])
                if code == EXIT_PORT_BUSY:
                    return code, f"Port {port} was busy."
                if code == EXIT_NO_GUI:
                    return code, ("COAL's window could not load its screens (the GUI package is "
                                  "missing from this installation)." + self._log_hint() + "\n\n" + tail)
                return code if code is not None else -1, (
                    f"COAL's local app stopped while starting (exit code {code})." + self._log_hint()
                    + ("\n\n" + tail if tail else ""))
            if line.startswith(READY_PREFIX):
                self.url = line.split(None, 1)[1].strip()
                self.port = port
                log.info("local app ready at %s (pid %s)", self.url, self.proc.pid)
                return None, ""

    def _pump(self, proc, lines):
        """Read the child's output for as long as it runs (so the pipe never fills and blocks it)."""
        try:
            for raw in iter(proc.stdout.readline, b""):
                text = raw.decode("utf-8", "replace").rstrip("\r\n")
                with self._lock:
                    self.output.append(text)
                    del self.output[:-200]
                    if self._log_fh:
                        try:
                            self._log_fh.write(text + "\n")
                            self._log_fh.flush()
                        except (OSError, ValueError):
                            pass
                lines.put(text)
        except (OSError, ValueError):
            pass
        finally:
            lines.put(None)

    def _log_hint(self) -> str:
        return f" Details are in {self.log_path}." if self.log_path else ""

    def _wait_exit(self, timeout):
        try:
            return self.proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            return None

    def _kill(self):
        if self.proc and self.proc.poll() is None:
            self.proc.kill()
            self._wait_exit(5)

    # --- stop
    def running(self) -> bool:
        return self.proc is not None and self.proc.poll() is None

    def stop(self) -> int | None:
        """Close the lifeline (stdin), wait; then terminate; then kill. Safe to call twice."""
        p = self.proc
        if p is None:
            self._close_log()
            return None
        if p.poll() is None:
            try:
                p.stdin.close()
            except OSError:
                pass
            if self._wait_exit(self.stop_timeout) is None:
                log.warning("local app did not stop by itself; terminating")
                p.terminate()
                if self._wait_exit(3) is None:
                    log.warning("local app did not terminate; killing")
                    p.kill()
                    self._wait_exit(5)
        self._close_pipes(p)
        self._close_log()
        log.info("local app stopped (exit code %s)", p.returncode)
        return p.returncode

    @staticmethod
    def _close_pipes(p):
        for s in (p.stdin, p.stdout):
            try:
                s.close()
            except (OSError, AttributeError, ValueError):
                pass

    def _close_log(self):
        with self._lock:
            if self._log_fh:
                try:
                    self._log_fh.close()
                except OSError:
                    pass
                self._log_fh = None

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, *exc):
        self.stop()


# ============================================================== single instance

class SingleInstance:
    """At most one COAL window per user.

    The first instance holds an OS file lock on <dir>/instance.lock for its lifetime (released by
    the OS if it crashes) and listens on a random 127.0.0.1 port for messages from later
    instances, which must present a random token from <dir>/instance.json. A later instance
    sends {"action": "focus", "files": [...]} and exits."""

    def __init__(self, directory: str, on_message=None):
        self.dir = directory
        self.lock_path = os.path.join(directory, "instance.lock")
        self.info_path = os.path.join(directory, "instance.json")
        self.on_message = on_message
        self._fh = None
        self._sock = None
        self._token = None
        self._thread = None

    # --- primary
    def acquire(self) -> bool:
        os.makedirs(self.dir, exist_ok=True)
        fh = open(self.lock_path, "a+b")
        try:
            if IS_WINDOWS:
                import msvcrt
                fh.seek(0)
                msvcrt.locking(fh.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            fh.close()
            return False
        self._fh = fh
        self._start_listener()
        return True

    def _start_listener(self):
        self._token = secrets.token_hex(16)
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind(("127.0.0.1", 0))
        s.listen(4)
        self._sock = s
        info = {"pid": os.getpid(), "port": s.getsockname()[1], "token": self._token}
        tmp = self.info_path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(info, f)
        os.replace(tmp, self.info_path)
        self._thread = threading.Thread(target=self._accept_loop, daemon=True)
        self._thread.start()

    def _accept_loop(self):
        while True:
            try:
                conn, _ = self._sock.accept()
            except OSError:
                return  # listener closed by release()
            with conn:
                try:
                    conn.settimeout(2.0)
                    data = b""
                    while not data.endswith(b"\n") and len(data) < 65536:
                        chunk = conn.recv(4096)
                        if not chunk:
                            break
                        data += chunk
                    msg = json.loads(data.decode("utf-8"))
                    if not hmac.compare_digest(str(msg.get("token", "")), self._token):
                        conn.sendall(b"denied\n")
                        continue
                    msg.pop("token", None)
                    conn.sendall(b"ok\n")
                except (OSError, ValueError, AttributeError):
                    continue
            if self.on_message:
                try:
                    self.on_message(msg)
                except Exception as e:  # noqa: BLE001 - a bad message must not kill the listener
                    log.warning("instance message handler failed: %s", e)

    def release(self):
        if self._sock:
            try:
                self._sock.close()
            except OSError:
                pass
            self._sock = None
        if self._fh:
            try:
                os.remove(self.info_path)
            except OSError:
                pass
            try:
                if IS_WINDOWS:
                    import msvcrt
                    self._fh.seek(0)
                    msvcrt.locking(self._fh.fileno(), msvcrt.LK_UNLCK, 1)
            except OSError:
                pass
            self._fh.close()
            self._fh = None

    # --- secondary
    def notify_primary(self, message: dict, timeout: float = 5.0) -> bool:
        """Send a message to the running instance. False if it could not be reached."""
        try:
            with open(self.info_path, encoding="utf-8") as f:
                info = json.load(f)
            payload = dict(message, token=info["token"])
            with socket.create_connection(("127.0.0.1", int(info["port"])), timeout=timeout) as c:
                c.settimeout(timeout)
                c.sendall(json.dumps(payload).encode("utf-8") + b"\n")
                return c.recv(16).startswith(b"ok")
        except (OSError, ValueError, KeyError):
            return False


# ============================================================== dropped files

_DNA = set("ACGTUNRYSWKMBDHV")


def sniff_sequence_file(path: str, max_bytes: int = 256 * 1024) -> str:
    """'protein', 'dna' or 'unknown', from the first `max_bytes` of a FASTA/CLUSTAL file.
    Aligned protein: has gaps or amino-acid letters. DNA: only nucleotide (IUPAC) letters."""
    try:
        with open(path, "rb") as f:
            head = f.read(max_bytes).decode("utf-8-sig", "replace")
    except OSError:
        return "unknown"
    lines = head.splitlines()
    if len(head) == max_bytes and lines:
        lines = lines[:-1]  # last line may be cut
    letters, gaps = [], 0
    clustal = bool(lines) and lines[0].upper().startswith(("CLUSTAL", "MUSCLE"))
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith((">", ";", "#")) or (clustal and (s.upper().startswith(("CLUSTAL", "MUSCLE"))
                                                                  or set(s) <= set("*:. "))):
            continue
        if clustal:
            parts = s.split()
            if len(parts) < 2:
                continue
            s = parts[1]
        gaps += s.count("-") + s.count(".")
        letters.extend(c for c in s.upper() if c.isalpha() or c == "*")
    if not letters:
        return "unknown"
    nonnuc = sum(1 for c in letters if c not in _DNA)
    if nonnuc == 0 and gaps == 0:
        return "dna"
    if clustal or gaps or nonnuc / len(letters) > 0.02:
        return "protein"
    return "dna"


def classify_dropped_files(paths) -> dict:
    """Sort dropped paths into {'protein': path|None, 'dna': path|None, 'other': [paths]}."""
    out = {"protein": None, "dna": None, "other": []}
    for p in paths:
        if not p or not os.path.isfile(p):
            out["other"].append(p)
            continue
        kind = sniff_sequence_file(p)
        if kind in ("protein", "dna") and out[kind] is None:
            out[kind] = p
        else:
            out["other"].append(p)
    return out


# ============================================================== engine role

def run_engine(args: list[str]) -> int:
    """`COAL.exe -m coal ARGS` -> the engine, exactly as `python -m coal ARGS`."""
    import runpy
    sys.argv = ["coal"] + list(args)
    try:
        runpy.run_module("coal", run_name="__main__", alter_sys=True)
    except SystemExit as e:
        code = e.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        print(code, file=sys.stderr)
        return 1
    return 0


# ============================================================== window (pywebview; untested here)

SPLASH = """<!doctype html><html><head><meta charset="utf-8"><title>COAL</title>
<style>
:root { color-scheme: light dark; }
body { font: 18px system-ui, "Segoe UI", sans-serif; display:flex; align-items:center;
       justify-content:center; height:100vh; margin:0; background:Canvas; color:CanvasText; }
main { text-align:center; } p { opacity:.8; }
</style></head><body><main role="status" aria-live="polite">
<h1>Starting COAL&hellip;</h1><p>This computer only. No internet needed.</p></main></body></html>"""


def _error_dialog(title: str, text: str, offer_browser: bool = False) -> bool:
    """Native message box (Windows) or stderr. Returns True if the person chose 'Yes'."""
    if IS_WINDOWS:
        try:
            import ctypes
            MB_ICONWARNING, MB_YESNO, MB_OK, IDYES = 0x30, 0x4, 0x0, 6
            flags = MB_ICONWARNING | (MB_YESNO if offer_browser else MB_OK)
            return ctypes.windll.user32.MessageBoxW(None, text, title, flags) == IDYES
        except Exception:  # noqa: BLE001
            pass
    print(f"{title}: {text}", file=sys.stderr)
    return False


def run_window(server: LocalServer, instance: SingleInstance, files: list[str]) -> int:
    """Show the native window. UNTESTED (needs Windows + pywebview); written from pywebview 6.x docs:
    create_window, window.events.closed, window.dom.document.events.drop with DOMEventHandler
    (full paths in dataTransfer.files[i].pywebviewFullPath), load_url, restore/show, on_top,
    webview.settings['ALLOW_DOWNLOADS'], webview.start(func, args, gui, private_mode, storage_path)."""
    import webview
    from webview.dom import DOMEventHandler

    # Gradio's download buttons need this; pywebview disables downloads by default.
    webview.settings["ALLOW_DOWNLOADS"] = True

    window = webview.create_window(WINDOW_TITLE, html=SPLASH, width=1200, height=820,
                                   min_size=(900, 600), text_select=True, zoomable=True,
                                   background_color="#FFFFFF")
    failure = {}

    def bring_to_front(msg: dict):
        log.info("message from a second launch: %s", {k: v for k, v in msg.items() if k != "files"})
        if msg.get("files"):
            log.info("files given to the second launch: %s", classify_dropped_files(msg["files"]))
        try:
            window.restore()
            window.show()
            window.on_top = True   # brief raise above other windows, then back to normal
            window.on_top = False
        except Exception as e:  # noqa: BLE001
            log.warning("could not bring the window to the front: %s", e)

    instance.on_message = bring_to_front

    def on_drag(_e):
        pass

    def on_drop(e):
        # Drops on Gradio's own drop zones are handled by Gradio (it uploads over 127.0.0.1).
        # This document-level handler mainly stops the window from navigating away to a file
        # dropped elsewhere, and records the full paths for the planned path-based input.
        try:
            paths = [f.get("pywebviewFullPath") for f in e.get("dataTransfer", {}).get("files", [])]
            log.info("files dropped on the window: %s", classify_dropped_files([p for p in paths if p]))
        except Exception as ex:  # noqa: BLE001
            log.warning("could not read the dropped files: %s", ex)

    def on_closed():
        server.stop()

    window.events.closed += on_closed

    def startup(win):
        try:
            win.dom.document.events.dragover += DOMEventHandler(on_drag, True, False, debounce=500)
            win.dom.document.events.drop += DOMEventHandler(on_drop, True, False)
        except Exception as ex:  # noqa: BLE001 - drag-and-drop is optional
            log.warning("drag and drop not available: %s", ex)
        try:
            url = server.start()
        except ServerError as ex:
            failure["error"] = str(ex)
            win.destroy()
            return
        if files:
            log.info("opened with files: %s", classify_dropped_files(files))
        win.load_url(url)

    gui = "edgechromium" if IS_WINDOWS else None  # never fall back to the IE11 engine
    try:
        webview.start(startup, window, gui=gui, private_mode=True,
                      storage_path=os.path.join(data_dir(), "webview"), debug=False)
    except Exception as ex:  # noqa: BLE001 - typically: WebView2 runtime missing
        log.exception("the window could not open")
        if _error_dialog("COAL",
                         "COAL could not open its window. This usually means the Microsoft Edge "
                         "WebView2 Runtime is missing (it is part of Windows 11 and of Edge updates "
                         "on Windows 10).\n\nOpen COAL in your web browser instead? It still runs "
                         f"only on this computer.\n\n(Details: {ex})", offer_browser=True):
            return _run_in_browser(server)
        return 1
    finally:
        server.stop()
    if failure:
        _error_dialog("COAL could not start", failure["error"])
        return 1
    return 0


def _run_in_browser(server: LocalServer) -> int:
    """Fallback when no native window is possible: default browser + a small 'Quit COAL' dialog."""
    import webbrowser
    try:
        url = server.url or server.start()
    except ServerError as ex:
        _error_dialog("COAL could not start", str(ex))
        return 1
    webbrowser.open(url)
    _error_dialog("COAL is running", f"COAL is open in your browser at {url}\n\n"
                  "Press OK to quit COAL.")
    server.stop()
    return 0


# ============================================================== entry point

def parse_args(argv):
    import argparse
    ap = argparse.ArgumentParser(prog="COAL", description="COAL desktop app")
    ap.add_argument("files", nargs="*", help="sequence files to open (from 'Open with')")
    ap.add_argument("--serve", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--port", type=int, default=0, help=argparse.SUPPRESS)
    ap.add_argument("--parent-pid", type=int, default=None, help=argparse.SUPPRESS)
    ap.add_argument("--mcp", action="store_true", help="also serve the MCP tools to AI agents on this computer")
    ap.add_argument("--browser", action="store_true", help="open in the default browser instead of a window")
    ap.add_argument("--version", action="store_true")
    return ap.parse_args(argv)


def dispatch(argv: list[str]) -> int:
    if argv[:2] == ["-m", "coal"]:
        return run_engine(argv[2:])
    if argv[:1] == ["-m"]:
        print(f"COAL: only '-m coal' is supported here, not '-m {' '.join(argv[1:2])}'", file=sys.stderr)
        return 2
    args = parse_args(argv)
    if args.version:
        from coal import __version__ as v  # noqa: PLC0415
        print(f"COAL desktop {v}")
        return 0
    if args.serve:
        if not args.port:
            print("COAL: --serve needs --port", file=sys.stderr)
            return 2
        return serve(args.port, mcp=args.mcp, parent_pid=args.parent_pid)
    return main_window(args)


def _hold_app_mutex():
    """Windows: create a named mutex for as long as COAL runs (the handle is released by Windows
    when the process ends). Only used so the installer/uninstaller can tell COAL is open.
    Untested here (Linux)."""
    if not IS_WINDOWS:
        return None
    try:
        import ctypes
        return ctypes.windll.kernel32.CreateMutexW(None, False, APP_MUTEX)
    except Exception as e:  # noqa: BLE001
        log.warning("could not create the app mutex: %s", e)
        return None


def main_window(args) -> int:
    log_path = setup_logging()
    _mutex = _hold_app_mutex()  # noqa: F841 - kept alive until main_window returns
    log.info("COAL desktop starting (frozen=%s, python %s)", FROZEN, sys.version.split()[0])
    instance = SingleInstance(data_dir())
    if not instance.acquire():
        files = [os.path.abspath(f) for f in args.files]
        if instance.notify_primary({"action": "focus", "files": files}):
            log.info("COAL is already open; asked it to come to the front")
            return 0
        _error_dialog("COAL", "COAL seems to be already running but did not respond. "
                      "Close it from the taskbar (or Task Manager) and try again.")
        return 1
    server = LocalServer(command=lambda port: default_server_command(port, mcp=args.mcp),
                         log_path=os.path.join(os.path.dirname(log_path), "local_app.log"))
    try:
        if args.browser:
            return _run_in_browser(server)
        try:
            import webview  # noqa: F401
        except ImportError:
            _error_dialog("COAL", "The window component (pywebview) is not installed. "
                          "Install desktop/requirements.txt, or start with --browser.")
            return 1
        return run_window(server, instance, [os.path.abspath(f) for f in args.files])
    finally:
        server.stop()
        instance.release()


def main(argv=None) -> int:
    multiprocessing.freeze_support()  # the engine's worker processes in a PyInstaller build
    setup_source_paths()
    return dispatch(list(sys.argv[1:] if argv is None else argv))


if __name__ == "__main__":
    sys.exit(main())
