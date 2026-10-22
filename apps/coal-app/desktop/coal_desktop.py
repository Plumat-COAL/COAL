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
from pathlib import Path
import queue
import re
import secrets
import socket
import subprocess
import sys
import threading
import time
import traceback

APP_NAME = "COAL"
WINDOW_TITLE = "COAL codon alignment"
READY_PREFIX = "COAL-READY"
EXIT_PORT_BUSY = 75          # serve(): the port was taken; the window process tries another one
EXIT_NO_GUI = 3              # serve(): the local app (Gradio) could not be imported
EXIT_CRASH = 1               # serve()/engine: unexpected error (traceback printed to the pipe)
# 180 s: the first start after install can be slow while Defender scans the bundle.
SERVER_START_TIMEOUT = 180.0
SERVER_STOP_TIMEOUT = 5.0
PORT_ATTEMPTS = 3
APP_MUTEX = "COAL-desktop-running"  # Inno Setup's AppMutex: setup/uninstall ask to close COAL first

FROZEN = bool(getattr(sys, "frozen", False))
HERE = os.path.dirname(os.path.abspath(__file__))
IS_WINDOWS = os.name == "nt"

log = logging.getLogger("coal_desktop")


# ============================================================== paths and environment

def source_dirs() -> list[str]:
    """In a source checkout (apps/coal-app/desktop/): the app's src/ (apps/coal-app/src, coal_app) and the engine's
    src/ (the repository root's src, coal and coal_compat), each only if it is there."""
    app = os.path.normpath(os.path.join(HERE, "..", "src"))
    engine = os.path.normpath(os.path.join(HERE, "..", "..", "..", "src"))
    return [d for d, pkg in ((app, "coal_app"), (engine, "coal")) if os.path.isdir(os.path.join(d, pkg))]


def setup_source_paths() -> None:
    """In a source checkout, use this repo's app and engine sources (not some other installed copy) for us and
    for every child process (the local app and the engine runs it starts)."""
    if FROZEN:
        return
    srcs = source_dirs()
    if srcs:
        for src in reversed(srcs):
            if src in sys.path:
                sys.path.remove(src)
            sys.path.insert(0, src)
        parts = [p for p in os.environ.get("PYTHONPATH", "").split(os.pathsep) if p and p not in srcs]
        os.environ["PYTHONPATH"] = os.pathsep.join(srcs + parts)


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
    env["NO_PROXY"] = ",".join(no_proxy)
    # no_proxy only off Windows: Windows variable names are case-insensitive, so two spellings in one environment
    # block are duplicates.
    if not IS_WINDOWS:
        env["no_proxy"] = env["NO_PROXY"]
    env["PYTHONUNBUFFERED"] = "1"
    env.setdefault("PYTHONIOENCODING", "utf-8")
    return env


def rotate_log(path: str, keep: int) -> None:
    """path.<keep> <- ... <- path.1 <- path. Only the instance that holds the lock calls this."""
    try:
        for i in range(keep - 1, 0, -1):
            a, b = f"{path}.{i}", f"{path}.{i + 1}"
            if os.path.exists(a):
                os.replace(a, b)
        if os.path.exists(path):
            os.replace(path, path + ".1")
    except OSError as e:
        log.warning("could not rotate %s: %s", path, e)


def setup_logging(name: str = "coal_desktop.log", keep: int = 5) -> str:
    """Log to <data_dir>/logs/<name>; keep the newest `keep` old logs. Never logs sequences."""
    d = os.path.join(data_dir(), "logs")
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, name)
    rotate_log(path, keep)
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


# "Port busy": EADDRINUSE / WinError 10048, and also WinError 10013 (WSAEACCES: Windows' answer for a port held
# with exclusive use or inside a reserved/excluded port range) and EACCES raised by a bind.
_PORT_BUSY = re.compile(r"empty port|address already in use|only one usage of each socket address|"
                        r"errno 98|errno 48|winerror 10048|winerror 10013|"
                        r"access a socket in a way forbidden by its access permissions", re.I)
_WSAEACCES, _WSAEADDRINUSE = 10013, 10048


def _raised_by_bind(e: BaseException) -> bool:
    """True if a function with 'bind' in its name (socket bind, socketserver's server_bind, uvicorn's
    bind_socket, ...) is on the traceback: an EACCES there is about the port, not about a file."""
    tb = e.__traceback__
    while tb is not None:
        if "bind" in tb.tb_frame.f_code.co_name.lower():
            return True
        tb = tb.tb_next
    return False


def is_port_busy_error(e: BaseException) -> bool:
    import errno
    if not isinstance(e, OSError):
        return bool(_PORT_BUSY.search(str(e)))
    if getattr(e, "winerror", None) in (_WSAEACCES, _WSAEADDRINUSE):
        return True
    if e.errno in (errno.EADDRINUSE, _WSAEADDRINUSE):
        return True
    if e.errno == errno.EACCES and _raised_by_bind(e):
        return True
    return bool(_PORT_BUSY.search(str(e)))


# ============================================================== child role: serve the local app

def launch_gradio(port: int, mcp: bool = False):
    """Start the existing Gradio GUI on 127.0.0.1:port without blocking; return (url, close)."""
    os.environ["GRADIO_ANALYTICS_ENABLED"] = "False"
    from coal_app import core, ui  # ImportError -> EXIT_NO_GUI
    ex = bundled_examples_dir()
    if ex:
        # core computes EXAMPLES_DIR as <coal_app>/../../examples, which doesn't exist in a bundle.
        # core reads COAL_EXAMPLES_DIR, so it is set here too.
        os.environ["COAL_EXAMPLES_DIR"] = ex
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
        # Installed app -> "reinstall"; source checkout -> pip with cmd-safe double quotes.
        if FROZEN:
            hint = "This installation of COAL is incomplete; please reinstall COAL."
        else:
            hint = ('In a source checkout: pip install -e . -e "apps/coal-app[gui]" '
                    '-r apps/coal-app/desktop/requirements.txt')
        print(f"coal_desktop: the local app could not start (missing {e.name or e}). {hint}",
              file=stdout, flush=True)
        return EXIT_NO_GUI
    except OSError as e:
        if is_port_busy_error(e):
            print(f"coal_desktop: port {port} is busy: {e}", file=stdout, flush=True)
            return EXIT_PORT_BUSY
        traceback.print_exc(file=stdout)
        return EXIT_CRASH
    except Exception:  # noqa: BLE001
        # Printed to the pipe (-> local_app.log), exit 1, so the windowed build never shows a modal
        # "Unhandled exception" box in a hidden process.
        traceback.print_exc(file=stdout)
        stdout.flush()
        return EXIT_CRASH
    print(f"{READY_PREFIX} {url}", file=stdout, flush=True)
    if mcp:
        print(f"coal_desktop: MCP tools for AI agents at {url}gradio_api/mcp/sse", file=stdout, flush=True)
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
    """The local app did not start. str() is a plain-English message for a dialog.
    ``cancelled``: start() gave up because stop() was called (the person closed COAL)."""

    def __init__(self, message: str = "", cancelled: bool = False):
        super().__init__(message)
        self.cancelled = cancelled


def default_server_command(port: int, mcp: bool = False) -> list[str]:
    base = [sys.executable] if FROZEN else [sys.executable, os.path.abspath(__file__)]
    cmd = base + ["--serve", "--port", str(port), "--parent-pid", str(os.getpid())]
    # The window's default is MCP off (mcp=False): the endpoint sits on a random port, so no MCP client could be registered with it. `COAL.exe --mcp` turns it on; the stdio server (desktop/setup_mcp.bat)
    # is the way to register COAL with Claude. The --serve child still serves MCP unless given --no-mcp.
    if not mcp:
        cmd.append("--no-mcp")
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
        self._stop_lock = threading.Lock()
        self._pump_thread = None
        self._stopped = False
        self._job_failed = False

    # --- start
    def start(self) -> str:
        last = ""
        for attempt in range(1, self.attempts + 1):
            # stop() during start(): raise ServerError(cancelled=True), so it doesn't look like a crash or start a
            # new attempt.
            if self._stopped:
                raise ServerError("COAL was closed while it was starting.", cancelled=True)
            port = find_free_port()
            code, last = self._start_once(port)
            if code is None:
                return self.url
            if self._stopped:
                self.stop()
                raise ServerError("COAL was closed while it was starting.", cancelled=True)
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
        argv = self.command(port)
        with self._stop_lock:  # WIN-W2: stop() runs either before this (no child) or after (sees it)
            if self._stopped:
                return -1, "COAL was closed while it was starting."
            self.proc = subprocess.Popen(
                argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, env=self.env, creationflags=flags, close_fds=True)
        if IS_WINDOWS:
            try:
                self._job = self._job or _WindowsJob()
                self._job.add(self.proc)
                self._job_failed = False
            except Exception as e:  # noqa: BLE001 - the stdin lifeline still works
                # stop() then kills the process tree (taskkill /T /F) so engine runs don't outlive COAL.
                log.warning("could not attach the local app to a job object (%s); "
                            "will stop it with taskkill /T", e)
                self._job_failed = True
        self._pump_thread = threading.Thread(target=self._pump, args=(self.proc, self._lines), daemon=True)
        self._pump_thread.start()
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

    def _kill_tree(self):
        """Windows without a job object: kill the local app AND the engine runs it started.
        Untested here (Linux)."""
        if not (IS_WINDOWS and self._job_failed and self.proc and self.proc.poll() is None):
            return
        try:
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(self.proc.pid)], capture_output=True,
                           timeout=15, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0x08000000))
        except (OSError, subprocess.TimeoutExpired) as e:
            log.warning("taskkill failed: %s", e)

    def _kill(self):
        if self.proc and self.proc.poll() is None:
            self._kill_tree()
            if self.proc.poll() is None:
                self.proc.kill()
            self._wait_exit(5)

    # --- stop
    def running(self) -> bool:
        return self.proc is not None and self.proc.poll() is None

    def stop(self) -> int | None:
        """Close the lifeline (stdin), wait; then terminate; then kill. Safe to call twice, and
        from another thread while start() is waiting (start() then raises, cancelled)."""
        self._stopped = True
        with self._stop_lock:
            return self._stop_locked()

    def _stop_locked(self) -> int | None:
        p = self.proc
        if p is None:
            self._close_log()
            return None
        if p.poll() is None and IS_WINDOWS and self._job_failed:
            self._kill_tree()  # WIN-R10: no graceful close, but no orphaned engine runs either
            self._wait_exit(5)
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

    def _close_pipes(self, p):
        # Let the reader thread see EOF first: closing stdout under a blocked readline() can block.
        t = self._pump_thread
        if t is not None and t is not threading.current_thread() and p.poll() is not None:
            t.join(5)
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

def lock_file(fh) -> None:
    """Non-blocking exclusive OS lock on an open file (byte 0 on Windows); OSError if another handle holds it.
    Separate from SingleInstance.acquire so tests lock the same way on every OS."""
    if IS_WINDOWS:
        import msvcrt
        fh.seek(0)
        msvcrt.locking(fh.fileno(), msvcrt.LK_NBLCK, 1)
    else:
        import fcntl
        fcntl.flock(fh.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)


def unlock_file(fh) -> None:
    """Release lock_file(fh). (Closing the file also releases it, asynchronously on Windows.)"""
    if IS_WINDOWS:
        import msvcrt
        fh.seek(0)
        msvcrt.locking(fh.fileno(), msvcrt.LK_UNLCK, 1)
    else:
        import fcntl
        fcntl.flock(fh.fileno(), fcntl.LOCK_UN)


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
            lock_file(fh)
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
        sock = self._sock  # release() sets self._sock to None after closing it
        while True:
            try:
                conn, _ = sock.accept()
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
                if IS_WINDOWS:  # POSIX: close() below releases the flock at once
                    unlock_file(self._fh)
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
            if IS_WINDOWS:
                # Besides toggling on_top, let the running COAL take the foreground (Windows' foreground
                # lock). Untested here.
                try:
                    import ctypes
                    ctypes.windll.user32.AllowSetForegroundWindow(int(info["pid"]))
                except Exception:  # noqa: BLE001
                    pass
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
    """Sort dropped paths into {'protein': path|None, 'dna': path|None, 'other': [paths], 'extra': [(kind, path)]}
    ('extra': a second protein or DNA file; only the first of each kind is used)."""
    # A second file of a kind goes in 'extra', not 'other', so files_to_load can say only one is used.
    out = {"protein": None, "dna": None, "other": [], "extra": []}
    for p in paths:
        if not p or not os.path.isfile(p):
            out["other"].append(p)
            continue
        kind = sniff_sequence_file(p)
        if kind in ("protein", "dna") and out[kind] is None:
            out[kind] = p
        elif kind in ("protein", "dna"):
            out["extra"].append((kind, p))
        else:
            out["other"].append(p)
    return out


# ============================================================== engine role

def run_engine(args: list[str], module: str = "coal") -> int:
    """`COAL.exe -m coal ARGS` -> the engine, exactly as `python -m coal ARGS`."""
    import runpy
    sys.argv = ["coal"] + list(args)
    try:
        runpy.run_module(module, run_name="__main__", alter_sys=True)
    except Exception:  # noqa: BLE001
        # Traceback on stderr (the GUI shows it), exit 1, like python -m coal; in the windowed build a
        # propagated error would be a modal traceback box blocking the run until its timeout.
        if sys.stderr is not None:
            traceback.print_exc(file=sys.stderr)
        return 1
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

# No reassurance line under the title ("This computer only. No internet needed."): the privacy facts are in
# desktop/README.md and on the site.
SPLASH = """<!doctype html><html><head><meta charset="utf-8"><title>COAL</title>
<style>
:root { color-scheme: light dark; }
body { font: 18px system-ui, "Segoe UI", sans-serif; display:flex; align-items:center;
       justify-content:center; height:100vh; margin:0; background:Canvas; color:CanvasText; }
main { text-align:center; } p { opacity:.8; }
</style></head><body><main role="status" aria-live="polite">
<h1>Starting COAL&hellip;</h1>
<p id="slow" hidden>First launch can take a minute while Windows checks the app.</p></main>
<script>setTimeout(function () { document.getElementById("slow").hidden = false; }, 10000);</script>
</body></html>"""


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


def _notice(title: str, text: str) -> None:
    """Information box that does not block the caller (Windows), or stderr."""
    if IS_WINDOWS:
        def show():
            try:
                import ctypes
                MB_ICONINFORMATION, MB_SETFOREGROUND, MB_TOPMOST = 0x40, 0x10000, 0x40000
                ctypes.windll.user32.MessageBoxW(None, text, title,
                                                 MB_ICONINFORMATION | MB_SETFOREGROUND | MB_TOPMOST)
            except Exception:  # noqa: BLE001
                pass
        threading.Thread(target=show, daemon=True).start()
        return
    print(f"{title}: {text}", file=sys.stderr)


def _names(classified: dict) -> dict:
    """classify_dropped_files() result with file names only (logs never hold full paths)."""
    b = lambda p: os.path.basename(p) if p else p  # noqa: E731
    return {"protein": b(classified["protein"]), "dna": b(classified["dna"]),
            "other": [b(p) for p in classified["other"]]}


def query_url(base: str, protein: str | None, dna: str | None) -> str:
    """The local app's URL with its ?protein=&dna= path hook filled in (URL-encoded); see core.paths_from_query."""
    from urllib.parse import urlencode
    q = {k: v for k, v in (("protein", protein), ("dna", dna)) if v}
    return base.rstrip("/") + "/" + ("?" + urlencode(q) if q else "")


def files_to_load(paths, current: dict) -> tuple[dict, str]:
    """Sort files from "Open with", a second launch or a drop into protein / DNA
    (by sniffing their content), keep the other box's file from `current` if only one kind is given, and
    return ({"protein": path|None, "dna": path|None}, message for the person or "")."""
    c = classify_dropped_files(paths)
    want = {"protein": c["protein"] or current.get("protein"), "dna": c["dna"] or current.get("dna")}
    lines = []
    for p in c["other"]:
        lines.append(f"{os.path.basename(p or '')}: COAL couldn't tell whether this is a protein alignment or "
                     "DNA. Type or paste its path into the box it belongs in.")
    # A second file of the same kind.
    for kind in ("protein", "dna"):
        if any(k == kind for k, _ in c["extra"]):
            what = "protein alignment" if kind == "protein" else "DNA file"
            lines.append(f"Only one {what} can be loaded at a time; using {c[kind]}.")
    if (c["protein"] or c["dna"]) and not (want["protein"] and want["dna"]):
        lines.append("Now add the coding DNA (box 2)." if want["protein"] else
                     "Now add the protein alignment (box 1).")
    return want, "\n".join(lines)


def open_with_message(paths) -> str:
    """Files are loaded by path, so this only explains what could not be loaded (kept for callers/tests)."""
    return files_to_load(paths, {})[1]


def run_window(server: LocalServer, instance: SingleInstance, files: list[str]) -> int:
    """Show the native window. UNTESTED on Windows (tested here only against a fake pywebview in
    test_wrapper.py); written from pywebview 6.x docs: create_window, window.events.loaded/closed,
    window.get_current_url, window.dom.document.events.drop with DOMEventHandler (full paths in
    dataTransfer.files[i].pywebviewFullPath), load_url, restore/show, on_top,
    webview.settings['ALLOW_DOWNLOADS'], webview.start(func, args, gui, private_mode, storage_path)."""
    import webview
    from webview.dom import DOMEventHandler

    # Gradio's download buttons need this; pywebview disables downloads by default.
    webview.settings["ALLOW_DOWNLOADS"] = True

    state = {"closing": False, "failure": None, "stopper": None, "loaded": {}}
    window = None

    def load_files(paths):
        # The files are loaded by path (read in place, never copied) through the GUI's ?protein=&dna= hook.
        want, msg = files_to_load(paths, state["loaded"])
        log.info("loading files: %s", _names(classify_dropped_files([p for p in want.values() if p])))
        if (want["protein"] or want["dna"]) and server.url:
            state["loaded"] = want
            window.load_url(query_url(server.url, want["protein"], want["dna"]))
        if msg:
            _notice("COAL", msg)

    def bring_to_front(msg: dict):
        log.info("message from a second launch: %s", {k: v for k, v in msg.items() if k != "files"})
        try:
            window.restore()
            window.show()
            window.on_top = True   # brief raise above other windows, then back to normal
            window.on_top = False
        except Exception as e:  # noqa: BLE001
            log.warning("could not bring the window to the front: %s", e)
        if msg.get("files"):
            load_files(msg["files"])

    instance.on_message = bring_to_front

    def on_drag(_e):
        pass

    def on_drop(e):
        # This document-level handler stops the window from navigating away to a dropped file.
        # Dropped files with a full path are loaded by path (read in place). A drop on one of Gradio's own boxes
        # also starts Gradio's upload, which the page reload replaces.
        try:
            paths = [f.get("pywebviewFullPath") for f in e.get("dataTransfer", {}).get("files", [])]
            paths = [p for p in paths if p]
        except Exception as ex:  # noqa: BLE001
            log.warning("could not read the dropped files: %s", ex)
            return
        if paths:
            load_files(paths)

    def on_loaded():
        # The drop guard is registered on every load of the local app's page (one on the splash page would be
        # lost at load_url()).
        try:
            current = window.get_current_url() or ""
            if not (server.url and current.startswith(server.url)):
                return
            doc = window.dom.document
            doc.events.dragover += DOMEventHandler(on_drag, True, False, debounce=500)
            doc.events.drop += DOMEventHandler(on_drop, True, False)
        except Exception as ex:  # noqa: BLE001 - drag-and-drop guard is optional
            log.warning("drag and drop guard not available: %s", ex)

    def on_closed():
        # Mark closing and stop in the background (server.stop() here would block, and a stop during start
        # would be reported as "could not start").
        if state["failure"] is None:
            state["closing"] = True
        t = threading.Thread(target=server.stop, daemon=True)
        t.start()
        state["stopper"] = t

    def startup(win):
        if state["closing"]:
            return
        try:
            url = server.start()
        except ServerError as ex:
            if state["closing"] or ex.cancelled:
                return  # the person closed COAL while it was starting: not an error
            state["failure"] = str(ex)
            try:
                win.destroy()
            except Exception:  # noqa: BLE001 - already closed
                pass
            return
        if state["closing"]:
            return
        if files:  # the files are loaded by path
            want, msg = files_to_load(files, {})
            log.info("opened with files: %s", _names(classify_dropped_files(files)))
            state["loaded"] = want
            win.load_url(query_url(url, want["protein"], want["dna"]))
            if msg:
                _notice("COAL", msg)
            return
        win.load_url(url)

    gui = "edgechromium" if IS_WINDOWS else None  # never fall back to the IE11 engine
    try:
        # create_window is inside the try: any exception opening the window offers the browser fallback.
        window = webview.create_window(WINDOW_TITLE, html=SPLASH, width=1200, height=820,
                                       min_size=(900, 600), text_select=True, zoomable=True,
                                       background_color="#FFFFFF")
        window.events.loaded += on_loaded
        window.events.closed += on_closed
        webview.start(startup, window, gui=gui, private_mode=True,
                      storage_path=os.path.join(data_dir(), "webview"), debug=False)
    except Exception as ex:  # noqa: BLE001 - typically: WebView2 runtime missing
        log.exception("the window could not open")
        if _error_dialog("COAL",
                         "COAL could not open its window. This usually means the Microsoft Edge "
                         "WebView2 Runtime is missing (it is part of Windows 11 and of Edge updates "
                         "on Windows 10).\n\nOpen COAL in your web browser instead? It still runs "
                         f"only on this computer.\n\n(Details: {ex})", offer_browser=True):
            return _run_in_browser(server, instance)
        return 1
    finally:
        t = state["stopper"]
        if t is not None:
            t.join(server.stop_timeout + 10)
        server.stop()
    if state["failure"] and not state["closing"]:
        _error_dialog("COAL could not start", state["failure"])
        return 1
    return 0


def _run_in_browser(server: LocalServer, instance: SingleInstance | None = None) -> int:
    """Fallback when no native window is possible: default browser + a small 'Quit COAL' dialog."""
    import webbrowser
    try:
        url = server.url or server.start()
    except ServerError as ex:
        _error_dialog("COAL could not start", str(ex))
        return 1
    webbrowser.open(url)
    if instance is not None:
        # A second launch in --browser mode reopens the page (and explains "Open with" files).
        loaded = {}

        def on_message(msg):
            if msg.get("files"):  # open the page with the files loaded by path
                want, text = files_to_load(msg["files"], loaded)
                loaded.update(want)
                webbrowser.open(query_url(url, want["protein"], want["dna"]))
                if text:
                    _notice("COAL", text)
                return
            webbrowser.open(url)
        instance.on_message = on_message
    _error_dialog("COAL is running", f"COAL is open in your browser at {url}\n\n"
                  "Press OK to quit COAL.")
    server.stop()
    return 0


# ============================================================== entry point

def _windowed() -> bool:
    """The PyInstaller build is windowed (console=False): nobody sees stderr."""
    return FROZEN and IS_WINDOWS


def parse_args(argv):
    import argparse

    class Parser(argparse.ArgumentParser):
        def error(self, message):
            # The windowed build would exit 2 silently: show a message box (or write to the pipe, for the
            # hidden --serve child).
            if _windowed() or sys.stderr is None:
                text = f"COAL could not understand its command line: {message}"
                if "--serve" in argv:
                    if sys.stdout is not None:
                        print(text, file=sys.stdout, flush=True)
                else:
                    _error_dialog("COAL", text)
                self.exit(2)
            super().error(message)

    ap = Parser(prog="COAL", description="COAL desktop app")
    ap.add_argument("files", nargs="*", help="sequence files to open (from 'Open with')")
    ap.add_argument("--serve", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--port", type=int, default=0, help=argparse.SUPPRESS)
    ap.add_argument("--parent-pid", type=int, default=None, help=argparse.SUPPRESS)
    # --mcp is opt in (the URL is only in local_app.log, on a random port)
    ap.add_argument("--mcp", action="store_true",
                    help="also serve COAL's MCP tools to AI agents on this computer (random port, see local_app.log)")
    ap.add_argument("--no-mcp", action="store_true",
                    help="don't serve COAL's MCP tools (the default for the window)")
    ap.add_argument("--browser", action="store_true", help="open in the default browser instead of a window")
    # The installer's "Connect COAL to Claude Desktop" box runs --connect-claude; uninstall
    # runs --disconnect-claude --quiet.
    ap.add_argument("--connect-claude", action="store_true",
                    help="let Claude Desktop run COAL (adds COAL to Claude Desktop's settings)")
    ap.add_argument("--disconnect-claude", action="store_true", help="remove COAL from Claude Desktop's settings")
    ap.add_argument("--quiet", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--version", action="store_true")
    return ap.parse_args(argv)


def dispatch(argv: list[str]) -> int:
    if argv[:2] == ["-m", "coal"]:
        # UTF-8, as the app reads them (coal_app.core): in cp1252 a file path outside cp1252 would crash the run,
        # and one inside it (é) would arrive garbled
        _utf8_streams()
        return run_engine(argv[2:])
    if argv[:3] == ["-m", "coal_app", "--mcp"]:
        return run_mcp_stdio()
    if argv[:1] == ["-m"]:
        print(f"COAL: only '-m coal' and '-m coal_app --mcp' are supported here, not '-m {' '.join(argv[1:2])}'",
              file=sys.stderr)
        return 2
    args = parse_args(argv)
    if args.connect_claude or args.disconnect_claude:
        return claude_desktop_command(connect=args.connect_claude, quiet=args.quiet)
    if args.version:
        from coal import __version__ as v  # noqa: PLC0415
        print(f"COAL desktop {v}")
        return 0
    if args.serve:
        if not args.port:
            print("COAL: --serve needs --port", file=sys.stderr)
            return 2
        _utf8_streams()
        return serve(args.port, mcp=not args.no_mcp, parent_pid=args.parent_pid)
    return main_window(args)


# ============================================================== Claude Desktop (SHIP-18)

MCP_ARGS = ["-m", "coal_app", "--mcp"]
CLAUDE_CONFIG = "claude_desktop_config.json"


def run_mcp_stdio() -> int:
    """`COAL.exe -m coal_app --mcp`: COAL's MCP server over stdio, started by Claude Desktop (as `python -m coal_app
    --mcp`), so the installed app can be used by Claude."""
    if sys.stdin is None or sys.stdout is None:  # double-clicked or started without pipes (windowed build)
        _error_dialog("COAL", "This part of COAL is started by Claude Desktop, not by hand.\n\nTo connect COAL to "
                              "Claude Desktop, run the COAL installer again and keep \"Connect COAL to Claude "
                              "Desktop\" ticked.")
        return 2
    _utf8_streams()
    from coal_app import mcp_stdio  # noqa: PLC0415
    return mcp_stdio.main()


def claude_config_dirs(env=None) -> list[Path]:
    """Folders where Claude Desktop keeps claude_desktop_config.json on this computer (existing ones only).
    Windows: %APPDATA%\\Claude (installer) and the Microsoft Store package's own copy of it; macOS: Application
    Support/Claude; Linux: ~/.config/Claude."""
    env = os.environ if env is None else env
    found = []
    if env.get("APPDATA"):
        found.append(Path(env["APPDATA"]) / "Claude")
    if env.get("LOCALAPPDATA"):
        found += sorted((Path(env["LOCALAPPDATA"]) / "Packages").glob("Claude_*/LocalCache/Roaming/Claude"))
    home = Path(env.get("HOME") or env.get("USERPROFILE") or Path.home())
    found += [home / "Library" / "Application Support" / "Claude", home / ".config" / "Claude"]
    return [d for d in found if d.is_dir()]


def coal_server_entry(exe: str | None = None) -> dict:
    """The mcpServers entry that starts this COAL: the app itself when frozen, else this Python."""
    return {"command": exe or sys.executable, "args": list(MCP_ARGS)}


def update_claude_config(path: Path, entry: dict | None) -> str:
    """Add (entry) or remove (None) mcpServers.coal in one claude_desktop_config.json, keeping everything else.
    Returns "added", "updated", "unchanged", "removed" or "absent". Raises ValueError if the file is not a JSON
    object (it is then left exactly as it was). A copy of the old file is kept as claude_desktop_config.json.bak-coal."""
    import json
    text = path.read_text(encoding="utf-8-sig") if path.exists() else ""
    if text.strip():
        try:
            data = json.loads(text)
        except json.JSONDecodeError as e:
            raise ValueError(f"{path} is not valid JSON ({e.msg}, line {e.lineno}); COAL did not change it") from e
        if not isinstance(data, dict):
            raise ValueError(f"{path} is not a JSON object; COAL did not change it")
    else:
        data = {}
    servers = data.get("mcpServers")
    if servers is not None and not isinstance(servers, dict):
        raise ValueError(f'"mcpServers" in {path} is not an object; COAL did not change it')
    servers = dict(servers or {})
    old = servers.get("coal")
    if entry is None:
        if old is None:
            return "absent"
        del servers["coal"]
        status = "removed"
    else:
        if old == entry:
            return "unchanged"
        servers["coal"] = entry
        status = "updated" if old is not None else "added"
    data["mcpServers"] = servers
    if text:
        import shutil
        shutil.copyfile(path, path.with_name(path.name + ".bak-coal"))  # byte for byte
    tmp = path.with_name(path.name + ".coal-tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)
    return status


def _info_box(title: str, text: str) -> None:
    if IS_WINDOWS:
        try:
            import ctypes
            ctypes.windll.user32.MessageBoxW(None, text, title, 0x40 | 0x10000)  # information, foreground
            return
        except Exception:  # noqa: BLE001
            pass
    print(f"{title}: {text}")


def claude_desktop_command(connect: bool, quiet: bool = False, env=None, exe: str | None = None) -> int:
    """--connect-claude / --disconnect-claude. Connect writes into every Claude Desktop settings file on this
    computer (creating %APPDATA%\\Claude\\claude_desktop_config.json when Claude Desktop is not installed yet, so it
    works once it is). Disconnect removes only an entry that starts this COAL (or any COAL.exe -m coal_app --mcp)."""
    dirs = claude_config_dirs(env)
    if connect and not dirs:
        e = os.environ if env is None else env
        base = e.get("APPDATA") or str(Path.home() / ".config")
        dirs = [Path(base) / "Claude"]
        dirs[0].mkdir(parents=True, exist_ok=True)
    entry = coal_server_entry(exe)
    results, errors = [], []
    for d in dirs:
        path = d / CLAUDE_CONFIG
        try:
            if connect:
                results.append(update_claude_config(path, entry))
            else:
                import json
                cur = None
                if path.exists():
                    try:
                        cur = (json.loads(path.read_text(encoding="utf-8-sig") or "{}").get("mcpServers") or {}).get("coal")
                    except (ValueError, AttributeError):
                        cur = None
                if isinstance(cur, dict) and cur.get("args") == MCP_ARGS and \
                        Path(str(cur.get("command", ""))).name.lower() == Path(entry["command"]).name.lower():
                    results.append(update_claude_config(path, None))
        except (OSError, ValueError) as e:
            errors.append(str(e))
    if quiet:
        return 1 if errors else 0
    if errors:
        _error_dialog("COAL", "COAL could not connect to Claude Desktop:\n\n" + "\n".join(errors))
        return 1
    if connect:
        _info_box("COAL is connected to Claude Desktop",
                  "Claude Desktop can now run COAL for you.\n\n"
                  "If Claude Desktop is open, quit it completely (right-click its icon by the clock, then Quit) and "
                  "open it again.\n\n"
                  "Then ask Claude, for example: \"Use COAL to make a codon alignment from the files in "
                  "C:\\Users\\you\\Documents\\genes\".")
    else:
        _info_box("COAL", "COAL was removed from Claude Desktop." if "removed" in results
                  else "COAL was not connected to Claude Desktop.")
    return 0


def _utf8_streams() -> None:
    """The --serve child writes to a pipe; in the frozen app (PyInstaller) that pipe uses the ANSI code page
    (cp1252) whatever PYTHONIOENCODING says, so Gradio's own start-up line with an emoji ('🔨 Launching MCP
    server') would raise UnicodeEncodeError and kill the local app.
    The parent decodes the pipe as UTF-8 (LocalServer), so use UTF-8 here and never fail on a character."""
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            try:
                reconfigure(encoding="utf-8", errors="replace")
            except (OSError, ValueError):
                pass


def _hold_app_mutex():
    """Windows: create a named mutex that exists for as long as COAL runs. The handle is never
    closed: Windows closes it when the process exits, which is when the mutex disappears for
    Inno Setup's AppMutex check. Untested here (Linux)."""
    if not IS_WINDOWS:
        return None
    try:
        import ctypes
        from ctypes import wintypes
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        # HANDLE restype: the default int would truncate it to 32 bits.
        k32.CreateMutexW.restype = wintypes.HANDLE
        k32.CreateMutexW.argtypes = [ctypes.c_void_p, wintypes.BOOL, wintypes.LPCWSTR]
        return k32.CreateMutexW(None, False, APP_MUTEX)
    except Exception as e:  # noqa: BLE001
        log.warning("could not create the app mutex: %s", e)
        return None


INSTANCE_WAIT = 3.0  # seconds to wait for a lock still being released / a primary still starting


def claim_instance(instance: SingleInstance, message: dict, wait: float = INSTANCE_WAIT) -> str:
    """'primary' (we hold the lock), 'notified' (the running COAL got `message`) or 'unreachable'.

    Retries for `wait` seconds: on Windows a crashed holder's lock is released asynchronously, and
    a new primary takes the lock before writing instance.json."""
    deadline = time.monotonic() + wait
    while True:
        if instance.acquire():
            return "primary"
        if instance.notify_primary(message, timeout=2.0):
            return "notified"
        if time.monotonic() >= deadline:
            return "unreachable"
        time.sleep(0.25)


def main_window(args) -> int:
    instance = SingleInstance(data_dir())
    files = [os.path.abspath(f) for f in args.files]
    # Only the lock holder touches the logs (rotating before this check would rotate the running
    # COAL's log); a second launch writes none.
    role = claim_instance(instance, {"action": "focus", "files": files})
    if role == "notified":
        return 0
    if role == "unreachable":
        _error_dialog("COAL", "COAL seems to be already running but did not respond. "
                      "Close it from the taskbar (or Task Manager) and try again.")
        return 1
    server = None
    try:
        log_path = setup_logging()
        local_log = os.path.join(os.path.dirname(log_path), "local_app.log")
        rotate_log(local_log, keep=3)  # rotated per launch
        if IS_WINDOWS and not _hold_app_mutex():
            log.warning("app mutex not created: setup/uninstall won't see that COAL is open")
        log.info("COAL desktop starting (frozen=%s, python %s)", FROZEN, sys.version.split()[0])
        server = LocalServer(command=lambda port: default_server_command(port, mcp=args.mcp and not args.no_mcp),
                             log_path=local_log)
        if args.browser:
            return _run_in_browser(server, instance)
        try:
            import webview  # noqa: F401
        except ImportError:
            _error_dialog("COAL", "The window component (pywebview) is not installed. "
                          "Install desktop/requirements.txt, or start with --browser.")
            return 1
        return run_window(server, instance, files)
    except Exception:  # noqa: BLE001
        # Logged and explained (in the windowed build a propagated error is an "Unhandled exception" box).
        log.exception("unexpected error")
        _error_dialog("COAL", "COAL stopped because of an unexpected error. Details are in "
                      f"{os.path.join(data_dir(), 'logs')}.")
        return 1
    finally:
        if server is not None:
            server.stop()
        instance.release()


def main(argv=None) -> int:
    multiprocessing.freeze_support()  # the engine's worker processes in a PyInstaller build
    setup_source_paths()
    return dispatch(list(sys.argv[1:] if argv is None else argv))


if __name__ == "__main__":
    # First thing in the frozen entry point, before anything else runs, as PyInstaller recommends.
    multiprocessing.freeze_support()
    sys.exit(main())
