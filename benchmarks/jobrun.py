"""Run one command inside a Windows Job Object and return its time, memory and CPU accounting.

Timing and memory come from the same run: wall time is time.perf_counter around the process's
lifetime; memory is JOBOBJECT_EXTENDED_LIMIT_INFORMATION.PeakJobMemoryUsed (peak *committed* memory of the
whole process tree, not RSS); user/kernel time from JOBOBJECT_BASIC_ACCOUNTING_INFORMATION. The job has a
memory cap (JobMemoryLimit): allocations beyond it fail, which is reported as status OOM.
Falls back to psutil peak_wset polling (flagged mem_source="psutil_peak_wset") if pywin32 is unavailable.
"""
from __future__ import annotations

import ctypes
import os
import subprocess
import time
from dataclasses import dataclass

try:
    import win32api
    import win32con
    import win32job
except ImportError:  # pragma: no cover
    win32job = None

CREATE_SUSPENDED = 0x00000004
DEFAULT_MEM_CAP = 24 * 2**30  # agreed 24 GB artificial cap
OOM_MARKERS = (b"Out of memory", b"MemoryError", b"bad_alloc", b"not enough memory", b"Cannot allocate")


@dataclass
class RunResult:
    status: str          # OK | ERROR | TIMEOUT | OOM
    exit_code: int | None
    wall_s: float
    user_s: float | None
    kernel_s: float | None
    peak_mem_mb: float | None
    mem_source: str
    reason: str = ""


def clean_env(extra: dict | None = None) -> dict:
    env = {k: v for k, v in os.environ.items() if k != "COAL_REF_SRC"}
    env.update(extra or {})
    assert "COAL_REF_SRC" not in env
    return env


def _stderr_tail(path, n=300) -> bytes:
    if not path or not os.path.exists(path):
        return b""
    with open(path, "rb") as fh:
        fh.seek(max(0, os.path.getsize(path) - n))
        return fh.read()


def run(cmd, stdout_path, stderr_path, timeout: float, mem_cap: int = DEFAULT_MEM_CAP,
        env: dict | None = None) -> RunResult:
    env = env if env is not None else clean_env()
    assert "COAL_REF_SRC" not in env, "COAL_REF_SRC must be unset"
    if win32job is None:
        return _run_psutil(cmd, stdout_path, stderr_path, timeout, env)

    job = win32job.CreateJobObject(None, "")
    info = win32job.QueryInformationJobObject(job, win32job.JobObjectExtendedLimitInformation)
    info["BasicLimitInformation"]["LimitFlags"] = (win32job.JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
                                                  | win32job.JOB_OBJECT_LIMIT_JOB_MEMORY)
    info["JobMemoryLimit"] = mem_cap
    win32job.SetInformationJobObject(job, win32job.JobObjectExtendedLimitInformation, info)

    timed_out = False
    with open(stdout_path, "wb") as out, open(stderr_path, "wb") as err:
        t0 = time.perf_counter()
        p = subprocess.Popen(cmd, stdout=out, stderr=err, stdin=subprocess.DEVNULL, env=env,
                             creationflags=CREATE_SUSPENDED)
        hproc = win32api.OpenProcess(win32con.PROCESS_ALL_ACCESS, False, p.pid)
        win32job.AssignProcessToJobObject(job, hproc)
        ctypes.windll.ntdll.NtResumeProcess(int(hproc))  # start only once it is inside the job
        try:
            rc = p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            win32job.TerminateJobObject(job, 1)
            rc = p.wait()
        wall = time.perf_counter() - t0
    # let stray children (multiprocessing workers) finish so accounting is complete
    for _ in range(500):
        if win32job.QueryInformationJobObject(job, win32job.JobObjectBasicAccountingInformation)["ActiveProcesses"] == 0:
            break
        time.sleep(0.01)
    acc = win32job.QueryInformationJobObject(job, win32job.JobObjectBasicAccountingInformation)
    ext = win32job.QueryInformationJobObject(job, win32job.JobObjectExtendedLimitInformation)
    win32job.TerminateJobObject(job, 1)
    hproc.Close()
    job.Close()
    peak = ext["PeakJobMemoryUsed"]
    res = RunResult("OK", rc, wall, acc["TotalUserTime"] / 1e7, acc["TotalKernelTime"] / 1e7,
                    peak / 2**20, "job_peak_committed")
    tail = _stderr_tail(stderr_path)
    if timed_out:
        res.status, res.reason = "TIMEOUT", f"timeout={timeout:g}s"
    elif rc != 0 and (peak >= 0.97 * mem_cap or any(m in tail for m in OOM_MARKERS)):
        res.status, res.reason = "OOM", f"cap={mem_cap / 2**30:g}GB; peak={peak / 2**20:.0f}MB"
    elif rc != 0:
        last = tail.decode("utf-8", "replace").strip().splitlines()
        res.status, res.reason = "ERROR", f"exit={rc}: " + (last[-1][:200] if last else "")
    return res


def _run_psutil(cmd, stdout_path, stderr_path, timeout, env) -> RunResult:  # pragma: no cover
    import psutil
    with open(stdout_path, "wb") as out, open(stderr_path, "wb") as err:
        t0 = time.perf_counter()
        p = psutil.Popen(cmd, stdout=out, stderr=err, env=env)
        peak, status = 0, "OK"
        while p.poll() is None:
            try:
                peak = max(peak, p.memory_info().peak_wset)
            except psutil.Error:
                pass
            if time.perf_counter() - t0 > timeout:
                p.kill()
                status = "TIMEOUT"
                break
            time.sleep(0.005)
        rc = p.wait()
        wall = time.perf_counter() - t0
    if status == "OK" and rc != 0:
        status = "ERROR"
    return RunResult(status, rc, wall, None, None, peak / 2**20, "psutil_peak_wset", "")
