#!/usr/bin/env python3
"""Windows only: run a command inside a Job Object with a memory limit (the `ulimit -v` of the stress test).

    python tests/tools/stress_winjob.py --limit-mb 256 -- coal pep.fasta nuc.fasta -o out.fa

Creates a job with JOB_OBJECT_LIMIT_PROCESS_MEMORY and JOB_OBJECT_LIMIT_JOB_MEMORY set to --limit-mb (committed
memory, the whole process tree), starts the command suspended, assigns it to the job, resumes it and waits. Prints
the exit code, wall time and the job's peak memory. Standard library (ctypes) only.
"""
from __future__ import annotations

import argparse
import ctypes
import subprocess
import sys
import time
from ctypes import wintypes

k32 = ctypes.WinDLL("kernel32", use_last_error=True)

JobObjectExtendedLimitInformation = 9
JOB_OBJECT_LIMIT_PROCESS_MEMORY = 0x00000100
JOB_OBJECT_LIMIT_JOB_MEMORY = 0x00000200
JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE = 0x00002000
CREATE_SUSPENDED = 0x00000004


class IO_COUNTERS(ctypes.Structure):
    _fields_ = [(n, ctypes.c_ulonglong) for n in ("ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
                                                  "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]


class JOBOBJECT_BASIC_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [("PerProcessUserTimeLimit", ctypes.c_longlong), ("PerJobUserTimeLimit", ctypes.c_longlong),
                ("LimitFlags", wintypes.DWORD), ("MinimumWorkingSetSize", ctypes.c_size_t),
                ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", wintypes.DWORD),
                ("Affinity", ctypes.c_size_t), ("PriorityClass", wintypes.DWORD), ("SchedulingClass", wintypes.DWORD)]


class JOBOBJECT_EXTENDED_LIMIT_INFORMATION(ctypes.Structure):
    _fields_ = [("BasicLimitInformation", JOBOBJECT_BASIC_LIMIT_INFORMATION), ("IoInfo", IO_COUNTERS),
                ("ProcessMemoryLimit", ctypes.c_size_t), ("JobMemoryLimit", ctypes.c_size_t),
                ("PeakProcessMemoryUsed", ctypes.c_size_t), ("PeakJobMemoryUsed", ctypes.c_size_t)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--limit-mb", type=int, default=256)
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args(argv)
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    job = k32.CreateJobObjectW(None, None)
    if not job:
        raise ctypes.WinError(ctypes.get_last_error())
    info = JOBOBJECT_EXTENDED_LIMIT_INFORMATION()
    info.BasicLimitInformation.LimitFlags = (JOB_OBJECT_LIMIT_PROCESS_MEMORY | JOB_OBJECT_LIMIT_JOB_MEMORY
                                             | JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE)
    info.ProcessMemoryLimit = a.limit_mb * 1024 * 1024
    info.JobMemoryLimit = a.limit_mb * 1024 * 1024
    if not k32.SetInformationJobObject(job, JobObjectExtendedLimitInformation, ctypes.byref(info), ctypes.sizeof(info)):
        raise ctypes.WinError(ctypes.get_last_error())
    t0 = time.monotonic()
    p = subprocess.Popen(cmd, creationflags=CREATE_SUSPENDED)
    hproc = wintypes.HANDLE(p._handle)
    if not k32.AssignProcessToJobObject(job, hproc):
        raise ctypes.WinError(ctypes.get_last_error())
    # resume the main thread: find it through the toolhelp snapshot
    TH32CS_SNAPTHREAD = 0x4

    class THREADENTRY32(ctypes.Structure):
        _fields_ = [("dwSize", wintypes.DWORD), ("cntUsage", wintypes.DWORD), ("th32ThreadID", wintypes.DWORD),
                    ("th32OwnerProcessID", wintypes.DWORD), ("tpBasePri", wintypes.LONG), ("tpDeltaPri", wintypes.LONG),
                    ("dwFlags", wintypes.DWORD)]
    snap = k32.CreateToolhelp32Snapshot(TH32CS_SNAPTHREAD, 0)
    te = THREADENTRY32()
    te.dwSize = ctypes.sizeof(te)
    ok = k32.Thread32First(snap, ctypes.byref(te))
    while ok:
        if te.th32OwnerProcessID == p.pid:
            h = k32.OpenThread(0x0002, False, te.th32ThreadID)  # THREAD_SUSPEND_RESUME
            k32.ResumeThread(h)
            k32.CloseHandle(h)
        ok = k32.Thread32Next(snap, ctypes.byref(te))
    k32.CloseHandle(snap)
    rc = p.wait()
    dt = time.monotonic() - t0
    out = JOBOBJECT_EXTENDED_LIMIT_INFORMATION()
    k32.QueryInformationJobObject(job, JobObjectExtendedLimitInformation, ctypes.byref(out), ctypes.sizeof(out), None)
    print(f"winjob: exit {rc} ({rc & 0xFFFFFFFF:#x}) in {dt:.2f}s, peak job memory {out.PeakJobMemoryUsed / 1e6:.0f} MB, "
          f"peak process {out.PeakProcessMemoryUsed / 1e6:.0f} MB, limit {a.limit_mb} MB", flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
