"""_resource_monitor.py -- reusable peak-RSS / CPU sampler for benchmark runs.

Extracted from the battle-tested ResourceSampler in run_dss_replication.py so
that every benchmark runner (phi-sweep, simulator, head-to-heads) can capture
**per-tool peak memory** with one helper, not just wall/CPU time.

Why a process-tree poller (not /usr/bin/time -v): methylKit's
``calculateDiffMeth(mc.cores=N)`` forks worker processes. ``/usr/bin/time``'s
``ru_maxrss`` only reports the main child's peak and undercounts the forked
workers. We poll the full psutil tree (root + recursive children) every
``interval`` seconds and sum RSS across all live processes, tracking the peak.

Public API:
    run_subprocess_monitored(cmd, *, interval=0.2, cwd=None, env=None,
                             stdout=None, stderr=None) -> dict
        Run ``cmd`` as a subprocess while sampling its process tree. Returns a
        dict with returncode, wall_s, cpu_s (parent process_time), and the
        memory/CPU summary (rss_peak_mb, rss_mean_mb, uss_peak_mb,
        cpu_percent_peak, threads_peak, num_processes_peak, samples_collected),
        plus ``samples`` (the raw per-interval time series).

    ResourceSampler(root_proc, interval)
        The low-level threading.Thread poller, if you need to wrap an existing
        psutil.Process (e.g. the current interpreter) yourself.

    summarize_samples(samples, wall_seconds) -> dict
        Pure-python peak/mean reduction over the sample list.
"""

from __future__ import annotations

import subprocess
import sys
import threading
import time
from datetime import datetime

import psutil

DEFAULT_INTERVAL_S = 0.2


class ResourceSampler(threading.Thread):
    """Polls a psutil.Process tree for memory + CPU every ``interval`` seconds.

    ``cpu_percent`` is normalized to a single core (a fully-utilized 8-core
    machine reports up to 800%). RSS is summed across the parent + all
    children; threads are summed too. USS (unique set size) is summed where
    available (some processes deny access; falls back to None for that sample).
    """

    def __init__(self, root_proc: psutil.Process, interval: float = DEFAULT_INTERVAL_S):
        super().__init__(daemon=True)
        self.root = root_proc
        self.interval = interval
        self.samples: list[dict] = []
        self.stop_flag = threading.Event()
        self.t0 = time.time()
        try:
            self.root.cpu_percent(interval=None)  # prime
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass

    def _walk(self) -> list[psutil.Process]:
        procs = [self.root]
        try:
            procs.extend(self.root.children(recursive=True))
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
        return procs

    def _sample(self) -> dict | None:
        procs = self._walk()
        if not procs:
            return None
        rss = uss = 0
        uss_avail = False
        threads = 0
        cpu = 0.0
        n_alive = 0
        for p in procs:
            try:
                with p.oneshot():
                    rss += p.memory_info().rss
                    try:
                        uss += p.memory_full_info().uss
                        uss_avail = True
                    except (psutil.AccessDenied, AttributeError):
                        pass
                    threads += p.num_threads()
                    cpu += p.cpu_percent(interval=None)
                    n_alive += 1
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        if n_alive == 0:
            return None
        return dict(
            timestamp=datetime.now().isoformat(timespec="seconds"),
            elapsed_s=round(time.time() - self.t0, 2),
            rss_mb=round(rss / 1024**2, 2),
            uss_mb=round(uss / 1024**2, 2) if uss_avail else None,
            cpu_percent=round(cpu, 1),
            num_threads=threads,
            num_processes=n_alive,
        )

    def run(self) -> None:
        # Sample immediately (don't sleep first) so very short runs still get
        # at least one reading near their peak.
        while not self.stop_flag.is_set():
            s = self._sample()
            if s is None:
                break
            self.samples.append(s)
            time.sleep(self.interval)


def summarize_samples(samples: list[dict], wall_seconds: float) -> dict:
    """Pure-python peak/mean reduction over the per-interval samples."""
    if not samples:
        return dict(samples_collected=0, wall_s=round(wall_seconds, 2),
                    rss_peak_mb=None, rss_mean_mb=None, uss_peak_mb=None,
                    cpu_percent_peak=None, cpu_percent_mean=None,
                    threads_peak=None, num_processes_peak=None)
    rss = [s["rss_mb"] for s in samples if s["rss_mb"] is not None]
    uss = [s["uss_mb"] for s in samples if s.get("uss_mb") is not None]
    cpu = [s["cpu_percent"] for s in samples if s["cpu_percent"] is not None]
    thr = [s["num_threads"] for s in samples if s.get("num_threads") is not None]
    nproc = [s["num_processes"] for s in samples if s.get("num_processes") is not None]
    cpu_mean = (sum(cpu) / len(cpu)) if cpu else 0.0
    core_seconds = (cpu_mean / 100.0) * wall_seconds
    return dict(
        samples_collected=len(samples),
        sample_interval_s=None,  # filled by caller
        wall_s=round(wall_seconds, 2),
        rss_peak_mb=round(max(rss), 1) if rss else None,
        rss_mean_mb=round(sum(rss) / len(rss), 1) if rss else None,
        uss_peak_mb=round(max(uss), 1) if uss else None,
        cpu_percent_peak=round(max(cpu), 1) if cpu else None,
        cpu_percent_mean=round(cpu_mean, 1),
        approx_core_seconds=round(core_seconds, 1),
        threads_peak=int(max(thr)) if thr else None,
        num_processes_peak=int(max(nproc)) if nproc else None,
    )


def run_subprocess_monitored(
    cmd: list[str], *, interval: float = DEFAULT_INTERVAL_S,
    cwd=None, env=None, stdout=None, stderr=None,
) -> dict:
    """Run ``cmd`` as a subprocess, sampling its process tree for peak RSS/CPU.

    Returns a dict: returncode, wall_s, cpu_s (parent process_time delta), the
    summarize_samples() fields, and ``samples`` (raw time series). ``stdout`` /
    ``stderr`` default to DEVNULL; pass subprocess.PIPE or a file handle to
    capture/redirect.
    """
    if stdout is None:
        stdout = subprocess.DEVNULL
    if stderr is None:
        stderr = subprocess.DEVNULL
    t0_wall = time.time()
    t0_proc = time.process_time()
    proc = subprocess.Popen([str(c) for c in cmd], cwd=cwd, env=env,
                            stdout=stdout, stderr=stderr)
    try:
        ps = psutil.Process(proc.pid)
    except psutil.NoSuchProcess:
        proc.wait()
        res = summarize_samples([], time.time() - t0_wall)
        res.update(returncode=proc.returncode, cpu_s=0.0, samples=[])
        return res
    sampler = ResourceSampler(ps, interval)
    sampler.start()
    rc = proc.wait()
    sampler.stop_flag.set()
    sampler.join(timeout=5)
    wall_s = time.time() - t0_wall
    cpu_s = time.process_time() - t0_proc
    res = summarize_samples(sampler.samples, wall_s)
    res["sample_interval_s"] = interval
    res.update(returncode=rc, cpu_s=round(cpu_s, 3), samples=sampler.samples)
    return res


if __name__ == "__main__":
    # Smoke test: monitor any command passed on argv.
    if len(sys.argv) < 2:
        print("usage: python _resource_monitor.py <cmd> [args...]", file=sys.stderr)
        sys.exit(2)
    res = run_subprocess_monitored(sys.argv[1:], stdout=None, stderr=None)
    for k in ("returncode", "wall_s", "rss_peak_mb", "uss_peak_mb",
              "cpu_percent_peak", "num_processes_peak", "samples_collected"):
        print(f"{k:20s} {res.get(k)}")
