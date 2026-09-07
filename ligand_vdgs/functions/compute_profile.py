"""Optional per-phase timing and counters for a vdG generation run.

On by default (disable with ``--no-profile-compute``); one small JSON sidecar is
written per fragment and the ordinary logfile gets a single line pointing at it,
so a production log stays readable while a slow filesystem node or a pathological
bucket is still diagnosable after the fact.

Disabled, every method is a no-op: the timing calls can sit in hot paths without
a flag check at each call site.
"""
import json
import os
import time
from contextlib import contextmanager


class ComputeProfile:
    """Phase wall/CPU times, counters, and nested sections."""

    __slots__ = ("enabled", "phases", "counters", "sections", "_meta")

    def __init__(self, enabled=False, **meta):
        self.enabled = bool(enabled)
        self.phases = {}
        self.counters = {}
        self.sections = {}
        self._meta = dict(meta)

    @contextmanager
    def phase(self, name):
        """Time a block, accumulating if the name repeats."""
        if not self.enabled:
            yield
            return
        wall, cpu = time.perf_counter(), time.process_time()
        try:
            yield
        finally:
            entry = self.phases.setdefault(name, {"wall_s": 0.0, "cpu_s": 0.0, "calls": 0})
            entry["wall_s"] += time.perf_counter() - wall
            entry["cpu_s"] += time.process_time() - cpu
            entry["calls"] += 1

    def add(self, name, n=1):
        if self.enabled:
            self.counters[name] = self.counters.get(name, 0) + n

    def merge(self, mapping, prefix=""):
        """Fold a worker's counter dict into this one."""
        if not self.enabled or not mapping:
            return
        for key, value in mapping.items():
            self.add(f"{prefix}{key}", value)

    def set(self, name, value):
        if self.enabled:
            self.counters[name] = value

    def max(self, name, value):
        """Keep the largest value seen for a counter (peaks, not totals)."""
        if self.enabled:
            prior = self.counters.get(name)
            self.counters[name] = value if prior is None else max(prior, value)

    def section(self, name, payload):
        """Attach a nested block, e.g. a subprocess's own profile."""
        if self.enabled and payload:
            self.sections[name] = payload

    def record_peak_rss(self, prefix=""):
        """Record peak RSS of this process and of its reaped children, in MiB.

        ``ru_maxrss`` for children is the largest single child's peak, not the
        sum, and only counts children already waited for -- so this must be
        called after a pool closes. It is the right number for sizing a job's
        ``mem_free`` only once multiplied by how many of those children run at
        once; the profile records the two separately rather than guessing.
        """
        if not self.enabled:
            return
        import resource
        # Linux reports ru_maxrss in KiB.
        self.max(f"{prefix}rss_peak_self_mb",
                 round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1))
        self.max(f"{prefix}rss_peak_child_mb",
                 round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, 1))

    def record_dir_size(self, name, path):
        """Record the on-disk size of a scratch tree in MiB, as a peak."""
        if not self.enabled or not path or not os.path.isdir(path):
            return 0.0
        total = 0
        for root, _dirs, files in os.walk(path):
            for fname in files:
                try:
                    total += os.path.getsize(os.path.join(root, fname))
                except OSError:
                    pass
        return self.record_size_bytes(name, total)

    def record_size_bytes(self, name, total_bytes):
        """Record an already-summed byte total in MiB, as a peak. For callers that
        are walking the tree anyway and should not walk it a second time."""
        if not self.enabled:
            return 0.0
        mib = round(total_bytes / (1024 * 1024), 1)
        self.max(name, mib)
        return mib

    def to_dict(self):
        return {"meta": self._meta, "phases": self.phases,
                "counters": self.counters, "sections": self.sections}

    def write(self, path):
        """Write the sidecar. Returns the path, or None when disabled."""
        if not self.enabled or not path:
            return None
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        with open(path, "w") as handle:
            json.dump(self.to_dict(), handle, indent=2, sort_keys=True)
        return path

    @staticmethod
    def load(path):
        """Read a sidecar written by another process; {} if absent or corrupt."""
        if not path or not os.path.isfile(path):
            return {}
        try:
            with open(path) as handle:
                return json.load(handle)
        except (OSError, ValueError):
            return {}
