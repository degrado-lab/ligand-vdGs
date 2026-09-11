"""A bucket clustered as phased tasks must write exactly what the one-task path
writes, and one failing bucket must not take the others down with it."""
import concurrent.futures
import multiprocessing as mp
import os
import tempfile
import unittest
from unittest import mock

import numpy as np

from ligand_vdgs.functions.clus_helpers import save_columns
from ligand_vdgs.functions.compute_profile import ComputeProfile
from ligand_vdgs.generate_vdgs import clus_and_deduplicate_vdgs as C


def _bucket_columns(rng, n, n_cg=4, n_res=1, flank=2, n_centers=6):
    n_total = n_cg + 3 * n_res
    n_flank = n_res * (2 * flank + 1)
    centers = rng.normal(scale=4.0, size=(n_centers, n_total, 3))
    per = -(-n // n_centers)
    cgvdmbb = np.concatenate([c + rng.normal(scale=0.2, size=(per, n_total, 3))
                              for c in centers])[:n].astype(np.float32)
    tokens = np.array(["ALA", "GLY", "SER", "-", "!"])
    seq = tokens[rng.integers(0, 5, size=(n, n_flank))].astype("U4")
    seq[:, flank::(2 * flank + 1)] = "vdm"
    return {
        "cgvdmbb": cgvdmbb,
        "flank_ca": rng.normal(scale=3.0, size=(n, n_flank, 3)).astype(np.float32),
        "flank_seq": seq,
        "biounit": np.array([f"1a{i % 23:02d}" + ("_2" if i % 7 == 0 else "")
                             for i in range(n)], dtype="U32"),
        "scrr_seg": np.full((n, n_res), "", dtype="U8"),
        "scrr_chain": np.full((n, n_res), "A", dtype="U2"),
        "scrr_resnum": np.tile(np.arange(n_res, dtype=np.int32) + 10, (n, 1)),
        "scrr_resname": np.full((n, n_res), "ALA", dtype="U4"),
        "cg_names": np.tile(np.array([f"C{k}" for k in range(n_cg)], dtype="U4"), (n, 1)),
        "cg_elements": np.full((n, n_cg), "C", dtype="U2"),
        "cg_seg": np.full(n, "", dtype="U8"),
        "cg_chain": np.full(n, "L", dtype="U2"),
        "cg_resnum": np.full(n, 1, dtype=np.int32),
        "cg_resname": np.full(n, "LIG", dtype="U4"),
        "slot_flags": np.zeros((n, n_res), dtype=np.int8),
        "quality": np.ones((n, 4), dtype=np.float32),
        "vdm_o": np.full((n, n_res, 3), np.nan, dtype=np.float32),
    }


def _run(tmp, tag, split_min, buckets_spec, break_bucket=None):
    """Cluster `buckets_spec` ({aa_label: columns}) into <tmp>/<tag>; returns
    (failed labels, profile counters, library dir)."""
    lib = os.path.join(tmp, tag)
    stream = os.path.join(tmp, f"stream_{tag}", "1")
    buckets = []
    for label, cols in buckets_spec.items():
        bucket_dir = os.path.join(stream, label)
        save_columns(bucket_dir, cols)
        if label == break_bucket:
            os.remove(os.path.join(bucket_dir, "cgvdmbb.npy"))
        buckets.append(C.Bucket(key=(1, label), size_subset=1,
                                aa_parts=tuple(label.split("_")),
                                n=len(cols["biounit"]), bucket_dir=bucket_dir))
    buckets.sort(key=lambda b: -b.n)
    run = C.Run(cg_automorphisms=((0, 1, 2, 3), (0, 2, 1, 3)), seq_sim_thresh=0.4,
                vdglib_dir=lib, logfile=os.path.join(tmp, f"log_{tag}"),
                parent_pdb_dir="/db")
    profile = ComputeProfile(enabled=True)
    pool = concurrent.futures.ProcessPoolExecutor(
        max_workers=2, mp_context=mp.get_context("spawn"))
    try:
        with mock.patch.object(C, "_SPLIT_MIN_RECORDS", split_min), \
                mock.patch.object(C, "_BLOCKS_PER_PROC", 3):
            failed = C._cluster_buckets(pool, 2, run, buckets, profile)
    finally:
        pool.shutdown(wait=True, cancel_futures=True)
    return failed, profile.counters, lib


def _load(lib, label):
    with np.load(os.path.join(lib, "nr_vdgs", "1", f"{label}.npz")) as z:
        return {k: z[k] for k in z.files}


class PhasedBucketTests(unittest.TestCase):
    def test_split_bucket_writes_the_same_library_as_the_one_task_path(self):
        rng = np.random.default_rng(3)
        spec = {"ALA": _bucket_columns(rng, 150), "GLY": _bucket_columns(rng, 20)}
        with tempfile.TemporaryDirectory() as tmp:
            failed_a, counters_a, lib_a = _run(tmp, "split", 40, spec)
            failed_b, counters_b, lib_b = _run(tmp, "single", 10_000, spec)
            self.assertEqual((failed_a, failed_b), ([], []))
            self.assertEqual(counters_a["buckets_split"], 1)
            self.assertEqual(counters_b["buckets_split"], 0)
            # Both paths ran the same stages over the same pairs.
            self.assertEqual(counters_a["stage1.edges"], counters_b["stage1.edges"])
            self.assertEqual(counters_a["stage1.stage1_clusters"],
                             counters_b["stage1.stage1_clusters"])
            for label in spec:
                a, b = _load(lib_a, label), _load(lib_b, label)
                self.assertEqual(sorted(a), sorted(b))
                for key in a:
                    np.testing.assert_array_equal(a[key], b[key], err_msg=f"{label}:{key}")
                self.assertEqual(int(a["cluster_size"].sum()), len(spec[label]["biounit"]))
                self.assertGreater(len(a["cluster_id"]), 1)
            # The scratch intermediates of the split bucket are gone.
            self.assertEqual(sorted(os.listdir(os.path.join(tmp, "stream_split", "1", "ALA"))),
                             ["cgvdmbb.npy", "columns.npz"])

    def test_one_failing_bucket_is_marked_and_the_rest_still_finish(self):
        rng = np.random.default_rng(4)
        spec = {"ALA": _bucket_columns(rng, 120), "GLY": _bucket_columns(rng, 30),
                "SER": _bucket_columns(rng, 90)}
        with tempfile.TemporaryDirectory() as tmp:
            failed, _counters, lib = _run(tmp, "fail", 40, spec, break_bucket="SER")
            self.assertEqual(failed, ["SER"])
            self.assertEqual(sorted(os.listdir(os.path.join(lib, "nr_vdgs", "1"))),
                             ["ALA.npz", "GLY.npz", "SER.FAILED"])
            marker = open(os.path.join(lib, "nr_vdgs", "1", "SER.FAILED")).read()
            self.assertIn("task stage1_block", marker)


if __name__ == "__main__":
    unittest.main()
