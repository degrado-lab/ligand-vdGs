"""--workers changes wall time only: post-build bucket errors are identical for 1 and 3 workers.

Falsifier: a pool that drops or misattributes a worker's errors reports a different error list
(or none) for the fragment whose bucket lacks an H-class field.
"""
import os
import tempfile
import unittest
import numpy as np
from ligand_vdgs.functions import clus_helpers, utils, vdg_npz_utils
from ligand_vdgs.generate_vdgs import clus_and_deduplicate_vdgs as clus
from ligand_vdgs.generate_vdgs.post_build_buckets import inspect_library
from tests.test_bucket_schema_pass import _record

# key -> (nr rows, CG automorphisms)
FRAGS = {"OPO": (10, [(0, 1, 2), (2, 1, 0)]), "OPS": (12, [(0, 1, 2)]),
         "[O-]P[O-]": (14, [(0, 1, 2), (2, 1, 0)])}
FIELD = "nr_cg_num_h"

def _build(lib):
    for n, (key, (n_rows, perms)) in enumerate(FRAGS.items()):
        frag_dir = os.path.join(lib, utils.smiles_to_filename(key))
        clus._write_bucket_npz(frag_dir, 2, "neg", ("ALA", "GLY"), clus_helpers.records_to_columns(
            [_record(cg_coords=xyz, cg_num_h=[i % 2, 0, (i + 1) % 2], biounit=f"{n}x{i:02d}_1") for i, xyz
             in enumerate(np.random.default_rng(n).uniform(-4, 5, (n_rows, 3, 3)).astype(np.float32))]),
                               [clus.Subgroup(1, 1, i, np.array([i], dtype=np.int32), 0.25)
                                for i in range(n_rows)], "/db")
        vdg_npz_utils.write_cg_symmetry(frag_dir, key, perms)
        with open(os.path.join(frag_dir, f"{utils.smiles_to_filename(key)}_log"), "w") as fh:
            fh.write("Job completed.\n")

def _bucket(lib, key): return vdg_npz_utils.vdg_npz_path(lib, utils.smiles_to_filename(key), 2, "neg", "ALA_GLY")

def _break(lib, key):
    """Drop one H-class field from the fragment's bucket."""
    with np.load(_bucket(lib, key)) as z:
        np.savez_compressed(_bucket(lib, key), **{k: z[k] for k in z.files if k != FIELD})

class WorkersInvariance(unittest.TestCase):
    def test_inspect_library_identical_across_worker_counts(self):
        with tempfile.TemporaryDirectory() as lib:
            _build(lib)
            frags = sorted(utils.smiles_to_filename(k) for k in FRAGS)
            self.assertEqual(inspect_library(lib, frags, workers=3), (len(FRAGS), []))
            _break(lib, "OPS")
            serial, parallel = inspect_library(lib, frags, workers=1), inspect_library(lib, frags, workers=3)
            self.assertEqual(serial, parallel)
            self.assertEqual(serial[0], len(FRAGS))
            self.assertTrue(serial[1] and all(_bucket(lib, "OPS") in e and FIELD in e for e in serial[1]),
                            serial[1])

if __name__ == "__main__":
    unittest.main()
