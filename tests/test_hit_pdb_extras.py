"""write_vdg_hit_pdbs adds the same parent-derived display atoms as the library
writer, and degrades the same way when no parent PDB mirror is reachable.

The hit's superposition is the hit finder's, over CG + vdM N/CA/C; the extras
ride the stored rigid transform like every other atom.
"""
import os
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np
import prody as pr

from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.tools import write_vdg_hit_pdbs as whp


FRAG = "CCO"
BIOUNIT = "1abc"
# ALA: backbone triplet, carbonyl, CB. Stored library frame == parent frame here,
# so an extra atom must land on its parent coordinate before the hit transform.
RES_NAMES = ["N", "CA", "C", "O", "CB"]
RES_ELEMS = ["N", "C", "C", "O", "C"]
RES_COORDS = np.array([[0.0, 0.0, 0.0], [1.46, 0.0, 0.0], [2.0, 1.4, 0.0],
                       [1.4, 2.4, 0.4], [2.0, -0.8, 1.2]])
CG_COORDS = np.array([[5.0, 0.0, 0.0], [6.0, 0.5, 0.0], [7.0, 0.0, 0.5]])

THETA = 0.4
R_HIT = np.array([[np.cos(THETA), -np.sin(THETA), 0.0],
                  [np.sin(THETA), np.cos(THETA), 0.0],
                  [0.0, 0.0, 1.0]])
T_HIT = np.array([2.0, -1.0, 3.0])


def _write_bucket(lib_dir, parent_pdb_dir):
    path = vdg_npz.vdg_npz_path(lib_dir, FRAG, 1, "ALA")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez_compressed(
        path,
        nr_cg_coords=CG_COORDS[None], cg_elements=np.array(["C", "C", "O"]),
        nr_cg_names=np.array([["C1", "C2", "O1"]]), nr_cg_seg=np.array([""]),
        nr_cg_chain=np.array(["L"]), nr_cg_resnum=np.array([1]),
        nr_cg_resname=np.array(["LIG"]),
        nr_vdm_bb_coords=RES_COORDS[:3][None][None].reshape(1, 1, 3, 3),
        nr_scrr_seg=np.array([[""]]), nr_scrr_chain=np.array([["A"]]),
        nr_scrr_resnum=np.array([[10]]), nr_scrr_resname=np.array([["ALA"]]),
        nr_parent_biounit=np.array([BIOUNIT]),
        parent_pdb_dir=np.array(parent_pdb_dir),
        aa_bucket_parts=np.array(["ALA"]),
        cluster_id=np.array([0]), cluster_size=np.array([3]))
    return path


def _write_mirror(root):
    sub = os.path.join(root, BIOUNIT[1:3].lower())
    os.makedirs(sub, exist_ok=True)
    ag = pr.AtomGroup("res")
    ag.setCoords(RES_COORDS)
    ag.setNames(np.array(RES_NAMES))
    ag.setResnames(np.array(["ALA"] * len(RES_NAMES)))
    ag.setResnums(np.array([10] * len(RES_NAMES)))
    ag.setChids(np.array(["A"] * len(RES_NAMES)))
    ag.setSegnames(np.array([""] * len(RES_NAMES)))
    ag.setElements(np.array(RES_ELEMS))
    ag.setOccupancies(np.ones(len(RES_NAMES)))
    pr.writePDB(os.path.join(sub, BIOUNIT + ".pdb"), ag)


def _write_tsv(path):
    cols = ["pdbfile", "frag", "bsr_combo", "subset_size", "aa_bucket", "vdg_index",
            "aa_perm_idx", "vdg_rmsd"]
    vals = ["query.pdb", FRAG, ":A:10", 1, "ALA", 0, 0, 0.42]
    for i in range(3):
        for j in range(3):
            cols.append(f"R{i}{j}")
            vals.append(R_HIT[i, j])
    for i in range(3):
        cols.append(f"t{i}")
        vals.append(T_HIT[i])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.write("\t".join(str(v) for v in vals) + "\n")
    return path


class HitPdbExtrasTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.lib = os.path.join(self.tmp.name, "frag_lib")
        self.mirror = os.path.join(self.tmp.name, "mirror")
        self.out = os.path.join(self.tmp.name, "out")
        self.tsv = _write_tsv(os.path.join(self.tmp.name, "hits.tsv"))
        # The mirror check is memoised process-wide; a stale entry would let a
        # deliberately-missing mirror pass.
        vdg_npz._checked_pdb_dirs.clear()
        vdg_npz._warned_missing_parent_db.clear()
        self.env = mock.patch.dict(os.environ, {}, clear=False)
        self.env.start()
        os.environ.pop(vdg_npz.PDB_DIR_ENV_VAR, None)
        self.addCleanup(self.env.stop)

    def _run(self, *extra):
        argv = ["write_vdg_hit_pdbs.py", "--hits-tsv", self.tsv,
                "--vdg-lib-dir", self.lib, "--outdir", self.out, *extra]
        with mock.patch.object(sys, "argv", argv), mock.patch("builtins.print") as p:
            whp.main()
        said = " ".join(str(c) for c in p.call_args_list)
        written = [os.path.join(d, f) for d, _, fs in os.walk(self.out) for f in fs]
        self.assertEqual(len(written), 1, written)
        return pr.parsePDB(written[0]), said

    def test_sidechain_rides_the_hit_transform(self):
        _write_bucket(self.lib, self.mirror)
        _write_mirror(self.mirror)
        ag, said = self._run()
        names = list(ag.getNames())
        self.assertEqual(names, ["C1", "C2", "O1", "N", "CA", "C", "CB"])
        np.testing.assert_allclose(ag.getCoords()[-1],
                                   RES_COORDS[4] @ R_HIT + T_HIT, atol=1e-2)
        self.assertNotIn("WARNING", said)

    def test_carbonyl_is_opt_in_and_lands_with_the_sidechain(self):
        _write_bucket(self.lib, self.mirror)
        _write_mirror(self.mirror)
        ag, _ = self._run("--carbonyl")
        self.assertEqual(list(ag.getNames()),
                         ["C1", "C2", "O1", "N", "CA", "C", "O", "CB"])
        np.testing.assert_allclose(ag.getCoords()[-2],
                                   RES_COORDS[3] @ R_HIT + T_HIT, atol=1e-2)

    def test_no_mirror_still_writes_backbone_only_hits(self):
        """The case most users are in: the library's recorded mirror is not here."""
        _write_bucket(self.lib, os.path.join(self.tmp.name, "not_a_mirror"))
        ag, said = self._run("--carbonyl")
        self.assertEqual(list(ag.getNames()), ["C1", "C2", "O1", "N", "CA", "C"])
        self.assertIn("--sidechain", said)
        self.assertIn("--carbonyl", said)

    def test_the_no_database_message_comes_from_the_shared_function(self):
        """Both PDB writers must print the same block; the way that is kept true is
        that neither owns the wording."""
        _write_bucket(self.lib, os.path.join(self.tmp.name, "not_a_mirror"))
        with mock.patch.object(vdg_npz, "missing_parent_db_message",
                               return_value="SHARED-BLOCK") as msg:
            _, said = self._run()
        self.assertIn("SHARED-BLOCK", said)
        self.assertEqual(msg.call_args.args[0], ["--sidechain"])

    def test_no_sidechain_writes_the_same_thing_without_reading_parents(self):
        _write_bucket(self.lib, self.mirror)
        _write_mirror(self.mirror)
        with mock.patch.object(vdg_npz, "parse_pdb_or_none",
                               side_effect=AssertionError("parent read")) as _:
            ag, said = self._run("--no-sidechain")
        self.assertEqual(list(ag.getNames()), ["C1", "C2", "O1", "N", "CA", "C"])
        self.assertNotIn("WARNING", said)

    def test_pdb_dir_without_a_consumer_is_rejected(self):
        _write_bucket(self.lib, self.mirror)
        with self.assertRaises(SystemExit):
            self._run("--no-sidechain", "--pdb-dir", self.mirror)


if __name__ == "__main__":
    unittest.main()
