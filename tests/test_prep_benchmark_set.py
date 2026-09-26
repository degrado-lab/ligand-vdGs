"""prep_benchmark_set altloc handling: highest-occupancy pick, near-ties, no dropped protein
sites, and the ProDy default-parse falsifier the rule exists to fix."""
import prody as pr
from ligand_vdgs.tools import prep_benchmark_set as p

# ALA(1) no altloc; SER(2) A=0.6/C=0.4; GLY(3) B=0.7/C=0.3 (modelled only in B/C, not A);
# LIG(101) label B only -> highest summed occupancy is unambiguous.
PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00 10.00           C
ATOM      3  N  ASER A   2       2.000   0.000   0.000  0.60 10.00           N
ATOM      4  N  CSER A   2       2.100   0.000   0.000  0.40 10.00           N
ATOM      5  CA ASER A   2       3.000   0.000   0.000  0.60 10.00           C
ATOM      6  CA CSER A   2       3.100   0.000   0.000  0.40 10.00           C
ATOM      7  N  BGLY A   3       4.000   0.000   0.000  0.70 10.00           N
ATOM      8  N  CGLY A   3       4.100   0.000   0.000  0.30 10.00           N
ATOM      9  CA BGLY A   3       4.500   0.000   0.000  0.70 10.00           C
ATOM     10  CA CGLY A   3       4.600   0.000   1.000  0.30 10.00           C
HETATM   11  C1 BLIG A 101       5.000   0.000   0.000  0.55 10.00           C
HETATM   12  C2 BLIG A 101       6.000   0.000   0.000  0.55 10.00           C
END
"""

def _fixture(tmp_path, pdb=PDB):
    path = tmp_path / "alt.pdb"
    path.write_text(pdb)
    return p.parse_structure(str(path)), str(path)

def test_altloc_units_picks_highest_occupancy_and_flags_near_ties(tmp_path):
    ag, _ = _fixture(tmp_path)
    occ, units = p.altloc_units(ag, p.ligand_atoms(ag, "LIG"))
    assert occ == {"B": 1.1}
    assert units == [("B", False)]

def test_altloc_units_near_tie_yields_one_unit_per_tied_label(tmp_path):
    tied_pdb = PDB[:-4] + (
        "HETATM   13  C1 ALIG A 102       5.000   0.000   0.000  0.52 10.00           C\n"
        "HETATM   14  C2 ALIG A 102       6.000   0.000   0.000  0.52 10.00           C\nEND\n")
    ag, _ = _fixture(tmp_path, tied_pdb)
    occ, units = p.altloc_units(ag, p.ligand_atoms(ag, "LIG"), tie=0.2)
    assert occ == {"B": 1.1, "A": 1.04}
    assert dict(units) == {"B": True, "A": True}

def test_resolve_altlocs_never_drops_a_protein_site(tmp_path):
    ag, _ = _fixture(tmp_path)
    keep = p.resolve_altlocs(ag, "B")
    kept = sorted(zip(ag.getResnums()[keep], ag.getNames()[keep], ag.getAltlocs()[keep]))
    assert kept == [(1, "CA", " "), (1, "N", " "), (2, "CA", "A"), (2, "N", "A"),
                    (3, "CA", "B"), (3, "N", "B"), (101, "C1", "B"), (101, "C2", "B")]
    assert p.protein_sites(ag[keep].select("protein")) == p.protein_sites(ag)

def test_resolve_altlocs_falsifier_default_prody_parse_drops_a_residue(tmp_path):
    # This is the bug the rule exists to fix: ProDy's default (single-label) parse takes
    # altloc A, which GLY(3) never carries, so the whole residue vanishes from `protein`.
    ag, path = _fixture(tmp_path)
    default = pr.parsePDB(path, altloc="A")
    assert p.protein_sites(default) != p.protein_sites(ag)
    assert len(p.protein_sites(ag)) - len(p.protein_sites(default)) == 2  # GLY N, CA lost

def test_bsr_altloc_audit_flags_divergent_backbone(tmp_path):
    ag, _ = _fixture(tmp_path)
    lig_idx = p.ligand_atoms(ag, "LIG")
    n_bsr, n_alt, divergent = p.bsr_altloc_audit(ag, lig_idx, cutoff=5.0, bb_tol=0.5)
    assert n_bsr == 3  # ALA(1) at exactly 5.0A, SER(2), GLY(3) all within cutoff
    assert n_alt == 2  # only SER, GLY carry altlocs
    assert divergent == ["A:3"]  # GLY's CA moves 1.0A between B/C; SER's does not move at all
