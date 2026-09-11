"""The one place ligand chemistry is decided.

Every bond order, aromatic flag, formal charge and hydrogen count the pipeline
uses comes from here, so the CCD-template swap is a change to this module only.

Two entry points, for the two shapes of chemistry the pipeline needs:

* ``perceive_ligand_graph``   -- the ligand *type*: a coordinate-free RDKit Mol,
  used to enumerate fragments and write their keys.
* ``perceive_ligand_instance`` -- one observed *copy*: an OBMol from a PDB
  block, used to match those keys back onto a structure.

They must agree, or a key's ``H0``/``!H0`` (bond-order dependent, unlike the
graph-only ``D<n>``) can silently mine zero sites when perception misses a
bond. Both sides now read the same CCD template by atom name, so they agree by
construction (round-trip pinned by tests/test_ccd_templates.py).

Perception is the fallback; which one produced an observation is recorded per
observation so a read path can drop the perceived ones. Measured on 200
prepwizard structures (616 ligand residues): 100% templated, 41.4% of ligands
have >=1 unobserved heavy atom -- the case that makes perception read a
truncated degree.
"""
from collections import namedtuple

from openbabel import openbabel as ob
from rdkit import Chem

from ligand_vdgs.functions import ccd_templates

# Provenance of an observation's chemistry, so a read path can drop anything
# that did not come from the template. Non-negative, disjoint from the negative
# "unreadable" sentinel the int8 annotation columns use.
PERCEPTION_CCD_TEMPLATE = 0   # not reachable yet; task 4
PERCEPTION_OPENBABEL = 1      # coordinates + OB's bond-order perception
PERCEPTION_SMILES = 2         # the ligand table's SMILES string
PERCEPTION_ATOM_NAME_TABLE = 3  # hard-coded CG atom-name list (constants.cg_atoms);
                                # no ligand graph, so per-atom annotations are
                                # all the unreadable sentinel

PERCEPTION_LABELS = {PERCEPTION_CCD_TEMPLATE: 'ccd',
                     PERCEPTION_OPENBABEL: 'openbabel',
                     PERCEPTION_SMILES: 'smiles',
                     PERCEPTION_ATOM_NAME_TABLE: 'atom_name_table'}

PerceivedGraph = namedtuple('PerceivedGraph', ('mol', 'provenance'))
PerceivedInstance = namedtuple('PerceivedInstance', ('obmol', 'provenance'))

_CONVERSION = None


def _pdb_conversion():
    """One OBConversion per process; constructing it per ligand is measurable."""
    global _CONVERSION
    if _CONVERSION is None:
        _CONVERSION = ob.OBConversion()
        _CONVERSION.SetInFormat('pdb')
    return _CONVERSION


def perceive_ligand_graph(resname, smiles=None):
    """Full-chemistry Mol for a ligand *type*, or None if it cannot be built.

    CCD template first (keyed by `resname`), `smiles` only as fallback -- the
    template is also what the instance side uses, and the two must agree on
    bond orders or a key's `H0`/`!H0` won't match at mining time.

    Heavy atoms only, including atoms no structure ever observes: degree and H
    count come from the whole molecule. Skipping a fragment with an unobserved
    atom is the instance side's decision, not this one's.

    Deliberately unsanitized in both branches: for the SMILES fallback, None
    means a genuine syntax error, and element-filtering callers must see a
    metal complex that fails valence checks as non-druglike, not a parse failure.
    """
    template = _try_template(resname)
    if template is not None:
        mol = mol_from_template(template)
        if mol is not None:
            return PerceivedGraph(mol, PERCEPTION_CCD_TEMPLATE)
    if smiles is None:
        return PerceivedGraph(None, None)
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return PerceivedGraph(None, None)
    return PerceivedGraph(mol, PERCEPTION_SMILES)


def require_template_store():
    """Raise now if the CCD template store is unusable.

    `_try_template` degrades to perception on purpose, so a missing or stale
    store would otherwise produce a wholly OpenBabel-perceived library with no
    error anywhere -- exactly what the template swap exists to prevent,
    visible only afterwards in the provenance column. Job entry points call
    this once, before forking, so the failure is a crash at minute zero.
    """
    ccd_templates.connection()


def _try_template(resname):
    """The CCD template for `resname`, or None -- never raising.

    A missing or unbuildable template database must degrade to perception, not
    stop a mining job: `smarts_to_cgs.py` retries any exception 100 times with
    100-second sleeps.
    """
    if not resname:
        return None
    try:
        return ccd_templates.get_template(resname)
    except Exception:
        return None


def mol_from_template(template):
    """RDKit Mol of a template's heavy atoms, or None if it has none.

    Aromaticity comes from the CCD's own flags, not from RDKit's perception:
    pooling pyridine with pyrrole (and splitting furanose from pyranose) is a
    decision about the chemistry, and it should not change because RDKit's
    aromaticity model does. Not sanitized, for the same reason fragments are
    not: they inherit the parent's perception (docs/pitfalls.md).
    """
    heavy = [atom for atom in template.atoms if atom.element not in ('H', 'D')]
    if not heavy:
        return None
    h_counts = ccd_templates.template_h_counts(template)
    mol = Chem.RWMol()
    index = {}
    for atom in heavy:
        rd_atom = Chem.Atom(0)
        try:
            rd_atom.SetAtomicNum(Chem.GetPeriodicTable().GetAtomicNumber(
                atom.element))
        except Exception:
            # An element RDKit does not know (the CCD has a few): keep it as a
            # dummy rather than dropping the atom, so degrees stay right.
            rd_atom.SetAtomicNum(0)
        rd_atom.SetFormalCharge(int(atom.charge))
        rd_atom.SetIsAromatic(bool(atom.aromatic))
        # Explicit + noImplicit, so the count is the template's and RDKit never
        # adds its own from valence.
        rd_atom.SetNumExplicitHs(int(h_counts.get(atom.name, 0)))
        rd_atom.SetNoImplicit(True)
        rd_atom.SetProp('_ccdAtomName', atom.name)
        index[atom.name] = mol.AddAtom(rd_atom)
    for bond in template.bonds:
        i, j = index.get(bond.a), index.get(bond.b)
        if i is None or j is None:      # a bond to hydrogen
            continue
        if mol.GetBondBetweenAtoms(i, j) is not None:
            continue
        order = (Chem.BondType.AROMATIC if bond.aromatic
                 else _RDKIT_BOND_ORDER.get(bond.order, Chem.BondType.SINGLE))
        mol.AddBond(i, j, order)
        mol.GetBondBetweenAtoms(i, j).SetIsAromatic(bool(bond.aromatic))
    mol = mol.GetMol()
    mol.UpdatePropertyCache(strict=False)
    # Ring perception only -- ring_query_for_atom reads the bond graph, and
    # GetSymmSSSR is what makes `r<n>` mean *smallest* ring.
    Chem.GetSymmSSSR(mol)
    return mol


_RDKIT_BOND_ORDER = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
                     3: Chem.BondType.TRIPLE, 4: Chem.BondType.QUADRUPLE}


def perceive_ligand_instance(block, resname=None):
    """Hydrogen-free OBMol for one ligand's PDB block, or None if unreadable.

    Chemistry comes from the CCD template matched by atom name; OpenBabel's own
    perception is the fallback, recorded in `provenance` so a read path can
    drop the perceived observations.

    Hydrogen-free because both toolkits count explicit H in the SMARTS `D`
    primitive, so a protonated hydroxyl O reads `D2` and every degree-annotated
    key silently misses (ethanol/obabel 3: with H `[O;D2]`=1, `[C;D1]`=0; after
    DeleteHydrogens `[O;D1]`=1, `[C;D1]`=1).

    DeleteHydrogens runs after PerceiveBondOrders (which needs the H positions)
    and renumbers atoms, so every index a caller takes -- `pdb_atom_names`,
    SMARTS match tuples -- must come from the returned mol.
    """
    mol = ob.OBMol()
    if not _pdb_conversion().ReadString(mol, block) or mol.NumAtoms() == 0:
        return None
    mol.PerceiveBondOrders()
    mol.DeleteHydrogens()
    if mol.NumAtoms() == 0:
        return None
    template = _try_template(resname)
    if template is not None and apply_template_to_obmol(mol, template):
        return PerceivedInstance(mol, PERCEPTION_CCD_TEMPLATE)
    return PerceivedInstance(mol, PERCEPTION_OPENBABEL)


def apply_template_to_obmol(mol, template):
    """Overwrite `mol`'s chemistry from the template. True if it applied.

    False (mol untouched) only when an *observed* atom has no template name --
    a real mapping failure. A template atom that is not observed is normal
    (leaving atoms like OXT/HXT are absent under covalent attachment; partial
    density is common) and does not block the template.

    An unobserved heavy atom is added as a *phantom*: a real graph atom with no
    PDB atom name. Without it, `D<n>` would read one degree short next to
    unmodelled density and silently miss the fragment rather than skip it.
    With it, degree/H-count come from the template, and `find_cg_matches`
    already refuses any match that reaches an unnamed atom -- so a fragment
    with an unobserved atom is skipped, while one merely adjacent still mines.
    """
    by_name = ccd_templates.index_by_name(template)
    h_counts = ccd_templates.template_h_counts(template)
    atom_of_name, name_of_idx = {}, {}
    for i in range(1, mol.NumAtoms() + 1):
        atom = mol.GetAtom(i)
        residue = atom.GetResidue()
        name = residue.GetAtomID(atom).strip() if residue is not None else ''
        entry = by_name.get(name)
        if entry is None:
            return False
        name_of_idx[i] = entry.name
        # Two observed atoms mapping to one template atom is a duplicate-name
        # residue (docs: prepwizard re-emits orphan atoms under a ligand's
        # resname); the template cannot say which is which.
        if entry.name in atom_of_name:
            return False
        atom_of_name[entry.name] = i

    # Phantoms first, so every template heavy atom has an index before bonds
    # are laid down. AddAtom appends, so the observed atoms keep their indices
    # and every map built above stays valid.
    phantoms = set()
    for entry in template.atoms:
        if entry.element in ('H', 'D') or entry.name in atom_of_name:
            continue
        atom = mol.NewAtom()
        try:
            atom.SetAtomicNum(ob.GetAtomicNum(entry.element))
        except Exception:
            atom.SetAtomicNum(0)
        atom.SetVector(0.0, 0.0, 0.0)
        atom_of_name[entry.name] = atom.GetIdx()
        phantoms.add(atom.GetIdx())

    wanted = {}
    for bond in template.bonds:
        i, j = atom_of_name.get(bond.a), atom_of_name.get(bond.b)
        if i is None or j is None:      # to a hydrogen or an unobserved atom
            continue
        wanted[frozenset((i, j))] = bond

    # Drop what perception invented, then lay down the template's bonds. Done
    # in this order so a bond whose order is wrong is corrected rather than
    # duplicated.
    for bond in [mol.GetBond(k) for k in range(mol.NumBonds())]:
        pair = frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        if pair not in wanted:
            mol.DeleteBond(bond)
    for pair, spec in wanted.items():
        i, j = tuple(pair)
        existing = mol.GetBond(i, j)
        order = 5 if spec.aromatic else spec.order
        if existing is None:
            mol.AddBond(i, j, order)
            existing = mol.GetBond(i, j)
        else:
            existing.SetBondOrder(order)
        if spec.aromatic:
            existing.SetAromatic()
        else:
            existing.SetAromatic(False)

    for name, i in atom_of_name.items():
        atom = mol.GetAtom(i)
        entry = by_name[name]
        atom.SetFormalCharge(int(entry.charge))
        # The template's nominal H count, not the placed hydrogens: this is
        # what makes `H0`/`!H0` and cg_num_h independent of how the structure
        # was protonated.
        atom.SetImplicitHCount(int(h_counts.get(name, 0)))
        atom.SetAromatic(bool(entry.aromatic))

    # Otherwise OpenBabel re-perceives on the first Match and overwrites all of
    # the above with its own guess from the geometry.
    mol.SetAromaticPerceived()
    mol.SetHybridizationPerceived()
    return True


def phantom_atom_indices(obmol):
    """Indices of atoms `apply_template_to_obmol` added for unobserved atoms.

    Derived rather than stored: an atom is a phantom exactly when it carries no
    PDB atom name, which is the same test `find_cg_matches` uses to drop a
    match that reaches one.
    """
    named = pdb_atom_names(obmol)
    return {i for i in range(1, obmol.NumAtoms() + 1) if i not in named}


def pdb_atom_names(obmol):
    """{1-based OB atom index: PDB atom name} for the atoms that have one.

    Separate from perception because the matchers that only need presence
    (`count_matching_structures`, cost estimation) never build it. Read from
    the OBMol rather than from the nth HETATM line: obabel's atom order need
    not track the input lines, and a positional map fails silently when it
    does not. Atoms obabel perceived rather than read are absent.
    """
    names = {}
    for i in range(1, obmol.NumAtoms() + 1):
        atom = obmol.GetAtom(i)
        residue = atom.GetResidue()
        if residue is not None:
            names[i] = residue.GetAtomID(atom).strip()
    return names
