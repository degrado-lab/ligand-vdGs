# utils.py

import itertools
import os
import re
import shutil
from rdkit import Chem
from rdkit.Chem import rdMolAlign as MA
import numpy as np

import hashlib


def file_sha256(path):
    """SHA-256 of a file, streamed. Used to pin derived artifacts (the cost
    estimate TSV, the library provenance record) to the fragment dict they were
    computed from."""
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def _int_or_none(v):
    if v.lower() == 'none':
        return None
    return int(v)

def set_up_outdir(outdir, overwrite=False):
    '''Create outdir if it does not exist, and require overwrite if it does, so that there are '
    no stale files.'''
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise ValueError(f'[ERROR] The filename you designated as the output directory, {outdir}, '
                             'already exists and is not a directory.')
        if overwrite:
            print(f'[WARNING] Overwriting existing output directory {outdir} because '
                  'overwrite_existing was set to True.')
            # Another process might remove/create concurrently; guard with try/except.
            try:
                shutil.rmtree(outdir)
            except FileNotFoundError:
                pass
            os.makedirs(outdir, exist_ok=True)
        else:
            # Allow empty dirs; error only if files are present.
            with os.scandir(outdir) as entries:
                not_empty = any(entries)
            if not_empty:
                raise ValueError(f'[ERROR] The output directory {outdir} is not empty. Remove files or '
                                 'set overwrite_existing to True to prevent accidental overwriting.')
    else:
        parent_dir = os.path.dirname(outdir)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)
        os.makedirs(outdir, exist_ok=True)

def valid_database_subdir_format(input_dir):
    ''' Ensure that the pdb database dir has subdirs formatted similarly to the RCSB
    mirror format (see docs/database_generation_guide.md). '''
    
    def error_message():
        print('Database structure must be similar to the RCSB mirror format; '
              'see docs/database_generation_guide.md')
    
    # Input path must be a dir
    if not os.path.exists(input_dir):
        error_message()
        return False
    # Input dir must not be empty
    if len(os.listdir(input_dir)) == 0:
        error_message()
        return False
    # Ensure that there is at least one correctly formatted subdir
    for p in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, p)
        if os.path.isdir(subdir_path):
            # Subdirs should be the inner 2 characters of a 4-char pdb file name
            if len(p) == 2 and p == p.lower():
                for pdb in os.listdir(subdir_path):
                    if len(pdb) >= 3 and pdb[1:3].lower() == p:
                        return True
    # No valid subdir was found
    error_message()
    return False

def smiles_equiv(existingfrag, sub_smiles):
    """Whether two fragment keys denote the same substructure.

    Delegates to fragment_keys_equivalent: comparing ring-annotated keys
    directly is not even reflexive, because `r<n>` cannot be satisfied on the
    acyclic graph a fragment key parses to.
    """
    return fragment_keys_equivalent(existingfrag, sub_smiles)

# Elements whose terminal atoms are treated as interchangeable positions.
_TERMINAL_RESONANCE_ELEMENTS = frozenset({7, 8, 16})  # N, O, S
# Centers those terminal atoms may hang off: B, C, N, O, P, S, Cl, Br, I.
# The halogens cover perchlorate/periodate-style oxyanions; B covers boronates.
_TERMINAL_RESONANCE_CENTERS = frozenset({5, 6, 7, 8, 15, 16, 17, 35, 53})


def _has_deliberate_charge_assignment(mol, bond_indices):
    """True if a terminal set's drawn charges should be trusted as-is.

    See the carve-out described in ``_find_resonance_terminal_groups``: four or
    more terminal atoms, every bond single, and charges that are not uniform.
    """
    if len(bond_indices) < 4:
        return False
    if any(mol.GetBondWithIdx(i).GetBondTypeAsDouble() != 1.0 for i in bond_indices):
        return False

    charges = set()
    for bond_idx in bond_indices:
        bond = mol.GetBondWithIdx(bond_idx)
        for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
            if atom.GetDegree() == 1:
                charges.add(atom.GetFormalCharge())
    return len(charges) > 1


def _find_resonance_terminal_groups(mol):
    """Find terminal X-N/O/S groups whose drawn bond order and charge may move.

    Two or more terminal atoms of the same element on a common center are
    treated as one set of interchangeable positions. In these hydrogen-free
    fragments a drawn ``-OH``, ``=O``, and ``[O-]`` on the same center are the
    same physical position recorded in different protonation or resonance
    states, so the drawing must not decide which vdG atom maps to which.

    Each element forms its own set, so a terminal O is never exchanged with a
    terminal S (see ``OP(O)(=S)[S-]``: the O pair and the S pair each permute,
    but not across). Aromatic terminal atoms are excluded -- a truncated ring
    atom is a real ring position, not a resonance form of an exocyclic
    substituent, so the ring ``o`` and the hydroxyl ``O`` of ``cc(o)O`` stay
    distinct. Degree, not implicit H count, defines "terminal" here because
    these fragments carry no hydrogens.

    Substituted and bridging atoms represented in the fragment are untouched.

    One carve-out: a saturated center carrying four or more terminal atoms of an
    element, all single-bonded but not all the same charge, is left alone. Such
    a drawing has no bond-order ambiguity left to resolve, so its explicit
    charge assignment is taken as deliberate. This keeps the bisulfate drawing
    ``[O-]S([O-])([O-])O`` at three interchangeable anionic oxygens plus a fixed
    hydroxyl, rather than merging all four. It fires on that fragment alone in
    the current dictionary; three-terminal centers such as ``cS([O-])(O)O`` are
    unaffected and still merge fully.

    Returns: dict {bond_idx: (center_atomic_num, group_id)}
    """
    resonance_groups = {}
    next_group_id = 0

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _TERMINAL_RESONANCE_CENTERS:
            continue

        terminal_bonds_by_element = {}
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            atomic_num = other.GetAtomicNum()
            if atomic_num not in _TERMINAL_RESONANCE_ELEMENTS:
                continue
            if other.GetIsAromatic() or other.GetDegree() != 1:
                continue
            terminal_bonds_by_element.setdefault(atomic_num, []).append(bond.GetIdx())

        for atomic_num, bond_indices in sorted(terminal_bonds_by_element.items()):
            if len(bond_indices) < 2:
                continue
            if _has_deliberate_charge_assignment(mol, bond_indices):
                continue
            for bond_idx in bond_indices:
                resonance_groups[bond_idx] = (atomic_num, next_group_id)
            next_group_id += 1

    return resonance_groups


def _find_resonance_center_groups(mol, center_z, neighbor_z, charge_rule):
    """Shared body of the C-N / N-O / S-N finders below.

    Recognizes a non-aromatic ``center_z`` atom carrying >=2 ``neighbor_z``
    neighbors, accepted when it has a double bond to one of them or when
    ``charge_rule`` (``None``, ``'positive'`` or ``'nonzero'``) admits its formal
    charge. Group ids restart at 0 per call; the tag's ``center_z`` is what keeps
    the finders' outputs distinct when the caller merges them.
    """
    if charge_rule not in (None, 'positive', 'nonzero'):
        # Unrecognized rules would otherwise fall through both tests below and
        # silently accept every center.
        raise ValueError(f"unknown charge_rule {charge_rule!r}")

    resonance_groups = {}
    next_group_id = 0

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != center_z or atom.GetIsAromatic():
            continue

        neighbor_bonds = []
        has_double_bond = False
        for bond in atom.GetBonds():
            if bond.GetOtherAtom(atom).GetAtomicNum() != neighbor_z:
                continue
            neighbor_bonds.append(bond.GetIdx())
            if bond.GetBondTypeAsDouble() == 2.0:
                has_double_bond = True

        if len(neighbor_bonds) < 2:
            continue
        if not has_double_bond:
            charge = atom.GetFormalCharge()
            if charge_rule is None:
                continue
            if charge_rule == 'positive' and charge <= 0:
                continue
            if charge_rule == 'nonzero' and charge == 0:
                continue

        for bond_idx in neighbor_bonds:
            resonance_groups[bond_idx] = (center_z, next_group_id)
        next_group_id += 1

    return resonance_groups


def _find_resonance_CN_groups(mol):
    """Find amidine/guanidine C-N groups whose drawn form is ignored.

    A group is recognized at a non-aromatic carbon with at least two nitrogen
    neighbors and either a C=N double bond or positive formal charge on the
    carbon. The latter admits all-single-bond ``[C+]`` resonance drawings. All
    C-N bonds at a recognized center are normalized; exact graph matching still
    preserves every substituent attached beyond each nitrogen.

    Ordinary aminals and nitrile-like C#N groups are deliberately excluded.
    They have neither the amidine/guanidine C=N bond nor the charged-carbon
    resonance representation.

    Returns: dict {bond_idx: (center_atomic_num, group_id)}
    """
    return _find_resonance_center_groups(mol, 6, 7, 'positive')


def _find_resonance_NO_groups(mol):
    """Find N-O groups whose bond/charge drawing may move among oxygens.

    The rule covers a non-aromatic nitrogen with at least two oxygen neighbors
    when either an N=O bond or a formal charge on nitrogen marks a resonance
    representation. It therefore includes both ``cN(=O)O`` and
    ``C=[N+]([O-])O``. Exact graph matching still preserves every represented
    substituent beyond each oxygen.

    Returns: dict {bond_idx: (center_atomic_num, group_id)}
    """
    return _find_resonance_center_groups(mol, 7, 8, 'nonzero')


def _find_resonance_SN_groups(mol):
    """Find S-N groups whose single/double-bond placement may exchange.

    A group is recognized at a non-aromatic sulfur with at least two nitrogen
    neighbors and at least one S=N double bond. This covers S(VI)-N resonance
    drawings such as ``CS(=N)(N)=O`` and ``N=S(N)(=O)F`` without treating
    ordinary all-single-bond S-N groups as equivalent.

    Returns: dict {bond_idx: (center_atomic_num, group_id)}
    """
    return _find_resonance_center_groups(mol, 16, 7, None)


def _find_resonance_NN_groups(mol):
    """Find conjugated non-aromatic N chains whose bond placement may move.

    Components must contain at least three nitrogens joined by N-N bonds and at
    least one N=N double bond. Normalizing only bonds inside that component lets
    exact graph matching recover the reversal of ``cN=NNc`` while retaining all
    represented attachments at the two ends.

    Returns: dict {bond_idx: (center_atomic_num, group_id)}
    """
    nitrogen_indices = {
        atom.GetIdx() for atom in mol.GetAtoms()
        if atom.GetAtomicNum() == 7 and not atom.GetIsAromatic()
    }
    visited = set()
    resonance_groups = {}
    next_group_id = 0

    for start in sorted(nitrogen_indices):
        if start in visited:
            continue

        component = set()
        component_bonds = set()
        stack = [start]
        visited.add(start)
        while stack:
            atom_idx = stack.pop()
            component.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for bond in atom.GetBonds():
                neighbor_idx = bond.GetOtherAtomIdx(atom_idx)
                if neighbor_idx not in nitrogen_indices:
                    continue
                component_bonds.add(bond.GetIdx())
                if neighbor_idx not in visited:
                    visited.add(neighbor_idx)
                    stack.append(neighbor_idx)

        has_double_bond = any(
            mol.GetBondWithIdx(bond_idx).GetBondTypeAsDouble() == 2.0
            for bond_idx in component_bonds
        )
        if len(component) < 3 or not has_double_bond:
            continue

        for bond_idx in component_bonds:
            resonance_groups[bond_idx] = (7, next_group_id)
        next_group_id += 1

    return resonance_groups


def _find_resonance_aromatic_N_atoms(mol):
    """Find aromatic N atoms whose charge/proton placement may be ignored.

    Only aromatic nitrogens in an aromatic-bond-connected component containing
    at least two aromatic nitrogens participate. Exact graph matching still
    decides whether any of those atoms can exchange. Requiring a shared aromatic
    component prevents unrelated sites such as the two ends of [n+]CCCn from
    becoming equivalent merely because their formal charges differ.
    """
    aromatic_atom_indices = {
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()
    }
    visited = set()
    resonance_atoms = set()

    for start in sorted(aromatic_atom_indices):
        if start in visited:
            continue

        component = set()
        stack = [start]
        visited.add(start)
        while stack:
            atom_idx = stack.pop()
            component.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for bond in atom.GetBonds():
                if not bond.GetIsAromatic():
                    continue
                neighbor_idx = bond.GetOtherAtomIdx(atom_idx)
                if neighbor_idx not in aromatic_atom_indices:
                    continue
                if neighbor_idx in visited:
                    continue
                visited.add(neighbor_idx)
                stack.append(neighbor_idx)

        nitrogen_atoms = {
            atom_idx for atom_idx in component
            if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7
        }
        if len(nitrogen_atoms) >= 2:
            resonance_atoms.update(nitrogen_atoms)

    return resonance_atoms


def mol_from_fragment(fragment):
    """Parse the pipeline's hydrogen-free fragment notation as unsanitized SMARTS.

    These strings describe extracted substructures, not complete molecules. SMARTS
    parsing preserves their query atoms and does not infer hydrogens, run valence
    checks, or sanitize them as though they were standalone SMILES molecules.
    """
    return init_query_ring_info(Chem.MolFromSmarts(fragment))


# Ring primitives encoded as isotopes so that mutual substructure matching can
# see WHICH atom carries WHICH annotation. `!R`/`R`/`r<n>` are unsatisfiable on
# the acyclic query graphs fragment keys parse to, so a direct match of two
# annotated keys is always False; an isotope is an ordinary matchable atom
# property and preserves the atom-to-annotation correspondence that a bare
# skeleton comparison throws away.
_RING_PRIMITIVE = re.compile(r'\[([^\];]*);(!R|R|r\d+)\]')
_RING_ISOTOPE = {'!R': 1, 'R': 2}


def _ring_annotation_as_isotope(fragment):
    """A fragment key rewritten with its ring primitives as isotope labels."""
    def _encode(match):
        inner, annotation = match.group(1), match.group(2)
        isotope = _RING_ISOTOPE.get(annotation)
        if isotope is None:
            isotope = 10 + int(annotation[1:])
        return f'[{isotope}{inner}]'
    return _RING_PRIMITIVE.sub(_encode, fragment)


def fragment_keys_equivalent(key_a, key_b):
    """Whether two fragment keys describe the same annotated substructure.

    Canonical SMILES alone does not settle this: fragments inherit the parent's
    unsanitized perception flags (docs/pitfalls.md), so two isomorphic submols
    can canonicalize to different strings -- `[C;!R][C;!R][O;!R][C;!R]` and
    `[C;!R][O;!R][C;!R][C;!R]` are one fragment written two ways.

    Matching the key against the *submol* it came from cannot answer it either:
    the submol has no RingInfo (PathToSubmol output is never sanitized), so any
    ring primitive raises, and PathToSubmol has cut the rings anyway, so an
    `r5` key would not match its own fragment even with RingInfo supplied. Both
    keys are therefore compared to each other, with the annotation carried as an
    isotope so the correspondence survives.
    """
    return fragment_query_mols_equivalent(fragment_key_query_mol(key_a),
                                          fragment_key_query_mol(key_b))


def fragment_key_query_mol(key):
    """The isotope-encoded, ring-info'd query mol that key comparison consumes.

    Callers comparing one key against many (the O(bucket^2) scan in
    fragment_database_ligs) should cache this per key and use
    ``fragment_query_mols_equivalent`` rather than re-parsing both sides per pair.
    """
    return init_query_ring_info(Chem.MolFromSmarts(_ring_annotation_as_isotope(key)))


def fragment_query_mols_equivalent(mol_a, mol_b):
    """Whether two mols from ``fragment_key_query_mol`` are the same substructure."""
    if mol_a is None or mol_b is None:
        return False
    if mol_a.GetNumAtoms() != mol_b.GetNumAtoms():
        return False
    if Chem.GetFormalCharge(mol_a) != Chem.GetFormalCharge(mol_b):
        return False
    return bool(mol_a.HasSubstructMatch(mol_b) and mol_b.HasSubstructMatch(mol_a))


def init_query_ring_info(mol):
    """Give a SMARTS query mol the RingInfo that ring primitives need.

    ``MolFromSmarts`` does not perceive rings, so a query mol used as the
    *target* of a match raises ``RingInfo not initialized`` as soon as either
    side carries `R` or `r<n>` -- which every fragment key does now (see
    Frags.ring_query_for_atom). ``FastFindRings`` supplies it without
    sanitizing, which matters because these are query mols with deliberately
    incomplete valences.
    """
    if mol is None:
        return mol
    # Called unconditionally: RingInfo exposes no portable "is initialized"
    # predicate across RDKit versions, and FastFindRings is idempotent and cheap
    # on a 4-5 atom query.
    Chem.FastFindRings(mol)
    return mol


# Bracket atom whose symbol is an aromatic (lower-case) element and which
# carries an explicit total-H count: [nH], [nH1], [cH], ...
_AROMATIC_H_ATOM = re.compile(r'\[([a-z][a-z]?)(H\d?)([^\]]*)\]')


def aromatic_h_constrained_atoms(fragment):
    """Bracket atoms in a fragment SMARTS that pin an H count on an aromatic atom.

    ``identify_mol_automorphisms`` normalizes H (and charge) among the N atoms of
    one aromatic component on purpose, so it treats an azole's tautomers as the
    same graph. The *matching* layer does not: an unbracketed aromatic atom in
    SMARTS carries no H constraint, but ``[nH]`` demands exactly one, and the
    structures being mined are protonated by prepwizard rather than by the
    depositor. A query that names the tautomer therefore silently keeps only the
    ligands modeled in it -- measured over the mirror, ``c1nnn[nH]1`` finds 7% of
    the tetrazoles ``cnnnn`` does and ``c1cnc[nH]1`` 15% of the imidazoles.

    Fragment keys written by ``fragment_database_ligs.py`` never carry an H count,
    so this only fires on a hand-written ``-s``. Returns the offending bracket
    texts, empty when there are none.
    """
    return [match.group(0) for match in _AROMATIC_H_ATOM.finditer(fragment or '')]


def _require_hydrogen_free(mol, context):
    """Reject explicit hydrogen atom *nodes*, in a molecule or a SMARTS query.

    The test is atomic number, which covers plain H atoms and the SMARTS tokens
    that pin one (`[H]`, `[#1]`, `[2H]` all report GetAtomicNum() == 1). Query
    atoms that merely *allow* hydrogen report 0 and slip through (`[#1,#6]`,
    `[$([H])]`, `[*H]`); the pipeline's fragment SMARTS never contain those. An
    H *count* on a heavy atom (`[nH]`, `[CH3]`) is a property, not a node, and
    is deliberately allowed.
    """
    if any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms()):
        raise ValueError(
            f"{context} requires a hydrogen-free molecular graph; strip explicit "
            "hydrogen atom nodes upstream"
        )


def _automorphism_graph(mol):
    """Return atom labels and a labeled adjacency matrix for symmetry matching.

    The graph preserves the molecular atom and bond attributes used for atom
    correspondence, except that bonds and incident-atom charges are normalized
    in recognized resonance groups: sets of terminal N/O/S atoms sharing a
    center, plus C-N, N-O, S-N, and conjugated N-N groups. Formal charge and
    declared-H distinctions are also normalized for aromatic nitrogens in the
    same aromatic component. Other SMARTS query labels are preserved outside
    recognized groups. This makes alternate drawings of those groups equivalent
    without treating arbitrary same-element atoms as interchangeable.
    """
    num_atoms = mol.GetNumAtoms()
    resonance_bonds = _find_resonance_terminal_groups(mol)
    resonance_bonds.update(_find_resonance_CN_groups(mol))
    resonance_bonds.update(_find_resonance_NO_groups(mol))
    resonance_bonds.update(_find_resonance_SN_groups(mol))
    resonance_bonds.update(_find_resonance_NN_groups(mol))
    aromatic_resonance_atoms = _find_resonance_aromatic_N_atoms(mol)
    resonance_atoms = set(aromatic_resonance_atoms)
    for bond_idx in resonance_bonds:
        bond = mol.GetBondWithIdx(bond_idx)
        resonance_atoms.add(bond.GetBeginAtomIdx())
        resonance_atoms.add(bond.GetEndAtomIdx())

    atom_labels = []
    for atom in mol.GetAtoms():
        # Query atoms from MolFromSmarts do not necessarily have a computed
        # implicit valence, so only properties that are safe on both QueryAtom
        # and Atom are used here.
        atom_idx = atom.GetIdx()
        if atom_idx in aromatic_resonance_atoms:
            query_label = ("aromatic_N_resonance",)
        elif atom.HasQuery() and atom_idx not in resonance_atoms:
            query_label = atom.GetSmarts()
        else:
            query_label = None
        atom_labels.append((
            atom.GetAtomicNum(),
            atom.GetIsAromatic(),
            atom.GetIsotope(),
            0 if atom.GetIdx() in resonance_atoms else atom.GetFormalCharge(),
            atom.GetNumRadicalElectrons(),
            int(atom.GetChiralTag()),
            query_label,
        ))

    adjacency = [[None] * num_atoms for _ in range(num_atoms)]
    for bond in mol.GetBonds():
        if bond.GetIdx() in resonance_bonds:
            # The shared label is deliberate: group membership is already
            # encoded by incidence on the same center atom. A per-center ID
            # would incorrectly prevent automorphisms that exchange two
            # otherwise equivalent resonance centers.
            bond_label = (1,)
        else:
            bond_label = (
                0,
                float(bond.GetBondTypeAsDouble()),
                bond.GetIsAromatic(),
                bond.GetIsConjugated(),
                int(bond.GetStereo()),
                int(bond.GetBondDir()),
                bond.GetSmarts() if bond.HasQuery() else None,
            )
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        adjacency[i][j] = bond_label
        adjacency[j][i] = bond_label

    return atom_labels, adjacency


def _wl_colors(atom_labels, adjacency):
    """WL refinement used only to prune exact automorphism enumeration."""
    num_atoms = len(atom_labels)
    colors = list(atom_labels)

    for _ in range(max(1, num_atoms)):
        signatures = []
        for i in range(num_atoms):
            neighbors = [
                (adjacency[i][j], colors[j])
                for j in range(num_atoms)
                if adjacency[i][j] is not None
            ]
            neighbors.sort(key=repr)
            signatures.append((colors[i], tuple(neighbors)))

        signature_ids = {}
        refined = []
        for signature in signatures:
            if signature not in signature_ids:
                signature_ids[signature] = len(signature_ids)
            refined.append(signature_ids[signature])

        # Including the previous color means refinement only splits classes.
        # Equal class counts therefore imply that the partition is stable.
        if len(set(refined)) == len(set(colors)):
            return refined
        colors = refined

    return colors


def validate_atom_permutations(permutations, num_atoms=None):
    """Validate and normalize index permutations to a tuple of unique tuples."""
    if permutations is None:
        return None

    normalized = []
    seen = set()
    inferred_num_atoms = num_atoms
    for permutation in permutations:
        try:
            perm = tuple(int(i) for i in permutation)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid atom permutation: {permutation!r}") from exc
        if inferred_num_atoms is None:
            inferred_num_atoms = len(perm)
        if len(perm) != inferred_num_atoms:
            raise ValueError(
                f"Atom permutation has {len(perm)} entries; expected {inferred_num_atoms}"
            )
        if set(perm) != set(range(inferred_num_atoms)):
            raise ValueError(f"Not a permutation of 0..{inferred_num_atoms - 1}: {perm}")
        if perm not in seen:
            seen.add(perm)
            normalized.append(perm)

    if not normalized:
        raise ValueError("At least one atom permutation is required")

    identity = tuple(range(inferred_num_atoms))
    if identity not in seen:
        raise ValueError("Atom permutations must include the identity mapping")

    normalized.sort()
    normalized.remove(identity)
    normalized.insert(0, identity)
    return tuple(normalized)


def group_preserving_permutations(labels):
    """Return all index permutations that only swap positions sharing a label.

    Shared by ``vdg_fp_utils.slot_orders`` (generation) and
    ``vdg_npz_utils.aa_perm_indices`` (read path) -- both need the same
    slot-permutation quotient, and keeping one implementation is what makes
    that guaranteed rather than coincidental.

    E.g.  ['ALA', 'ALA', 'SER'] -> [[0,1,2], [1,0,2]]
          ['ASP', 'HIS']        -> [[0,1]]
          ['bb', 'bb', 'bb']    -> all 6 permutations of [0,1,2]
    """
    labels = list(labels)
    n = len(labels)

    groups = {}
    for i, label in enumerate(labels):
        groups.setdefault(label, []).append(i)
    groups_list = list(groups.values())
    group_iperms = [list(itertools.permutations(g)) for g in groups_list]

    result = []
    for combo in itertools.product(*group_iperms):
        perm = [None] * n
        for pos_list, perm_of_pos in zip(groups_list, combo):
            for pos, orig in zip(pos_list, perm_of_pos):
                perm[pos] = orig
        result.append(perm)
    return result


def identify_mol_automorphisms(mol, max_automorphisms=10000):
    """Enumerate exact automorphisms of a protonation-normalized molecular graph.

    WL colors are used as candidate partitions only. Every returned mapping is
    then checked by a full labeled adjacency-preservation test, so equal WL
    colors can never introduce a non-automorphic atom permutation.

    Explicit hydrogen atom nodes are rejected. SMARTS-declared H constraints,
    bond order, and formal charge are preserved normally but normalized within
    recognized resonance groups: terminal N/O/S sets on a shared center, and
    C-N, N-O, S-N, and conjugated N-N groups. Charge/H distinctions are also
    normalized among N atoms in one aromatic component.

    Each returned tuple is the mapping itself: ``perm[i]`` is the atom that
    atom ``i`` is sent to. That is the inverse of NumPy's
    ``coords[permutation]`` gather order. The returned set is a group, hence
    closed under inversion, so a consumer that enumerates all of it -- which
    every consumer must, see docs/symmetry_edge_cases.md -- sees the same set
    of reorderings either way; only code applying a single tuple on its own
    needs the direction.
    """
    if mol is None:
        raise ValueError("Cannot identify automorphisms of a null molecule")
    _require_hydrogen_free(mol, "Automorphism enumeration")
    if not isinstance(max_automorphisms, int) or isinstance(max_automorphisms, bool) \
            or max_automorphisms < 1:
        raise ValueError("max_automorphisms must be a positive integer")

    atom_labels, adjacency = _automorphism_graph(mol)
    num_atoms = len(atom_labels)
    if num_atoms == 0:
        return ((),)

    colors = _wl_colors(atom_labels, adjacency)
    color_members = {}
    for atom_idx, color in enumerate(colors):
        color_members.setdefault(color, []).append(atom_idx)

    # Assign the most constrained vertices first. Trying the same index first
    # finds the identity immediately and makes low caps deterministic.
    order = sorted(
        range(num_atoms),
        key=lambda i: (len(color_members[colors[i]]),
                       -sum(edge is not None for edge in adjacency[i]), i),
    )
    mapping = [-1] * num_atoms
    used_targets = [False] * num_atoms
    automorphisms = []

    def search(depth):
        if len(automorphisms) > max_automorphisms:
            return
        if depth == num_atoms:
            automorphisms.append(tuple(mapping))
            return

        source = order[depth]
        candidates = sorted(
            color_members[colors[source]], key=lambda target: target != source
        )
        for target in candidates:
            if used_targets[target]:
                continue

            preserves_graph = True
            for other_source, other_target in enumerate(mapping):
                if other_target < 0:
                    continue
                if adjacency[source][other_source] != adjacency[target][other_target]:
                    preserves_graph = False
                    break
            if not preserves_graph:
                continue

            mapping[source] = target
            used_targets[target] = True
            search(depth + 1)
            used_targets[target] = False
            mapping[source] = -1
            if len(automorphisms) > max_automorphisms:
                return

    search(0)
    if len(automorphisms) > max_automorphisms:
        raise ValueError(
            f"Molecule has more than max_automorphisms={max_automorphisms}; "
            "refusing to use an incomplete symmetry mapping set"
        )
    return validate_atom_permutations(automorphisms, num_atoms)


def _proper_kabsch_rotations(H):
    """Return proper rotations for a batch of 3x3 cross-covariance matrices."""
    if not np.isfinite(H).all():
        raise ValueError("Non-finite Kabsch cross-covariance matrix")

    U, _, Vt = np.linalg.svd(H, full_matrices=False)
    UVt = np.matmul(U, Vt)
    det_uvt = _det3(UVt)   # same hand-expanded 3x3 determinant kabsch_ssd uses
    if not np.isfinite(det_uvt).all():
        raise ValueError("Non-finite determinant from Kabsch SVD")

    # A singular H has zero singular values, but U and Vt remain orthogonal, so
    # det(U @ Vt) is still +/-1. The comparison form defensively keeps D proper
    # even if a numerical backend ever reports signed zero.
    d = np.where(det_uvt < 0.0, -1.0, 1.0).astype(H.dtype, copy=False)
    D = np.zeros_like(H)
    D[..., 0, 0] = 1.0
    D[..., 1, 1] = 1.0
    D[..., 2, 2] = d
    return np.matmul(U, np.matmul(D, Vt))


def _det3(H):
    """Determinant of a batch of 3x3 matrices.

    ``np.linalg.det`` goes through LAPACK and the ``__array_function__``
    protocol; on 3x3 inputs called hundreds of thousands of times that dominates
    the arithmetic. Only the sign is used here, and it is exact at the only
    place it matters -- a singular H pairs it with a zero singular value.
    """
    a, b, c = H[..., 0, 0], H[..., 0, 1], H[..., 0, 2]
    d, e, f = H[..., 1, 0], H[..., 1, 1], H[..., 1, 2]
    g, h, i = H[..., 2, 0], H[..., 2, 1], H[..., 2, 2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def kabsch_ssd(X, Y, chunk_size=30000):
    """Squared deviations for the optimal superposition, without building R or t.

    Same inputs, same convention and the same ``ssd`` as ``kabsch``, but every
    caller that only thresholds an RMSD can skip most of the work: the optimal
    SSD is available from the singular *values* of the cross-covariance alone,

        ssd = |Xc|^2 + |Yc|^2 - 2 * (s1 + s2 + d * s3),   d = sign(det H)

    so ``compute_uv=False`` suffices and the rotation, translation and explicit
    residual all disappear. ``d`` enforces a proper rotation exactly as
    ``_proper_kabsch_rotations`` does: ``sign(det(U @ Vt)) == sign(det(H))``
    whenever ``H`` is nonsingular, and when it is singular ``s3`` is zero, so the
    term it multiplies vanishes either way.

    Clustering calls this on 7x3 arrays hundreds of thousands of times per
    bucket, where numpy's dispatch overhead outweighs the arithmetic, so the
    reductions are written as ``np.add.reduce`` and the determinant expanded by
    hand rather than routed through ``np.sum``/``np.linalg.det``.

    Algebraically identical to ``kabsch``'s ``ssd``, but not bit-identical: the
    precision differs by design -- this accumulates in float64, while ``kabsch``
    runs its SVD and residual in float32 because its ``R``/``t`` are consumed as
    float32 coordinates. Measured worst-case disagreement is ~2e-7 A in RMSD
    against a 0.5 A clustering threshold. Use ``kabsch`` when you need R or t.
    """
    X = np.asarray(X, dtype=np.float32)
    Y = np.asarray(Y, dtype=np.float32)

    if Y.ndim != 3:
        raise ValueError(f"kabsch_ssd: Y must have ndim=3, got Y.ndim={Y.ndim}, Y.shape={Y.shape}")
    if X.ndim not in (2, 3):
        raise ValueError(f"kabsch_ssd: X must have ndim 2 or 3, got X.ndim={X.ndim}")
    # Checked rather than left to the matmul: the centering below divides by Y's
    # point count, so a mismatch mis-centers X before anything raises.
    if X.shape[-2] != Y.shape[1]:
        raise ValueError(
            f"kabsch_ssd: point-count mismatch, X has {X.shape[-2]} points per structure, "
            f"Y has {Y.shape[1]}")
    if X.ndim == 3 and X.shape[0] != Y.shape[0]:
        raise ValueError(
            f"kabsch_ssd: batch mismatch, X has {X.shape[0]} structures, Y has {Y.shape[0]}")
    if not np.isfinite(X).all() or not np.isfinite(Y).all():
        raise ValueError(
            "Non-finite values detected in kabsch_ssd input. "
            f"X.shape={X.shape}; Y.shape={Y.shape}")

    M = Y.shape[0]
    if M == 0:
        return np.empty((0,), dtype=np.float64)

    n_pts = Y.shape[1]
    inv_n = 1.0 / n_pts
    fixed_X = X.ndim == 2
    if fixed_X:
        Xc = (X - np.add.reduce(X, axis=0) * inv_n).astype(np.float64)
        XcT = np.ascontiguousarray(Xc.T)
        x_norm = float(np.add.reduce(np.add.reduce(Xc * Xc)))

    chunks = []
    for start in range(0, M, chunk_size):
        Yc = Y[start:start + chunk_size].astype(np.float64)
        Yc -= np.add.reduce(Yc, axis=1)[:, None, :] * inv_n

        if fixed_X:
            H = np.matmul(XcT[None, :, :], Yc)
            xn = x_norm
        else:
            Xc_b = X[start:start + chunk_size].astype(np.float64)
            Xc_b -= np.add.reduce(Xc_b, axis=1)[:, None, :] * inv_n
            H = np.matmul(np.transpose(Xc_b, (0, 2, 1)), Yc)
            xn = np.add.reduce(np.add.reduce(Xc_b * Xc_b, axis=2), axis=1)

        sv = np.linalg.svd(H, compute_uv=False)
        d = np.where(_det3(H) < 0.0, -1.0, 1.0)
        y_norm = np.add.reduce(np.add.reduce(Yc * Yc, axis=2), axis=1)
        trace = sv[:, 0] + sv[:, 1] + d * sv[:, 2]
        # Clamped because the closed form, unlike an explicit residual, can land
        # a hair below zero when the two structures superimpose exactly.
        chunks.append(np.maximum(xn + y_norm - 2.0 * trace, 0.0))

    return np.concatenate(chunks, axis=0) if len(chunks) > 1 else chunks[0]


def kabsch(X, Y, chunk_size=30000):
    """Chunked Kabsch alignment under ``Y ~= X @ R + t``.

    ``X`` may be ``[N, 3]`` or ``[M, N, 3]``; ``Y`` must be ``[M, N, 3]``.
    Returns proper rotations ``R [M, 3, 3]``, translations ``t [M, 3]``, and
    squared deviations ``ssd [M]``.

    Rank-2 cross-covariance matrices are valid and common for three-point fits.
    At rank 0 or 1 the minimum SSD is still valid, but the optimal rotation is
    non-unique; callers must not extrapolate that rotation to points outside the
    fit unless they provide another non-collinear anchor.

    The SVD and residual run in float32, unlike ``kabsch_ssd``'s float64 -- see
    that docstring for why the two ``ssd`` values differ in the last digits.
    """
    X = np.asarray(X, dtype=np.float32)
    Y = np.asarray(Y, dtype=np.float32)

    if Y.ndim != 3:
        raise ValueError(f"kabsch: Y must have ndim=3, got Y.ndim={Y.ndim}, Y.shape={Y.shape}")

    if X.ndim not in (2, 3):
        raise ValueError(f"kabsch: X must have ndim 2 or 3, got X.ndim={X.ndim}")
    if X.shape[-2] != Y.shape[1]:
        raise ValueError(
            f"kabsch: point-count mismatch, X has {X.shape[-2]} points per structure, "
            f"Y has {Y.shape[1]}")
    if X.ndim == 3 and X.shape[0] != Y.shape[0]:
        raise ValueError(
            f"kabsch: batch mismatch, X has {X.shape[0]} structures, Y has {Y.shape[0]}")

    x_finite = np.isfinite(X)
    y_finite = np.isfinite(Y)
    if not x_finite.all() or not y_finite.all():
        x_invalid = np.argwhere(~x_finite)
        y_invalid = np.argwhere(~y_finite)
        x_msg = ("none" if x_invalid.size == 0
                 else f"first invalid X value at {tuple(x_invalid[0])}")
        y_msg = ("none" if y_invalid.size == 0
                 else f"first invalid Y value at {tuple(y_invalid[0])}")
        raise ValueError(
            "Non-finite values detected in kabsch input. "
            f"{x_msg}; {y_msg}; "
            f"X.shape={X.shape}; Y.shape={Y.shape}"
        )

    M = Y.shape[0]
    R_chunks = []
    t_chunks = []
    ssd_chunks = []

    # Fast path: one fixed X against many Y rows
    if X.ndim == 2:
        Xbar = X.mean(axis=0, keepdims=True)
        Xc = X - Xbar
        XcT = Xc.T

        for start in range(0, M, chunk_size):
            stop = min(start + chunk_size, M)
            Y_chunk = Y[start:stop]

            Ybar = Y_chunk.mean(axis=1, keepdims=True)
            Yc = Y_chunk - Ybar

            H = np.matmul(XcT[None, :, :], Yc)

            R = _proper_kabsch_rotations(H).astype(np.float32, copy=False)

            XRbar = np.matmul(Xbar[None, :, :], R)[:, 0, :]
            t = (Ybar[:, 0, :] - XRbar).astype(np.float32, copy=False)

            XR = np.matmul(Xc[None, :, :], R)
            diff = XR - Yc
            ssd = np.sum(diff * diff, axis=(1, 2), dtype=np.float64)

            R_chunks.append(R)
            t_chunks.append(t)
            ssd_chunks.append(ssd)

    # Fallback: batched X and batched Y
    else:
        for start in range(0, M, chunk_size):
            stop = min(start + chunk_size, M)
            X_chunk = X[start:stop]
            Y_chunk = Y[start:stop]

            Xbar = X_chunk.mean(axis=1, keepdims=True)
            Ybar = Y_chunk.mean(axis=1, keepdims=True)

            Xc = X_chunk - Xbar
            Yc = Y_chunk - Ybar

            H = np.matmul(np.transpose(Xc, (0, 2, 1)), Yc)

            R = _proper_kabsch_rotations(H).astype(np.float32, copy=False)
            t = (Ybar - np.matmul(Xbar, R)).reshape(-1, 3).astype(np.float32, copy=False)

            diff = np.matmul(Xc, R) - Yc
            ssd = np.sum(diff * diff, axis=(1, 2), dtype=np.float64)

            R_chunks.append(R)
            t_chunks.append(t)
            ssd_chunks.append(ssd)

    if R_chunks:
        R = np.concatenate(R_chunks, axis=0)
        t = np.concatenate(t_chunks, axis=0)
        ssd = np.concatenate(ssd_chunks, axis=0)
    else:
        R = np.empty((0, 3, 3), dtype=np.float32)
        t = np.empty((0, 3), dtype=np.float32)
        ssd = np.empty((0,), dtype=np.float64)

    return R, t, ssd

def convert_time_elapsed(seconds):
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = round(seconds % 60, 2)
    return h, m, s

def best_inplace_symmetry_rmsd(ref_mol, query_mol):
    """Return RDKit's symmetry-aware in-place RMSD without alignment.

    ``CalcRMS`` generates atom mappings internally and treats conjugated terminal
    groups symmetrically. In RDKit, ``maxMatches=0`` requests every mapping, so a
    match cap cannot silently omit the lowest-RMSD mapping. Hydrogen-free inputs
    avoid the most common source of combinatorial explosion.

    Those mappings are RDKit's own and deliberately do **not** follow
    ``identify_mol_automorphisms``. That policy merges terminal N/O/S because a
    *fragment*'s degree-1 atom may have a bond omitted by truncation; these inputs
    are whole ligands, where a terminal atom really is terminal, so importing it
    here would merge genuinely distinct atoms. Don't "reconcile" the two -- see
    docs/pitfalls.md, "vdG automorphisms and hit-finder CalcRMS use different
    policies, deliberately".
    """
    if ref_mol.GetNumAtoms() != query_mol.GetNumAtoms():
        raise ValueError(f"[ERROR] Atom count mismatch: ref={ref_mol.GetNumAtoms()}, "
            f"query={query_mol.GetNumAtoms()}")
    _require_hydrogen_free(ref_mol, "In-place symmetry RMSD")
    _require_hydrogen_free(query_mol, "In-place symmetry RMSD")
    return MA.CalcRMS(
        query_mol,
        ref_mol,
        maxMatches=0,
        symmetrizeConjugatedTerminalGroups=True,
    )


# Underscore is not a valid SMILES character, so _XX_ encodings are unambiguous.
_SMILES_TO_FILENAME = {'/': '_fs_', '\\': '_bs_'}
_SMILES_TO_JOB_NAME = {'#': '_tp_', '/': '_fs_', '\\': '_bs_'}  # these truncate or break SGE directives

def normalize_rmsd(num_atoms, atoms):
    '''Return a size-normalized RMSD threshold (Å) for the given atom set
    ('cgvdmbb' or 'flankbb'). The threshold scales linearly with the number
    of atoms between 8 and 15:
        flankbb: 0.5 Å → 1.5 Å
        cgvdmbb: 0.5 Å → 1.0 Å
    Below 8 atoms, use the minimum; above 15, use the maximum.
    '''
    if atoms == 'flankbb':
        max_threshold, min_threshold = 1.5, 0.5
    elif atoms == 'cgvdmbb':
        max_threshold, min_threshold = 1.0, 0.5
    else:
        raise ValueError(f"Unknown atom set for normalize_rmsd: {atoms}")
    min_atoms, max_atoms = 8, 15
    if num_atoms < min_atoms:
        return min_threshold
    if num_atoms > max_atoms:
        return max_threshold
    return min_threshold + (num_atoms - min_atoms) / (max_atoms - min_atoms) * (max_threshold - min_threshold)


# One atom per match, in SMILES order. The bracket branch comes first so that an
# isotope/charge/chirality/H spec inside [] never yields extra tokens: the symbol
# is the leading letter pair, and a trailing uppercase H ([nH], [SeH]) is a
# hydrogen count, not an atom. Outside brackets only the organic subset is legal,
# with Br/Cl before the single letters so they are not split.
_SMILES_ATOM_RE = re.compile(
    r"\[\d*(?P<bracket>[A-Za-z][a-z]?|\*)"
    r"|(?P<organic>Br|Cl|[BCNOFPSI]|[bcnosp])")


def extract_elements(smiles: str):
    """Return element symbols in SMILES order, one per atom.

    Aromatic atoms keep their lowercase form ('c', 'se'); callers that compare
    against RDKit symbols capitalize first. A bracketed query token carrying no
    element symbol (SMARTS '[+]', '[#6]') yields nothing, so this is a SMILES
    tokenizer, not a SMARTS one.
    """
    return [m.group('bracket') or m.group('organic')
            for m in _SMILES_ATOM_RE.finditer(smiles)]

def smiles_to_filename(smiles):
    '''Encode SMILES into a string safe for use as a filename (encodes / and \\).'''
    return ''.join(_SMILES_TO_FILENAME.get(c, c) for c in smiles)

def filename_to_smiles(name):
    '''Inverse of smiles_to_filename(). Injective because '_' is not a SMILES
    character, so an encoded library name round-trips exactly.'''
    for ch, enc in _SMILES_TO_FILENAME.items():
        name = name.replace(enc, ch)
    return name

def smiles_to_job_name(smiles):
    '''Encode SMILES into a string safe for SGE job names (encodes # which truncates directives).'''
    return ''.join(_SMILES_TO_JOB_NAME.get(c, c) for c in smiles)
