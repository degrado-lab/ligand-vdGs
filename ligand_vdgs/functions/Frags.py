import os
import itertools
import json
import re
import io
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import defaultdict
import prody as pr


def get_fragments(bond_radius, mol, min_frag_size=4, max_frag_size=5, quiet=True):
    """
    Decompose a ligand into fragments and group by binding site.

    Returns dict of {sub_smiles: [[site group 1], [site group 2], ...]} where each
    site group contains permutations of the same atoms at the same ligand site.
    """
    # Decompose the ligand into fragments and store the fragment SMILES. Use SMILES
    # instead of SMARTS b/c only SMILES (from rdkit) differentiates aliphatic and
    # aryl (C,c vs. [#6]). Fragment on bond radii `bond_radius` AND the postive
    # integers less than `bond_radius`, because for example, drugs containing
    # sulfonamide might produce only 6-atom sulfonamides and not CS(N)(=O)=O. 
    filtered_frags = {} # key=sub_smiles, value = list of Mol objs
    substructs = []
    # Computed once and shared across radii: it is a property of the parent mol,
    # which does not change here.
    ring_queries = parent_ring_queries(mol)
    for rad in list(range(1, bond_radius + 1)):
        substructs += fragment_on_bond_d(mol, rad, min_frag_size, max_frag_size,
                                         ring_queries)
    # Add substructs to `filtered_frags`.
    for orig_sub, orig_mol_inds in substructs: # contains H's that need to be scrubbed.
        # apply size threshold.
        frag_size = orig_sub.GetNumHeavyAtoms()
        if frag_size < min_frag_size or frag_size > max_frag_size:
            continue
        _results = manually_remove_Hs(orig_sub, return_single_mol_or_perms='perms')
        if _results is None:
            continue

        sub_perms_mols_inds, sub_smiles = _results
        for sub, perm_inds in sub_perms_mols_inds: 
            # skip if the heavy atoms are all carbons, or if there are none at all
            # (an all-hydrogen graph is not a CG either). Count atoms, not characters
            # of the concatenated symbols, so 2-letter elements (Cl, Ca, Sc, ...)
            # aren't miscounted as carbons.
            heavy_atoms = [a for a in sub.GetAtoms() if a.GetAtomicNum() != 1]
            if not heavy_atoms or all(a.GetAtomicNum() == 6 for a in heavy_atoms):
                continue
            # add to dict
            substruct_data = (sub, perm_inds, orig_mol_inds)
            if sub_smiles not in filtered_frags:
                filtered_frags[sub_smiles] = [substruct_data]
            else:
                # duplicates within same lig obj occur for unknown reasons; skip them.
                already_present = False
                for existing_sub, existing_sub_inds, existing_orig_mol_inds in filtered_frags[sub_smiles]:
                    if (perm_inds, orig_mol_inds) == (existing_sub_inds, existing_orig_mol_inds):
                        already_present = True
                        break
                if not already_present:
                    filtered_frags[sub_smiles].append(substruct_data)

    # Group substructs by instances (sites) in the ligand. For example, if a substruct has
    # orig_mol_inds [14, 15, 16, 17] and another has [16, 17, 18, 19], they are the same 
    # site. The purpose of grouping is to prevent overcounting of matches when there are diff 
    # permutations of the same atoms or slight overlap of sites.
    grouped_frags = {}
    for sub_smiles, substruct_data in filtered_frags.items():
        groups = group_lig_sites_by_overlap(substruct_data)
        grouped_frags[sub_smiles] = groups
        if not quiet: 
            print(f"Fragment: {sub_smiles}, # perms: {len(substruct_data)}, "
                  f"# sites: {len(groups)}", flush=True)
    return grouped_frags

def group_lig_sites_by_overlap(data, key_index=2, threshold=0.5):
    """
    Input: list of tuples (substruct Mol obj, perm_inds, orig_mol_inds) that describe 
        instances of a frag in a lig. Determine whether each site has >1 instances, as 
        determined by sharing (overlapping) >=1/2 of the atoms b/n one frag instance and 
        another. For example, if one instance has orig_mol_inds [14, 15, 16, 17] and 
        another has [16, 17, 18, 19], they are the considered the same site on the lig.
        This avoids overcounting of matches when there are diff permutations of the CG.
    Output: list of groups.
    Method: Two frag instances are related if |A ∩ B| >= threshold * max(|A|, |B|).
    """
    
    sets = [set(item[key_index]) for item in data]
    n = len(data)
    # Union–find
    parent = list(range(n))
    rank = [0]*n
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb: 
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1
    # Connect pairs that overlap 
    for i in range(n):
        for j in range(i+1, n):
            inter = len(sets[i] & sets[j])
            denom = max(len(sets[i]), len(sets[j])) # >= half
            if inter >= threshold * denom:
                union(i, j)
    # Collect groups
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(data[i])

    return list(groups.values())

# The property `fragment_on_bond_d` writes on each submol atom, carrying the
# ring class the atom had *in its parent*. It has to travel as a property
# because PathToSubmol cuts ring bonds: a pyridine-derived fragment is a chain
# in the submol, so by the time the key is written the ring is gone.
RING_QUERY_PROP = '_vdgRingQuery'


def ring_query_for_atom(atom):
    """SMARTS ring primitive for an atom's ring context, or None to leave it bare.

    Three regimes, each decided by measurement on a built library:

    * Aromatic atoms: None. Lower case already implies a ring, and 5- vs 6-ring
      aromatics are deliberately pooled -- the protein-side divergence for
      imidazole vs pyrimidine is only 2-3x the same-population null, the same
      donor/acceptor blur already exists inside imidazole through tautomer
      normalization, and encoding size there shatters every fused aromatic
      (indole's fusion carbons go r5 and r6) for ~150 extra jobs.
    * Saturated ring atoms: `r<n>` for the smallest ring containing them. Size
      is kept here because the largest effect in the whole library is aliphatic:
      furanose vs pyranose (`CCCOC`, ribose vs glucose) diverges 26x null in
      protein partners -- LEU/ILE/VAL nucleotide pockets vs the TRP/ASN/ASP
      carbohydrate triad -- and membership alone puts both in one bucket. Both
      halves are enormous, so the split costs nothing in sparsity. Fused
      aliphatic shattering is small (26 mixed keys at threshold) because the
      all-carbon filter already removes ring skeletons; steroids yield only
      their heteroatom-bearing positions.
    * Acyclic atoms: `!R`.

    Recorded at all because RMSD cannot recover ring context: 100% of
    chain-context `CCNC` vdGs have a ring-context one inside the hit threshold,
    with the cross-population nearest-neighbour distance (0.037 A) inside the
    within-population one (0.026 A). No cutoff separates them.

    Purely topological -- it reads the bond graph, not the aromaticity model --
    so it is immune to kekulization and tautomer choice. RDKit (query side) and
    OpenBabel (mining side) were checked to agree on smallest ring size across
    12 bridged, fused and spiro systems, cubane and morphinan included.
    """
    if atom.GetIsAromatic():
        return None
    if not atom.IsInRing():
        return '!R'
    # Smallest ring, uncapped. There is deliberately no plain-`R` fallback for
    # large rings: SMARTS `R` matches *any* ring atom, so a macrocycle-derived
    # `[C;R][C;R][O;R][C;R][C;R]` key (42 CCD ligands) would, as a mining query,
    # pull in every THF and sugar in the database -- estimated 27,399 structures --
    # and 266 such keys cleared the build threshold as superset jobs. `r<n>` is
    # exact per ring size, so a macrocycle key mines only macrocycles, and the
    # long tail of large-ring keys simply falls below the structure threshold.
    ring_info = atom.GetOwningMol().GetRingInfo()
    return f'r{min(ring_info.AtomRingSizes(atom.GetIdx()))}'


# Compiled once: these run per bracket atom per fragment (~1M re module-level
# dispatches per 3k ligands otherwise), and re's internal cache lookup was
# measurable against the actual matching.
_BRACKET_ATOM = re.compile(r'\[([^\]]+)\]')
_BRACKET_ATOM_LOOSE = re.compile(r'\[[^\]]*\]')
_EXPLICIT_H = re.compile(r'H(?![a-z])\d*')
_NON_LETTER = re.compile(r'[^A-Za-z]')


def _remove_bracket_hydrogens_in_smiles(smiles: str, ring_queries=None) -> str:
    """
    Strip explicit hydrogens from bracket atoms in a SMILES string (string-only).
    - [CH3] -> C
    - [NH2] -> N
    - [NH3+] -> [N+]
    - [H]   -> '' (remove the token)
    Brackets are dropped only if the remaining content is a plain element/atom symbol
    (letters only, e.g., C, N, n). If charge/chirality/isotope/class remains,
    brackets are kept.

    `ring_queries`, when given, is one entry per bracket atom in written order
    (see ring_query_for_atom): a SMARTS ring primitive, or None to leave that
    atom bare. A primitive is appended inside its atom's brackets and the
    brackets are then kept -- `[C;!R][C;!R][N;R]` is a different query from
    `CCN` and must stay one. The caller owns the ordering;
    MolToSmiles(allHsExplicit=True) brackets every atom, which is what makes
    bracket order and `_smilesAtomOutputOrder` line up.
    """
    counter = itertools.count()

    # Rewrite every bracket atom content. Standalone [H] is dropped here, not in a
    # pre-pass: it is a bracket atom like any other, so it must consume its
    # ring_queries slot or every following atom gets the wrong ring context.
    def _rewrite(m):
        inner = m.group(1)
        idx = next(counter)

        if inner == 'H':
            return ''

        # Remove explicit H or H<number> that is not part of an element symbol like Hg/He
        cleaned = _EXPLICIT_H.sub('', inner)

        ring_query = None if ring_queries is None else ring_queries[idx]
        if ring_query is not None:
            return f'[{cleaned};{ring_query}]'

        # If anything special remains, we must keep the brackets (e.g., [N+], [C@H],
        # [13C]). One test covers all of them: charge, chirality, isotope/class digits
        # and ':' are all non-letters.
        if _NON_LETTER.search(cleaned):
            return f'[{cleaned}]'

        # Keep brackets for multi-letter element symbols (e.g. [Se], [Te], [Sn]).
        # Single-letter content (C, N, n, ...) can safely lose its brackets.
        if len(cleaned) > 1:
            return f'[{cleaned}]'
        return cleaned

    s = _BRACKET_ATOM.sub(_rewrite, smiles)
    return s


class _RingAnnotationUnmappable(Exception):
    """Raised when an annotated fragment's ring context can't be placed on its SMILES.

    Distinct from "this molecule was never annotated": the bare key is a strictly
    looser query and is string-identical to a legitimate aromatic-only key, so
    silently writing one would nest that fragment's library inside another's.
    Callers on the library path must drop the fragment instead.
    """


def _ring_queries_in_written_order(mol, smiles_w_H):
    """Ring primitives ordered to match the bracket atoms of `smiles_w_H`.

    Returns None only for molecules that carry no ring annotation at all -- the
    hit-finding query path, which legitimately writes an unannotated key. If the
    molecule *is* annotated but the correspondence can't be established, raises
    `_RingAnnotationUnmappable` rather than guessing at it or degrading to a bare
    key: a misaligned annotation would silently attach one atom's ring context to
    another, producing a plausible key for the wrong chemistry.
    """
    if not mol.GetNumAtoms() or not mol.GetAtomWithIdx(0).HasProp(RING_QUERY_PROP):
        return None
    props = mol.GetPropsAsDict(includePrivate=True, includeComputed=True)
    order = props.get('_smilesAtomOutputOrder')
    if order is None:
        raise _RingAnnotationUnmappable('MolToSmiles did not record '
                                        '_smilesAtomOutputOrder')
    if isinstance(order, str):
        order = json.loads(order)
    order = list(order)
    # allHsExplicit=True brackets every atom, so these must agree; if they do
    # not, the SMILES was not written the way this assumes.
    if len(order) != mol.GetNumAtoms():
        raise _RingAnnotationUnmappable(
            f'output order covers {len(order)} of {mol.GetNumAtoms()} atoms')
    n_brackets = len(_BRACKET_ATOM_LOOSE.findall(smiles_w_H))
    if n_brackets != len(order):
        raise _RingAnnotationUnmappable(
            f'{n_brackets} bracket atoms in {smiles_w_H!r} for {len(order)} atoms')
    queries = []
    for atom_idx in order:
        atom = mol.GetAtomWithIdx(atom_idx)
        if not atom.HasProp(RING_QUERY_PROP):
            raise _RingAnnotationUnmappable('fragment is only partially annotated')
        queries.append(atom.GetProp(RING_QUERY_PROP) or None)
    return queries


# Cap on the atom-order permutations enumerated per fragment. See the warning in
# manually_remove_Hs: hitting it drops the fragment rather than truncating it.
_MAX_CG_ATOM_PERMS = 1000


def manually_remove_Hs(orig_substruct, return_single_mol_or_perms):
    '''return_single_mol_or_perms must be 'single' (for processing a whole ligand) or 
    'perms' (when processing frags).

    Returns (perm(s), smiles_no_Hs), or None if the H-stripped SMILES can't be mapped
    back onto the input. `perm_inds` index into the H-stripped molecule (identical to
    the input when the input had no hydrogens in its graph).'''
    if return_single_mol_or_perms not in ('single', 'perms'):
        raise ValueError("return_single_mol_or_perms must be 'single' or 'perms', "
                         f"got {return_single_mol_or_perms!r}")
    if orig_substruct is None:
        return None

    # Drop hydrogens that are real atoms in the graph (SMILES written as [H]O..., PDB
    # blocks parsed with removeHs=False). They aren't part of the CG, and leaving them
    # in makes the SMARTS below match only a subset of the atoms, which RenumberAtoms
    # rejects -- silently discarding the whole molecule. RemoveAllHs (unlike RemoveHs)
    # also drops charged/isotopic H's, and preserves conformer coords of the atoms kept.
    try:
        substruct = Chem.RemoveAllHs(orig_substruct, sanitize=False)
    except Exception as e:
        # Don't fall back to the H-bearing graph: it defeats the guarantee above,
        # and the fragment would be dropped a few lines down anyway (the H-free
        # SMARTS then covers only a subset of the atoms) with no reason given.
        print(f'[WARNING] manually_remove_Hs: RemoveAllHs failed ({e}); '
              'dropping fragment.', flush=True)
        return None
    if substruct.GetNumAtoms() == 0:
        return None

    # Remove hydrogens. RDKit docs say that Chem.RemoveHs() implicit and explicit are removed, 
    # but this isn't true for [nH], [OH], [Ho], etc. so need to manually remove H's. 
    # First, export SMILES *with explicit Hs shown* so that patterns like [NH3+] are present, 
    # then regex-strip bracket hydrogens while keeping charge and other annotations.
    # Caveat: `substruct` was never sanitized (PathToSubmol output), so its
    # aromaticity and valence flags may be stale. MolToSmiles can then fail, or
    # write a string RDKit itself will not re-parse below -- both end in a
    # dropped fragment, which is why each is reported.
    try:
        smiles_w_H = Chem.MolToSmiles(substruct, allHsExplicit=True, 
                                  isomericSmiles=False) # shows [NH3+], [CH3], etc.
    except Exception as e:
        print(f'[WARNING] manually_remove_Hs: MolToSmiles failed on an unsanitized '
              f'fragment ({e}); dropping it.', flush=True)
        return None

    # Two strings, because they answer two different questions.
    #
    #   smiles_no_Hs   the plain H-stripped SMARTS, matched back against `substruct`
    #                  below to recover atom order.
    #   annotated_key  the same query plus each atom's ring context *in its parent*,
    #                  which is what gets returned as the fragment's identity.
    #
    # They must stay separate: PathToSubmol cuts ring bonds, so a pyridine-derived
    # fragment is an acyclic chain inside `substruct`. Matching `[c;r6]...` against
    # it would find nothing and drop every ring fragment in the library. The ring
    # context is real at mining time, where the query runs against the whole
    # ligand -- which is exactly where the key is used.
    smiles_no_Hs = _remove_bracket_hydrogens_in_smiles(smiles_w_H)
    try:
        ring_queries = _ring_queries_in_written_order(substruct, smiles_w_H)
    except _RingAnnotationUnmappable as e:
        # The fragment was annotated but the annotation can't be placed. Falling
        # back to `smiles_no_Hs` would file it under a bare key that is both a
        # looser query and indistinguishable from a legitimate aromatic-only one.
        print(f'[WARNING] manually_remove_Hs: ring annotation for {smiles_no_Hs!r} '
              f'could not be mapped onto the written SMILES ({e}); dropping '
              'fragment.', flush=True)
        return None
    annotated_key = (smiles_no_Hs if ring_queries is None
                     else _remove_bracket_hydrogens_in_smiles(smiles_w_H, ring_queries))

    # Complication: after converting the Mol obj to smiles, the Mol obj won't have the same 
    # atom order as the original Mol obj, which is important when extract CG coords. 
    # We need to rearrange the atom order of substruct by getting the atom indices in 
    # the orig full molecule.
    # -- Parse the SMILES back into a new molecule
    mol_from_smarts = Chem.MolFromSmarts(smiles_no_Hs)
    if mol_from_smarts is None or mol_from_smarts.GetNumAtoms() == 0:
        # Both halves are reachable: stale flags can produce a SMILES that will
        # not re-parse, and an all-hydrogen input ('[H][H]') strips to '', which
        # parses to an empty query. A lone '[H+]' counterion is *not* one of
        # these -- it strips to '[+]', which RDKit accepts as a charge-only query.
        print(f'[WARNING] manually_remove_Hs: H-stripped SMARTS {smiles_no_Hs!r} is '
              'unusable; dropping fragment.', flush=True)
        return None
    # -- Map original atoms to new atom order using substructure matching.
    #    'single' only ever uses the first perm, so don't enumerate the rest -- symmetric
    #    ligands can generate a huge number of matches.
    max_matches = 1 if return_single_mol_or_perms == 'single' else _MAX_CG_ATOM_PERMS
    try:
        matches_to_map_to_sub = substruct.GetSubstructMatches(mol_from_smarts,
            uniquify=False, # returns all permutations of order of atom inds
            maxMatches=max_matches)
    except Exception as e:
        print(f'[WARNING] manually_remove_Hs: substructure match failed for '
              f'{smiles_no_Hs!r} ({e}); dropping fragment.', flush=True)
        return None
    if (return_single_mol_or_perms == 'perms'
            and len(matches_to_map_to_sub) >= _MAX_CG_ATOM_PERMS):
        # A truncated match set is a silently incomplete permutation set, so drop
        # the fragment loudly instead of mining part of its symmetry. Unreachable
        # for the 4-5 heavy-atom fragments this pipeline mines (|Aut| tops out at
        # 24), so reaching it means an assumption broke, not that the cap is low.
        print(f'[WARNING] manually_remove_Hs: {smiles_no_Hs!r} hit the '
              f'{_MAX_CG_ATOM_PERMS}-match cap, so its permutation set would be '
              'incomplete; dropping fragment.', flush=True)
        return None
    # -- Reorder atoms in the H-stripped mol. RenumberAtoms needs a full permutation, so
    #    partial matches (SMARTS covering only some atoms) can't be used.
    n_atoms = substruct.GetNumAtoms()
    cg_atom_perms = []
    for perm_inds in list(matches_to_map_to_sub):
        if len(perm_inds) != n_atoms:
            continue
        mol_copy = Chem.Mol(substruct) # for immutable deep copy; don't use copy.deepcopy()
        try:
            renumbered_substruct = Chem.RenumberAtoms(mol_copy, perm_inds)
        except Exception:
            continue
        cg_atom_perms.append((renumbered_substruct, perm_inds))
    if len(cg_atom_perms) == 0:
        print(f'[WARNING] manually_remove_Hs: {smiles_no_Hs!r} produced no '
              'whole-molecule atom mapping back onto the fragment; dropping it.',
              flush=True)
        return None
    # The annotated key is what the caller stores and mines with; smiles_no_Hs was
    # only ever the vehicle for recovering atom order above.
    if return_single_mol_or_perms == 'single':
        return cg_atom_perms[0], annotated_key # return the first permutation only
    return cg_atom_perms, annotated_key # return all permutations. `cg_atom_perms` is a 
                                        # list of (substructure Mol objs, perm_inds)

def parent_ring_queries(mol):
    """`ring_query_for_atom` for every atom of `mol`, indexed by atom index.

    Hoisted out of `fragment_on_bond_d`: the answer depends only on the parent
    atom, but the naive version recomputed it once per (center atom, radius,
    submol atom) -- 708k calls for 3k ligands where 119k distinct ones exist.
    '' means "leave bare", matching the property convention below.
    """
    return [ring_query_for_atom(atom) or '' for atom in mol.GetAtoms()]


def fragment_on_bond_d(mol, radius, min_frag_size=None, max_frag_size=None,
                       ring_queries=None):
    # Code from https://iwatobipen.wordpress.com/2020/08/12/get-and-draw-molecular-fragment-with-user-defined-path-rdkit-memo/
    #
    # `min_frag_size`/`max_frag_size` apply the caller's heavy-atom window here
    # rather than after the fact, so out-of-window submols are dropped before
    # being annotated. Same predicate on the same object as the caller's own
    # check, so the surviving set is unchanged; it only avoids annotating
    # submols that were about to be discarded. None disables the filter.
    atoms = mol.GetAtoms()
    if ring_queries is None:
        ring_queries = parent_ring_queries(mol)
    submols = []
    for atom in atoms:
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom.GetIdx(), 
            enforceSize=False) # don't enforce size to also get frags that are < radius away
        amap = {}
        submol = Chem.PathToSubmol(mol, env, atomMap=amap)
        if min_frag_size is not None or max_frag_size is not None:
            frag_size = submol.GetNumHeavyAtoms()
            if ((min_frag_size is not None and frag_size < min_frag_size) or
                    (max_frag_size is not None and frag_size > max_frag_size)):
                continue
        # PathToSubmol cuts ring bonds, so the submol has already forgotten which
        # of its atoms came from a ring. Record it now, from the parent, keyed by
        # amap (original index -> submol index; the sorted key list below loses
        # that correspondence).
        for orig_idx, sub_idx in amap.items():
            # '' means "leave bare"; the property must be set on every atom so
            # _ring_queries_in_written_order can tell a real annotation set from
            # a mol built by a path that never annotated at all.
            submol.GetAtomWithIdx(sub_idx).SetProp(RING_QUERY_PROP,
                                                   ring_queries[orig_idx])
        # Store the submol and its atom indices in the original mol to ensure that if there 
        # are >1 instances of a frag in a single ligand, they won't get skipped
        orig_inds = sorted(amap.keys())
        submols.append((submol, orig_inds))
    return submols

def is_organic(mol):
    '''Requires a carbon atom. Checking for any of C/O/N (as before) let small
    non-carbon ions like sulfate/phosphate through, since they fit the default
    4-5 heavy-atom fragment size window and contain O.'''
    return any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms())

def get_query_ligand_mol(query_struct_or_path, query_lig_smiles):
    """
    Hit-finding only. Build the single query ligand's RDKit Mol from either a prody
    obj or a path to a pdb/pdb.gz/cif file for a query protein-ligand structure.

    This hard-requires exactly one ligand residue instance and raises otherwise --
    it is not a general multi-ligand-structure parser.

    Deliberately does *not* fragment: hit finding enumerates the library's own
    fragment keys against this whole ligand (match_library_frags_to_query in
    hit_finder_core) instead of re-running generation-time enumeration on the
    query, whose radius/size parameters the query side cannot know.

    Returns the H-stripped ligand Mol (single permutation, with coords), or None
    if it could not be built (the reason is printed).
    """

    if isinstance(query_struct_or_path, str):
        if query_struct_or_path.endswith('.pdb') or query_struct_or_path.endswith('.pdb.gz'):
            query_struct = pr.parsePDB(query_struct_or_path)
        elif query_struct_or_path.endswith('.cif') or query_struct_or_path.endswith('.cif.gz'):
            query_struct = pr.parseCIF(query_struct_or_path)
        else:
            raise ValueError(f"Unsupported file format: {query_struct_or_path}")
    else:
        # assume prody obj
        query_struct = query_struct_or_path

    # Identify the single query ligand and fragment it
    query_hetatms = query_struct.select(
        'hetatm and not (ion or water or resname SEP or resname TPO or resname MSE)')

    if query_hetatms is None:
        raise ValueError(
            "get_query_ligand_mol: could not find ligand HETATM atoms.")
    query_resnames = sorted(set(query_hetatms.getResnames()))
    if len(query_resnames) != 1:
        # Disambiguate by matching heavy-atom count per residue to the query SMILES.
        # This naturally handles crystallographic waters (1 heavy atom), buffer
        # molecules (GOL, EDO, etc.), and other co-solvent HETATMs.
        query_lig_template_for_count = Chem.MolFromSmiles(query_lig_smiles)
        if query_lig_template_for_count is not None:
            n_heavy = query_lig_template_for_count.GetNumAtoms()
            matches = []
            for name in query_resnames:
                sel = query_hetatms.select(f"resname {name}")
                if sel is None:
                    continue
                n_residues = len(set(sel.getResindices().tolist()))
                sel_heavy = sel.select("not element H D")
                n_atoms = sel_heavy.numAtoms() if sel_heavy is not None else 0
                if n_residues > 0 and n_atoms == n_heavy * n_residues:
                    matches.append(name)
            if len(matches) == 1:
                query_hetatms = query_hetatms.select(f"resname {matches[0]}")
                query_resnames = matches
    if len(query_resnames) != 1:
        raise ValueError(
            f"get_query_ligand_mol: expected exactly one ligand "
            f"residue, got: {query_resnames}")

    # A single resname can still span multiple residue instances (e.g. a
    # polysaccharide's repeated monomer, or an NCS-duplicated ligand chain) --
    # that's covalently or spatially linked HETATMs getting merged into one
    # ligand, not one ligand. Fail hard rather than silently building a
    # cross-residue "molecule" from unrelated atoms. This is a hit-finding-only
    # requirement: a query structure must resolve to exactly one ligand residue.
    n_query_ligand_residues = len(set(query_hetatms.getResindices().tolist()))
    if n_query_ligand_residues != 1:
        raise ValueError(
            f"get_query_ligand_mol: expected exactly one ligand "
            f"residue instance of {query_resnames[0]!r}, got {n_query_ligand_residues} "
            "residues.")

    # Convert from ProDy object to RDKit Mol via an in-memory PDB block
    buf = io.StringIO()
    pr.writePDBStream(buf, query_hetatms)
    pdb_block = buf.getvalue()
    # Caveat: PDB files carry no bond records for ligands, so RDKit infers the
    # graph from interatomic distances. AssignBondOrdersFromTemplate below only
    # raises when the inferred graph cannot be matched to the template at all --
    # a mis-bonded graph that is still isomorphic to it passes silently, with
    # atom identities swapped.
    query_pdb_mol = Chem.MolFromPDBBlock(pdb_block, removeHs=True)
    query_lig_template = Chem.MolFromSmiles(query_lig_smiles)
    if query_pdb_mol is None or query_lig_template is None:
        print(f"[ERROR] Processing {query_struct_or_path}: RDKit failed to parse molecule",
              flush=True)
        return None

    # AssignBondOrdersFromTemplate raises only when the template matches *nothing*.
    # A template smaller than the selection matches a subgraph and is accepted, so
    # only that subgraph's bonds get orders -- which is also what happens when the
    # selection holds several copies of the ligand. Warn rather than raise: the
    # count is informative, but the copies are legitimate structures.
    if query_pdb_mol.GetNumAtoms() != query_lig_template.GetNumAtoms():
        print(f"[WARNING] Processing {query_struct_or_path}: ligand selection has "
              f"{query_pdb_mol.GetNumAtoms()} heavy atoms but the query SMILES has "
              f"{query_lig_template.GetNumAtoms()}; bond orders are assigned only to the "
              "atoms the template matches, and the rest stay single-bonded.",
              flush=True)

    # Assign bond orders (valence + aromaticity). Catch broadly: RDKit signals
    # failure here through several unrelated types (its sanitization exceptions
    # subclass ValueError, but the C++ layers also raise bare RuntimeError), and
    # this function's contract is to report and skip a structure, not to abort
    # the caller's walk over a database.
    try:
        query_pdb_mol_assigned_bonds = AllChem.AssignBondOrdersFromTemplate(
            query_lig_template, query_pdb_mol)
        stripped = manually_remove_Hs(query_pdb_mol_assigned_bonds,
                                      return_single_mol_or_perms='single')
        if stripped is None:
            print(f"[ERROR] Processing {query_struct_or_path}: could not strip hydrogens "
                  "from the ligand (see the warning above)", flush=True)
            return None
        (query_pdb_mol_assigned_bonds_no_H, _pdb_mol_perm_inds), _smiles_no_H = stripped
        # Ring perception, not decoration: fragment keys carry `r<n>` primitives
        # whose n is the *smallest* ring size, taken at generation time from a
        # sanitized parent. Matching those keys against this ligand needs the
        # same notion of ring size, and manually_remove_Hs builds its result with
        # sanitize=False. GetSymmSSSR supplies it without re-sanitizing a mol
        # whose valences came from a PDB.
        Chem.GetSymmSSSR(query_pdb_mol_assigned_bonds_no_H)
        return query_pdb_mol_assigned_bonds_no_H
    except Exception as e:
        print(
            f"[ERROR] Processing {query_struct_or_path}: "
            f"ligand hydrogen removal failed: {type(e).__name__}: {e}",
            flush=True,
        )
        return None

def submol_from_match(mol, match):
    """Build the fragment Mol for one substructure match, in the match's atom order.

    Atom i of the result is ``mol``'s atom ``match[i]``, so the fragment already
    carries the atom order of the pattern that produced the match -- which is
    what the library's CG coords are stored in, so no post-hoc reordering (and
    no re-matching a ring primitive against a ring-cut submol) is needed.

    The parent's conformer positions come along; bonds are copied for the pairs
    of matched atoms that are bonded in the parent. The result is deliberately
    left unsanitized: it is a cut-out of a real molecule with open valences, and
    its only consumers read coordinates and element symbols.
    """
    conf_src = mol.GetConformer()
    em = Chem.RWMol()
    conf = Chem.Conformer(len(match))
    for new_idx, orig_idx in enumerate(match):
        em.AddAtom(Chem.Atom(mol.GetAtomWithIdx(orig_idx)))
        conf.SetAtomPosition(new_idx, conf_src.GetAtomPosition(orig_idx))
    positions = {orig_idx: new_idx for new_idx, orig_idx in enumerate(match)}
    for i, orig_i in enumerate(match):
        for bond in mol.GetAtomWithIdx(orig_i).GetBonds():
            other = bond.GetOtherAtomIdx(orig_i)
            j = positions.get(other)
            # `j > i` also skips the second visit of each in-fragment bond.
            if j is not None and j > i:
                em.AddBond(i, j, bond.GetBondType())
    sub = em.GetMol()
    sub.AddConformer(conf, assignId=True)
    return sub


def check_vdg_job_status(sub_smiles, vdg_lib_dir):
    # Check if the vdg generation job finished without issues
    vdg_log_file = os.path.join(vdg_lib_dir, sub_smiles, f'{sub_smiles}_log')
    if not os.path.exists(vdg_log_file):
        return False
    with open(vdg_log_file, 'r') as f:
        log_contents = f.read()
    return 'Job completed.' in log_contents

def summarize_frags(frags_in_lib, frags_to_exclude, frags_to_include, logfile_fh):
    groups = {
        "Excluded": [],
        "Not in include list": [],
        "In vdg lib but incomplete": [],
        "In vdg lib and in include list": []}

    for sub_smiles in frags_in_lib:
        if sub_smiles in frags_to_exclude:
            groups["Excluded"].append(sub_smiles)
        elif (frags_to_include != 'all' and isinstance(frags_to_include, list) and 
              sub_smiles not in frags_to_include):
            groups["Not in include list"].append(sub_smiles)
        elif not frags_in_lib[sub_smiles]:
            groups["In vdg lib but incomplete"].append(sub_smiles)
        else:
            groups["In vdg lib and in include list"].append(sub_smiles)

    # Log summary
    if groups:
        print("\n--- Fragment vdG Library Summary ---", file=logfile_fh)
        for category, frags in groups.items():
            if frags:
                print(f"{category}:", file=logfile_fh)
                print("   ", ", ".join(frags), file=logfile_fh)
        print("\n", file=logfile_fh)
