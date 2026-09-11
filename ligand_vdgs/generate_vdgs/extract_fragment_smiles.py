'''
Extract qualifying fragments from database_frags_dict.pkl as a scheduler-agnostic
one-column work list. The keys are the annotated fragment SMARTS written by
fragment_database_ligs.py -- canonical SMILES carrying a per-atom ring
primitive (`[C;r5]`, `[C;!R]`), so they are SMARTS and not valid SMILES;
downstream they are parsed with utils.mol_from_fragment. Exact automorphisms are intentionally not
serialized here: every execution path derives them from the fragment in the
clustering core.

Fragments qualify if the parent database offers at least --min-instances candidate
vdG sites for them (SMARTS matches summed over every ligand copy, not structures
and not CCD entries), they contain no more than --max-size heavy atoms, and they
are not crystallographic solvent artifacts. The count is measured on the fly by
sampling --pdb-dir rather than read from a stored table, so a work list depends on
nothing but the pdb parent db and the fragment dict. These are the same thresholds and
selection rule used by make_sge_scripts_for_frags.py.

Usage (--pdb-dir is required: the counts are sampled from the parent db):
    python ligand_vdgs/generate_vdgs/extract_fragment_smiles.py \
        --pdb-dir <path/to/pdb_database/>
    python ligand_vdgs/generate_vdgs/extract_fragment_smiles.py \
        --pdb-dir <path/to/pdb_database/> \
        --frags-dict <path/to/database_frags_dict.pkl> \
        --output <path/to/fragment_work_list.txt>

Consuming it (portable shell):
    while IFS= read -r smiles; do
        args=(-s "$smiles" -c "$smiles" -p "$PDB_DIR" -b "$PROBE_DIR" \
              -o "$OUT_DIR" --num-procs "$NPROCS" --subset-sizes 1 2)
        <your-submit-command> python ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py "${args[@]}"
    done < fragment_work_list.txt
'''

import os
import re
import argparse
import pickle as pkl

from ligand_vdgs.functions.utils import (
    fragment_key_query_mol, fragment_query_mols_equivalent, mol_from_fragment)
from ligand_vdgs.generate_vdgs.estimate_frag_cost import (
    DEFAULT_SAMPLE_SIZE, estimate_fragment_counts)

# Halogen oxyanions (perchlorate, chlorate, periodate, bromate, ...). No
# drug-like ligand carries a hypervalent halogen oxyanion, so any occurrence is
# a crystallization or cryoprotectant salt rather than a real binding moiety.
_SOLVENT_ARTIFACT_HALOGENS = frozenset({17, 35, 53})  # Cl, Br, I


def is_solvent_artifact(mol):
    '''True if the fragment contains a halogen oxyanion center.

    Per-atom test: a Cl/Br/I bearing >= 2 oxygens and no carbon neighbor. The
    carbon-neighbor exclusion keeps genuine organohalogen oxides (``C[I](=O)=O``),
    which the oxygen count alone would discard. An ester of the ion itself
    (``COCl(=O)(=O)=O``) has no carbon on the halogen and is still discarded.

    Detected structurally rather than by string, so alternate drawings of the same
    ion (``O=[Cl](=O)(=O)[O-]``, ``[O-][Cl]([O-])([O-])[O-]``, ...) are all caught.
    Residual limitation: a fragment cut small enough to drop the carbon off an
    organohalogen oxide is indistinguishable from the free ion.
    '''
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in _SOLVENT_ARTIFACT_HALOGENS:
            continue
        neighbor_nums = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
        if neighbor_nums.count(8) >= 2 and 6 not in neighbor_nums:
            return True
    return False


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract fragment SMILES as a one-column work list.")
    parser.add_argument('--frags-dict', default='resources/database_frags_dict.pkl',
                        help="Path to database_frags_dict.pkl. "
                             "Default: resources/database_frags_dict.pkl.")
    parser.add_argument('--output', default='resources/fragment_work_list.txt',
                        help="Path to write the one-column work list. "
                             "Default: resources/fragment_work_list.txt.")
    parser.add_argument('--pdb-dir', required=True,
                        help="Parent PDB database, sampled to count how many CG sites "
                             "each fragment has in it.")
    parser.add_argument('--min-instances', default=250, type=int,
                        help="Min candidate vdG sites (CG occurrences: SMARTS matches "
                             "summed over every ligand copy in the database) for a "
                             "fragment to qualify. Default: 250. See "
                             "make_sge_scripts_for_frags.py for the measured "
                             "fragment-count-vs-threshold table.")
    parser.add_argument('--sample-size', default=DEFAULT_SAMPLE_SIZE, type=int,
                        help=f"Structures sampled for that count. Default: "
                             f"{DEFAULT_SAMPLE_SIZE}.")
    parser.add_argument('--num-procs', default=10, type=int,
                        help="Worker processes for the count. Default: 10.")
    parser.add_argument('--max-size', default=5, type=int,
                        help="Max fragment heavy-atom count. Default: 5.")
    return parser.parse_args()


# Elements the SMILES organic subset lets you write without brackets. Aromatic
# forms are listed separately and in lower case on purpose: fragment keys match
# exactly on aromaticity, and nothing bridges a mismatch (docs/pitfalls.md).
_BARE_ATOMS = frozenset({'B', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I',
                         'b', 'c', 'n', 'o', 'p', 's'})
_BRACKET_ATOM = re.compile(r'\[([^\]]*)\]')
_CHARGE = re.compile(r'(?:\+\d+|-\d+|\++|-+)$')


def charge_normalized_fragment(fragment):
    """The fragment key with formal-charge constraints dropped, or None if none.

    Fragment keys are SMARTS, not SMILES (`utils.mol_from_fragment`), and that is
    what makes this sound rather than merely empirical: an unbracketed organic
    atom in SMARTS carries *no* charge constraint, so `N` matches what `[N+]`,
    `[N-]` and neutral N all match. Dropping the charge therefore yields a
    strictly looser query, and the neutral form's matches are a superset of the
    charged form's by construction. Verified on 13 alias/representative pairs
    from a real build, none of which contradicted it.

    Fragment selection is charge-strict while matching is not, which is why
    `cC(=O)O` and `cC(=O)[O-]` were recorded as two fragments and mined twice.

    Ring-annotated keys keep their annotation: only the charge is dropped, so
    `[O-;!R]` becomes `[O;!R]` and stays a different query from `[O;r6]`.

    This edits the SMARTS text rather than round-tripping through RDKit.
    Sanitizing and re-writing would re-perceive aromaticity and rewrite `c` as
    `C` on any fragment too small to close a ring, silently merging aromatic and
    aliphatic fragments -- a different and much larger claim than "same query,
    one fewer charge constraint". Keys match exactly on aromaticity and nothing
    bridges a mismatch (docs/pitfalls.md).
    """
    def _strip(match):
        inner = match.group(1)
        # A ring-context annotation (`[O-;!R]`, see Frags.ring_query_for_atom)
        # puts the charge in the middle of the bracket, where a `$`-anchored
        # pattern cannot see it. Strip only within the atom part, before the
        # first ';'. Getting this wrong is silent: every key stops normalizing,
        # protonation variants stop collapsing, and the library quietly mines
        # each one twice.
        head, sep, tail = inner.partition(';')
        stripped_head = _CHARGE.sub('', head)
        if stripped_head == head:
            return match.group(0)
        stripped = stripped_head + sep + tail
        return stripped if stripped in _BARE_ATOMS else f'[{stripped}]'

    normalized = _BRACKET_ATOM.sub(_strip, fragment)
    return normalized if normalized != fragment else None


def group_protonation_variants(qualifying):
    """Split ``{fragment SMARTS: ligand-name set}`` into representatives and aliases.

    Charged spellings that share a charge-free form are one group. The group's
    representative is that charge-free SMARTS: it is strictly looser, so it
    matches every member's structures, and the members' jobs would be duplicate
    compute mining nested subsets of the same structures. Two cases:

    * the charge-free form is itself a fragment in the dict (the *neutral twin*,
      e.g. ``CC(=O)O`` for ``CC(=O)[O-]``): it is the representative and the
      charged spellings become aliases;
    * no twin exists but two or more charged spellings share the form (aromatic
      nitro: ``c[N+;!R](=[O;!R])[O;!R]`` and ``c[N+;!R](=[O;!R])[O-;!R]``, which
      RDKit never writes neutral): the charge-free SMARTS is *promoted* to
      representative even though it appears in no CCD drawing. Every member
      becomes an alias; the ligand set is the union. Without this the 2026-09-06
      build mined both nitro spellings, one nested in the other.

    A lone charged spelling with no twin is left as its own representative:
    promoting it would change nothing about what gets mined and would only
    replace a key that exists in the dict with one that does not.

    Forms are compared by ``fragment_query_mols_equivalent``, not string
    equality. Charge-stripping edits the text in place, so it keeps each charged
    spelling's atom order, and both a twin and a sibling spelling are often
    written in a different but equivalent order: on the production fragment
    dict, exact matching collapses 431 groups and misses 60 more (e.g.
    ``[C-;!R][C;r5]([C;r5])[N+;r5]`` normalizes to ``[C;!R][C;r5]([C;r5])[N;r5]``
    while the dict holds ``[C;r5][C;r5]([C;!R])[N;r5]``). A missed group is not
    cosmetic: the variant gets its own job mining a nested subset of the same
    structures, and no alias row, so hit finding can never resolve it.

    Returns ``(representatives, aliases)``, aliases mapping each collapsed
    variant to the representative that now covers it. A promoted representative
    is recognisable as one absent from the input; ``resolve_fragment_key`` and
    ``scripts/lookup_fragment_key.py`` map it back to the dict keys it covers.
    """
    # Atom count is a cheap necessary condition for equivalence, so it keeps the
    # equivalence scan off the whole fragment list for each charged key.
    query_mols = {smiles: fragment_key_query_mol(smiles) for smiles in qualifying}
    by_size = {}
    for smiles, mol in query_mols.items():
        if mol is not None:
            by_size.setdefault(mol.GetNumAtoms(), []).append(smiles)

    def _qualifying_twin(normalized):
        """The qualifying fragment that is the same query as `normalized`, if any."""
        if normalized in qualifying:
            return normalized
        mol = fragment_key_query_mol(normalized)
        if mol is None:
            return None
        for cand in by_size.get(mol.GetNumAtoms(), ()):
            if fragment_query_mols_equivalent(mol, query_mols[cand]):
                return cand
        return None

    # Group charged spellings by their charge-free form, compared as queries so
    # that differently ordered spellings of one form meet. Each group:
    # [normalized SMARTS of the first member seen, its query mol, members].
    representatives, aliases = {}, {}
    charge_groups = []
    for smiles in sorted(qualifying):
        normalized = charge_normalized_fragment(smiles)
        mol = fragment_key_query_mol(normalized) if normalized is not None else None
        if mol is None:
            representatives[smiles] = qualifying[smiles]
            continue
        for group in charge_groups:
            if (group[1].GetNumAtoms() == mol.GetNumAtoms()
                    and fragment_query_mols_equivalent(group[1], mol)):
                group[2].append(smiles)
                break
        else:
            charge_groups.append([normalized, mol, [smiles]])

    for normalized, _, members in charge_groups:
        twin = _qualifying_twin(normalized)
        if twin is not None:
            rep, pooled = twin, [twin] + members
        elif len(members) > 1:
            # Promoted: sorted iteration above makes the spelling deterministic,
            # which matters because it becomes a directory name.
            rep, pooled = normalized, members
        else:
            representatives[members[0]] = qualifying[members[0]]
            continue
        representatives[rep] = set().union(*(qualifying[m] for m in pooled))
        for member in members:
            aliases[member] = rep
    return representatives, aliases


def prepare_fragments(frags_dict, max_size):
    '''(representatives, aliases) for `frags_dict`: the parse + size/solvent filter
    + protonation-variant pooling that select_fragments thresholds on.

    Split out because a run threshold-filters the same dict twice (once at 0 to get
    the counting candidates, once at --min-instances); walking and re-parsing every
    key twice is the dominant cost of both scheduler paths.'''
    qualifying = {}
    unparseable = []
    for _, smiles_dict in frags_dict.items():
        for smiles, lignames in smiles_dict.items():
            mol = mol_from_fragment(smiles)
            if mol is None:
                # Reported unconditionally. A key that will not parse cannot be
                # counted, so there is no threshold to hold it against, and the
                # count it used to be filtered by (CCD ligands) is not the quantity
                # this function thresholds on anyway.
                print(f'[WARNING] skipping unparseable fragment: {smiles}')
                unparseable.append(smiles)
                continue
            if mol.GetNumHeavyAtoms() > max_size:
                continue
            if is_solvent_artifact(mol):
                continue
            qualifying.setdefault(smiles, set()).update(lignames)
    if unparseable:
        print(f'[WARNING] {len(unparseable)} fragment(s) failed to parse and '
              'were excluded.')

    # Pool each variant group's ligands *before* thresholding: the counts come
    # from CCD SMILES drawings, where most acids are drawn neutral, so a split
    # count measures drawing convention rather than prevalence in the PDB.
    return group_protonation_variants(qualifying)


def select_fragments(frags_dict, counts_threshold, max_size, return_aliases=False,
                     instance_counts=None, prepared=None, include=None,
                     include_only=False):
    '''Qualifying fragment SMILES, sorted. Shared by scheduler paths.

    `instance_counts`, when given, replaces the CCD ligand count as the quantity
    thresholded: {fragment: candidate vdG sites in the parent database}, i.e. CG
    occurrences summed over every ligand copy (estimate_frag_cost). That is the
    quantity that decides whether a fragment library is usable, and it diverges
    badly from the ligand count -- against wall time (a proxy for vdG volume),
    log-log correlation is 0.63 for the structure count this is derived from and
    0.28 for the ligand count, and a fragment in 50 CCD ligands took 13.3 h while
    others at the same count finished in minutes. A ligand count measures how many
    *drawings* contain the fragment; this measures how many vdGs it can yield.

    Protonation-state variants are pooled before thresholding and reduced to one
    representative, so a variant group is one job rather than several mining
    nested subsets of the same structures. That is bookkeeping: no vdG's
    contents change, because the charge-free SMARTS already matches both states.

    Chemically distinct forms are *not* pooled -- a donor and an acceptor stay
    split. See docs/pitfalls.md.

    Unparseable keys are reported rather than dropped silently (in
    prepare_fragments), so they cannot be mistaken for fragments that simply missed
    the count threshold.
    '''
    # `prepared` is prepare_fragments' output, to reuse one walk across the two
    # thresholding passes a run makes.
    prepared = (prepared if prepared is not None
                else prepare_fragments(frags_dict, max_size))
    representatives, aliases = prepared

    included, unresolved = resolve_include_fragments(include, prepared)
    if unresolved:
        # Separate the two causes before blaming the spelling. prepare_fragments
        # drops anything over max_size or flagged as a solvent artifact, so a
        # perfectly well-spelled key the caller asked for can be missing simply
        # because this run's --max-size is smaller than the fragment. Telling
        # someone to check their SMARTS in that case sends them after the wrong bug.
        too_big, solventy, absent = [], [], []
        for frag in unresolved:
            mol = mol_from_fragment(frag)
            if mol is None:
                absent.append(frag)
            elif mol.GetNumHeavyAtoms() > max_size:
                too_big.append((frag, mol.GetNumHeavyAtoms()))
            elif is_solvent_artifact(mol):
                solventy.append(frag)
            else:
                absent.append(frag)
        problems = []
        if too_big:
            problems.append(
                "excluded by --max-size {}: {}".format(
                    max_size, ", ".join(f"{f} ({n} heavy atoms)" for f, n in too_big)))
        if solventy:
            problems.append(
                f"excluded as crystallographic solvent artifacts: {solventy}")
        if absent:
            problems.append(
                f"absent from the fragment dict even up to equivalent atom "
                f"ordering, so they have no vdG sites to mine -- check the "
                f"annotated-SMARTS spelling, ring primitives included, with "
                f"scripts/lookup_fragment_key.py: {absent}")
        raise ValueError(
            f"{len(unresolved)} requested fragment(s) could not be selected. "
            + "; ".join(problems))

    if include_only:
        # Top-up build: the threshold pass is skipped entirely, so an existing
        # library gains exactly the requested fragments and nothing else.
        selected = sorted(set(included.values()))
        kept = set(selected)
        aliases = {alias: rep for alias, rep in aliases.items() if rep in kept}
        if return_aliases:
            return selected, aliases
        return selected

    if instance_counts is None:
        selected = sorted(smiles for smiles, lignames in representatives.items()
                          if len(lignames) >= counts_threshold)
    else:
        # A representative covers its aliases' structures too, since its
        # charge-free SMARTS matches both states.
        missing = [s for s in representatives if s not in instance_counts]
        if missing:
            raise ValueError(
                f"instance_counts is missing {len(missing)} fragment(s), e.g. "
                f"{missing[:3]}. It must cover every fragment that survives the "
                "size and solvent filters, not just the ones already selected.")
        selected = sorted(smiles for smiles in representatives
                          if instance_counts[smiles] >= counts_threshold)
    selected = sorted(set(selected) | set(included.values()))
    kept = set(selected)
    aliases = {alias: rep for alias, rep in aliases.items() if rep in kept}
    if return_aliases:
        return selected, aliases
    return selected


def resolve_include_fragments(include, prepared):
    """``(representatives, unresolved)`` for fragments requested by name.

    The count threshold is a sampling estimate, not a property of the fragment:
    across five seeds of the 3000-structure sample, ~40% of the selected set
    changes membership while its size barely moves, so a fragment you actually
    need can miss the cut on the draw you happened to take. This is the way to
    ask for one regardless of its count.

    A requested key is resolved the same way any other key is compared --
    ``fragment_keys_equivalent``, not string equality -- so an equivalent atom
    ordering (`[C;!R][O;!R][C;!R][C;!R]` for `[C;!R][C;!R][O;!R][C;!R]`) still
    finds its fragment. A charged variant resolves to the representative that
    covers it, so asking for one never creates a second job mining a nested
    subset of an existing one.

    Unresolved requests are returned rather than raised on, so the caller can
    report every bad name at once instead of one per run.
    """
    representatives, aliases = prepared
    resolved, unresolved = {}, []
    if not include:
        return resolved, unresolved

    by_size = {}
    for key in list(representatives) + list(aliases):
        mol = fragment_key_query_mol(key)
        if mol is not None:
            by_size.setdefault(mol.GetNumAtoms(), []).append(key)

    for requested in include:
        if requested in representatives:
            resolved[requested] = requested
            continue
        if requested in aliases:
            resolved[requested] = aliases[requested]
            continue
        mol = fragment_key_query_mol(requested)
        if mol is None:
            unresolved.append(requested)
            continue
        match = None
        for cand in by_size.get(mol.GetNumAtoms(), ()):
            if fragment_query_mols_equivalent(mol, fragment_key_query_mol(cand)):
                match = cand
                break
        if match is None:
            unresolved.append(requested)
        else:
            resolved[requested] = aliases.get(match, match)
    return resolved, unresolved


def fragment_dict_keys(frags_dict):
    """Every fragment SMARTS in a ``{elements: {smarts: ligands}}`` frags dict."""
    return {smiles for smiles_dict in frags_dict.values() for smiles in smiles_dict}


def alias_kind(representative, dict_keys):
    """``neutral_twin`` if the representative is a key of the frags dict, else
    ``promoted`` (a charge-stripped SMARTS that appears in no CCD drawing)."""
    return 'neutral_twin' if representative in dict_keys else 'promoted'


def write_fragment_aliases(path, aliases, dict_keys):
    """Write ``alias<TAB>representative<TAB>kind`` rows, one per collapsed variant.

    Consumers that hold a charged fragment name -- saved hit-finder scripts,
    existing analysis paths -- resolve it through this file instead of finding
    no library directory. *dict_keys* (``fragment_dict_keys``) decides ``kind``:
    a ``promoted`` representative exists only as a library key, so looking it
    up in the frags dict finds nothing; ``scripts/lookup_fragment_key.py`` maps
    it back to the dict keys it covers.
    """
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(path, 'w') as f:
        f.write('# alias\trepresentative\tkind\n')
        f.write('# kind=promoted: the representative is the charge-stripped SMARTS of '
                'its aliases and is NOT a key of the fragment dict (no CCD ligand is '
                'drawn that way); use scripts/lookup_fragment_key.py to find the '
                'dict keys it covers.\n')
        for alias in sorted(aliases):
            rep = aliases[alias]
            f.write(f'{alias}\t{rep}\t{alias_kind(rep, dict_keys)}\n')


def resolve_fragment_key(key, candidates):
    """Which of *candidates* the fragment key *key* names, and how.

    Returns ``[(candidate, how), ...]`` with ``how`` one of ``exact`` (same
    string), ``equivalent`` (same annotated query, different atom order) or
    ``charge_variant`` (the candidate's charge-stripped form is equivalent to
    *key*'s, so *key* is the looser query that covers it). This is the search
    a promoted representative needs to be traced back to dict keys, and the one
    a human should run before concluding a key "is not in the dict".
    """
    query = fragment_key_query_mol(key)
    loose_key = charge_normalized_fragment(key) or key
    loose = fragment_key_query_mol(loose_key)
    found = []
    for cand in candidates:
        if cand == key:
            found.append((cand, 'exact'))
            continue
        cand_mol = fragment_key_query_mol(cand)
        if cand_mol is None:
            continue
        if query is not None and fragment_query_mols_equivalent(query, cand_mol):
            found.append((cand, 'equivalent'))
            continue
        if loose is None:
            continue
        cand_loose_key = charge_normalized_fragment(cand)
        cand_loose = (fragment_key_query_mol(cand_loose_key)
                      if cand_loose_key is not None else cand_mol)
        if fragment_query_mols_equivalent(loose, cand_loose):
            found.append((cand, 'charge_variant'))
    order = {'exact': 0, 'equivalent': 1, 'charge_variant': 2}
    return sorted(found, key=lambda pair: (order[pair[1]], pair[0]))


def main():
    args = parse_args()

    if not os.path.isfile(args.frags_dict):
        raise FileNotFoundError(f"--frags-dict path does not exist: {args.frags_dict}")

    with open(args.frags_dict, 'rb') as f:
        frags_dict = pkl.load(f)

    prepared = prepare_fragments(frags_dict, args.max_size)
    candidates = select_fragments(frags_dict, 0, args.max_size, prepared=prepared)
    print(f'Counting {len(candidates)} candidate fragments over '
          f'{args.sample_size} sampled structures...')
    _structure_counts, occurrences = estimate_fragment_counts(
        candidates, args.pdb_dir, sample_size=args.sample_size,
        num_procs=args.num_procs)
    smiles_to_run, aliases = select_fragments(
        frags_dict, args.min_instances, args.max_size, return_aliases=True,
        instance_counts=occurrences, prepared=prepared)

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.output, 'w') as f:
        for smiles in smiles_to_run:
            f.write(f"{smiles}\n")

    alias_path = os.path.splitext(args.output)[0] + '_aliases.tsv'
    dict_keys = fragment_dict_keys(frags_dict)
    write_fragment_aliases(alias_path, aliases, dict_keys)
    if aliases:
        reps = set(aliases.values())
        promoted = sorted(r for r in reps if alias_kind(r, dict_keys) == 'promoted')
        print(f'Collapsed {len(aliases)} protonation variant(s) into '
              f'{len(reps)} representative(s); wrote {alias_path}.')
        if promoted:
            print(f'{len(promoted)} representative(s) are promoted charge-stripped keys '
                  f'absent from the fragment dict (see {alias_path}): {promoted}')

    print(f'Wrote {len(smiles_to_run)} fragments to {args.output}.')


if __name__ == '__main__':
    main()
