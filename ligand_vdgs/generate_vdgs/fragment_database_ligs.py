'''
Generates fragments from all database ligands. Fragments are extracted at every bond radius
from 1 to `--bond-radius` around each atom (not just the max radius, so smaller frags nested
inside larger ones aren't lost), then filtered to `--min-frag-size`-`--max-frag-size` heavy
atoms with >=1 non-carbon atom.

Inputs:
    CCD_smiles: tab-delimited file with 2-3 columns:
                    col 1 (required): SMILES string
                    col 2 (required): short ligand identifier (e.g. CCD 3-letter code or
                                      internal compound ID)
                    col 3 (optional): ligand name; not used by this script

                Defaults are set up for PDB ligands: the CCD file is provided in
                ligand-vdGs under resources/Components-smiles-cactvs.smi, downloaded
                from https://www.wwpdb.org/data/ccd.

Usage (PDB defaults):
    >> python ligand_vdgs/generate_vdgs/fragment_database_ligs.py

Usage (parallel):
    >> python ligand_vdgs/generate_vdgs/fragment_database_ligs.py --nproc 8

Usage (custom ligand list):
    Specify --ccd and --outdir; --logfile defaults to logs/fragment_database_ligs.log.
    >> python ligand_vdgs/generate_vdgs/fragment_database_ligs.py \
           --ccd <path/to/your_ligands> \
           --outdir <output_dir>

Usage (qsub, SGE cluster):
    ligand_vdgs/generate_vdgs/run_fragment_database_ligs.qsub wraps the PDB-default
    invocation with --nproc 8 (-pe smp 8). Submit from the repo root, since the
    script uses -cwd and the wrapped --ccd/--outdir/--logfile paths are relative:
    >> qsub ligand_vdgs/generate_vdgs/run_fragment_database_ligs.qsub
    Edit the qsub script directly to point at a custom --ccd/--outdir.

Output (<outdir>/database_frags_dict.pkl):
    A pickled nested dict with structure:
        {
            elements_str (str): {
                smiles (str): [lig_resname (str), ...]
            }
        }
    where:
        elements_str  -- alphabetically sorted concatenation of non-H element symbols
                         for all atoms in the fragment (e.g. "CCNO"); used to group
                         structural isomers so duplicate-checking only scans isomers.
        smiles        -- annotated SMARTS key of the fragment (see Frags.ring_query_for_atom).
        lig_resname   -- identifier of each database ligand (e.g. CCD 3-letter code for PDB ligands).
    Within each elements_str bucket, fragments are sorted by number of unique ligands
    (descending).
'''

import os
import argparse
import csv
import multiprocessing as mp
import pickle as pkl
from collections import defaultdict
from rdkit import Chem, RDLogger
from ligand_vdgs.functions import Frags
from ligand_vdgs.functions import ligand_perception
from ligand_vdgs.functions.utils import (fragment_key_query_mol,
                                         fragment_query_mols_equivalent)

RDLogger.DisableLog('rdApp.*') 

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ccd', default='resources/Components-smiles-cactvs.smi',
                        help="Path to CCD file containing smiles of molecules. Defaults "
                             "to resources/Components-smiles-cactvs.smi.")
    parser.add_argument('--outdir', default='resources',
                        help="Path to output a pickled dict of database fragments. "
                        "Defaults to ligand-vdGs/resources.")
    parser.add_argument('--logfile', default='logs/fragment_database_ligs.log',
                        help="Path to output log file. Defaults to logs/fragment_database_ligs.log.")
    parser.add_argument('--bond-radius', default=2, type=int,
                        help="Number of bonds away from an atom to cut during "
                        "fragmentation. Defaults to 2.")
    parser.add_argument('--min-frag-size', default=4, type=int,
                        help="Minimum heavy atoms per fragment. Defaults to 4.")
    parser.add_argument('--max-frag-size', default=5, type=int,
                        help="Maximum heavy atoms per fragment. Defaults to 5.")
    parser.add_argument('--nproc', default=1, type=int,
                        help="Worker processes for the per-ligand fragmentation. "
                        "Defaults to 1 (in-process, no multiprocessing at all). "
                        "Results are merged in input order, so the output is "
                        "independent of this value.")
    return parser.parse_args()


# Set once per worker by _init_worker; read by process_ligand. A module global
# rather than a per-task argument so the fragmentation parameters aren't pickled
# and shipped 48k times.
_FRAG_PARAMS = None


def _init_worker(bond_radius, min_frag_size, max_frag_size):
    global _FRAG_PARAMS
    _FRAG_PARAMS = (bond_radius, min_frag_size, max_frag_size)
    RDLogger.DisableLog('rdApp.*')


def process_ligand(task):
    """Fragment one CCD row into plain data. Runs in a worker process.

    Returns only strings and ints -- never a Mol -- so the parallel path ships
    kilobytes per chunk instead of pickling RDKit molecules. `record_frag` needs
    exactly two things per fragment (its element string and its key), and both
    are permutation-invariant, so nothing is lost by discarding the Mol here.

    Warnings are returned rather than printed so the driver can emit them in
    input order; the ones raised inside Frags still print from the worker.
    """
    line_num, row = task
    out = {'stats': [], 'warnings': [], 'failed': None, 'unparseable': None,
           'frags': [], 'empty_site_groups': 0}

    if not (2 <= len(row) <= 3):
        out['warnings'].append(
            f"[WARNING] line {line_num}: skipping malformed row with {len(row)} cols: {row}")
        out['stats'].append('num_ligs_malformed_row')
        return out
    smiles, lig_resname = row[0], row[1]

    # Parse unsanitized first: None then means a genuine SMILES syntax error,
    # and the element filter below only needs atom symbols. Sanitizing before
    # that filter would reclassify metal complexes that fail valence checks
    # as parse failures instead of as non-druglike.
    orig_mol = ligand_perception.perceive_ligand_graph(lig_resname, smiles).mol

    if orig_mol is None:
        # The CCD writes dative bonds as '|', which RDKit rejects (it uses
        # '->'), so metal complexes die here before undesired_elements() ever
        # sees them and land in num_ligs_failed rather than not_druglike.
        # Reclassify by re-parsing without the '|' marks; anything still
        # unparseable stays a reported failure. The '|' test is only a
        # short-circuit, not a correctness guard. A misfire could at most move a
        # row between two counters -- it is discarded either way -- and the
        # SMILES is still logged below, so nothing can vanish silently.
        if '|' in smiles and undesired_elements_in_dative_smiles(smiles):
            out['stats'].append('num_ligs_not_druglike')
            out['unparseable'] = (line_num, lig_resname, smiles)
            return out
        out['warnings'].append(
            f'[WARNING] line {line_num} ({lig_resname}): invalid SMILES: {smiles}')
        out['failed'] = (line_num, lig_resname, 'invalid SMILES', smiles)
        out['stats'].append('num_ligs_failed')
        return out

    if undesired_elements(orig_mol):
        out['stats'].append('num_ligs_not_druglike')
        return out

    # Sanitize the survivors. Everything downstream (MolToSmiles, SMARTS
    # matching, FindAtomEnvironmentOfRadiusN) assumes ring perception and
    # aromaticity have been run, and the canonical SMILES that become
    # fragment identity for the whole pipeline are only stable if they have.
    # ~0.3% of the CCD fails here; drop those rather than mine an
    # un-perceived graph.
    #
    # Only the parent is sanitized -- fragments deliberately inherit its
    # perception rather than being sanitized as standalone molecules.
    # See "Fragments inherit the parent's perception" in docs/pitfalls.md before 
    # changing this.
    try:
        Chem.SanitizeMol(orig_mol)
    except Exception as e:
        out['warnings'].append(
            f'[WARNING] line {line_num} ({lig_resname}): sanitization failed '
            f'({e}): {smiles}')
        out['failed'] = (line_num, lig_resname, f'sanitization: {e}', smiles)
        out['stats'].append('num_ligs_sanitize_failed')
        return out

    results = Frags.manually_remove_Hs(orig_mol, 'single')
    if results is None:
        out['warnings'].append(
            f'[WARNING] line {line_num} ({lig_resname}): failed to process SMILES '
            f'(H-removal/canonicalization): {smiles}')
        out['failed'] = (line_num, lig_resname, 'H-removal/canonicalization', smiles)
        out['stats'].append('num_ligs_failed')
        return out
    mol_info, _ = results
    mol, _ = mol_info
    # Make sure there's at least a carbon (rule out ions, buffer components
    # like sulfate/phosphate, Fe-S clusters, etc.)
    if not Frags.is_organic(mol):
        out['stats'].append('num_ligs_not_druglike')
        return out
    out['stats'].append('num_ligs_parsed')

    # Decompose the ligand into fragments and store the fragment keys. Keys are
    # ANNOTATED SMARTS, not SMILES: each atom carries a ring primitive read from
    # the parent before PathToSubmol cuts the ring (r<n> for saturated ring atoms,
    # !R for acyclic, bare for aromatic), so THF and THP stay distinct while
    # pyridine and pyrrole pool. See Frags.ring_query_for_atom. Aromatic vs.
    # aliphatic is still carried by case (C vs. c), as it is in SMILES.
    # Fragment on bond radii `bond_radius` AND the positive
    # integers less than `bond_radius`, because for example, drugs containing
    # sulfonamide might produce only 6-atom sulfonamides and not CS(N)(=O)=O.
    bond_radius, min_frag_size, max_frag_size = _FRAG_PARAMS
    filtered_frags = Frags.get_fragments(bond_radius, mol, min_frag_size, max_frag_size)
    for sub_smiles, substruct_site_groups in filtered_frags.items():
        # get_fragments should never return an empty site group or an empty
        # perm list, but index defensively: the pickle is only written after
        # the whole file, so an IndexError here would discard the entire run.
        if not substruct_site_groups or not substruct_site_groups[0]:
            out['empty_site_groups'] += 1
            continue
        substruct_site = substruct_site_groups[0] # only need one site group
        substruct_perm = substruct_site[0] # only need one perm.
        sub = substruct_perm[0] # select Mol obj (Mol obj, perm inds, orig mol inds)
        out['frags'].append((fragment_elements(sub), sub_smiles))
    out['has_frags'] = bool(filtered_frags)
    out['lig_resname'] = lig_resname
    return out


def fragment_elements(substruct):
    """Alphabetically-sorted element string used to bucket a fragment.

    Skips hydrogens that are real atoms in the graph. Fragments arrive
    H-stripped (Frags.manually_remove_Hs -> RemoveAllHs), so this should never
    fire; implicit hydrogens never appear in GetAtoms() and are unaffected
    either way. Permutation-invariant, which is what lets the worker compute it
    and discard the Mol.
    """
    return "".join(sorted([a.GetSymbol() for a in substruct.GetAtoms() if
                           a.GetSymbol() != 'H']))


def main():
    args = parse_args()

    if not os.path.isfile(args.ccd):
        raise FileNotFoundError(f"CCD file {args.ccd} does not exist.")

    out_dict_path = os.path.join(args.outdir, 'database_frags_dict.pkl')
    if os.path.exists(out_dict_path):
        raise FileExistsError(f'Output file {out_dict_path} already exists. Exiting.')

    # Warned, not raised: unlike the miner, this script has a real fallback
    # (the ligand table's SMILES), and non-CCD ligands legitimately use it.
    # But a run where the store is simply missing produces a whole library
    # built on SMILES chemistry that will not agree with the miner's, and the
    # only other trace is the provenance counter at the end.
    try:
        ligand_perception.require_template_store()
    except Exception as exc:
        print(f'[WARNING] CCD template store unusable ({exc}); every ligand '
              'will fall back to its SMILES, which is NOT the chemistry the '
              'miner will use.', flush=True)

    os.makedirs(args.outdir, exist_ok=True)
    log_dir = os.path.dirname(args.logfile)
    if log_dir:
        os.makedirs(log_dir, exist_ok=True)

    # frag_dict: {elements_str: {smiles: [lig_resnames]}}
    # Keyed by alphabetically-sorted element string to group isomers together; duplicate-
    # checking via HasSubstructMatch() is expensive for large numbers of fragments, so we
    # only run it within isomer groups.
    frag_dict = {}
    stats = defaultdict(int)
    failed_smiles = []
    # Rows RDKit refused to parse that carry an undesired element anyway; see the '|'
    # branch below. Reported separately from failed_smiles.
    unparseable_undesired = []
    # Parsed-query cache for record_frag's duplicate check, keyed by fragment SMILES.
    # Kept out of frag_dict, which gets pickled and must stay Mol-free.
    smarts_cache = {}
    # (elements, fragment key) -> the already-recorded key it was found equivalent
    # to, so a degenerate key rescans its bucket once per run instead of once per
    # ligand carrying it. Also Mol-free and also kept out of frag_dict.
    alias_cache = {}
    # {exception message: [count, (smiles, existing_smiles)]}, written by report_stats.
    comparison_errors = {}

    def rows():
        with open(args.ccd, mode="r", newline="") as file:
            reader = csv.reader(file, delimiter="\t")
            for line_num, row in enumerate(reader, start=1):
                yield line_num, row

    frag_params = (args.bond_radius, args.min_frag_size, args.max_frag_size)
    if args.nproc > 1:
        # imap, not imap_unordered: results must be merged in input order.
        # record_frag resolves a degenerate key to the *first* equivalent key
        # already in its bucket, so bucket insertion order decides which spelling
        # of a fragment is the recorded one. Unordered merging would make the
        # output depend on worker scheduling.
        pool = mp.Pool(args.nproc, initializer=_init_worker, initargs=frag_params)
        results = pool.imap(process_ligand, rows(), chunksize=64)
    else:
        # No Pool at all at --nproc 1, so the default path keeps its current
        # single-process semantics (and its tracebacks).
        pool = None
        _init_worker(*frag_params)
        results = (process_ligand(task) for task in rows())

    try:
        for out in results:
            stats['num_ligs_total'] += 1
            for warning in out['warnings']:
                print(warning)
            for key in out['stats']:
                stats[key] += 1
            if out['failed'] is not None:
                failed_smiles.append(out['failed'])
            if out['unparseable'] is not None:
                unparseable_undesired.append(out['unparseable'])
            stats['num_empty_site_groups'] += out['empty_site_groups']
            if not out['frags'] and not out.get('has_frags'):
                continue
            lig_resname = out['lig_resname']
            for elements, sub_smiles in out['frags']:
                record_frag(elements, frag_dict, sub_smiles, lig_resname,
                            smarts_cache, comparison_errors, alias_cache)
            if out.get('has_frags'):
                stats['num_ligs_with_frags'] += 1
    finally:
        if pool is not None:
            pool.terminate()
            pool.join()

    # Once per key that would not parse, not once per comparison against it: the
    # bucket scan revisits the same bad key for every candidate in its bucket, which
    # inflated this figure by the traversal count.
    stats['num_smarts_parse_failed'] = sum(
        1 for mol in smarts_cache.values() if mol is None)

    final_num_frags = sum(len(frags) for frags in frag_dict.values())
    # Written before the exit check, so the log explaining an empty run survives it.
    report_stats(frag_dict, stats, failed_smiles, unparseable_undesired,
                 comparison_errors, args.logfile)
    # An empty result is a bad --ccd, not a fragment-free ligand set: a wrong
    # delimiter, a header row, or a non-SMILES first column makes every row
    # "malformed", which only warns. Without this the script writes a valid,
    # empty pickle and exits 0, and the whole downstream build then quietly
    # produces nothing.
    if not stats['num_ligs_parsed'] or not final_num_frags:
        raise SystemExit(
            f'ERROR: parsed {stats["num_ligs_parsed"]} of {stats["num_ligs_total"]} '
            f'input molecules and recorded {final_num_frags} fragments; refusing to '
            f'write an empty {out_dict_path}. Check that {args.ccd} is tab-delimited '
            f'with SMILES in column 1 and an identifier in column 2. '
            f'See {args.logfile} for the per-row reasons.')
    output_results(frag_dict, out_dict_path)

def output_results(frag_dict, out_pkl):
    with open(out_pkl, 'wb') as f:
        pkl.dump(sort_frag_dict(frag_dict), f)


def sort_frag_dict(frag_dict):
    sorted_dict = {}
    for elements, fragments in frag_dict.items():
        sorted_dict[elements] = dict(
            sorted(fragments.items(), key=lambda x: len(set(x[1])), reverse=True)
        )
    return sorted_dict

def report_stats(frag_dict, stats, failed_smiles, unparseable_undesired,
                 comparison_errors, logfile):
    final_num_frags = sum(len(frags) for frags in frag_dict.values())
    stats['num_comparison_errors'] = sum(c for c, _ in comparison_errors.values())
    with open(logfile, 'w') as log:
        log.write(f'{"="*30} Report {"="*30}\n'
                  f'Number of mols in input file: {stats["num_ligs_total"]}\n'
                  f'Number of unique fragments recorded: {final_num_frags}\n'
                  f'Number of mols parsed: {stats["num_ligs_parsed"]}\n'
                  f'Number of mols that produced >=1 fragment: {stats["num_ligs_with_frags"]}\n'
                  f'Number of mols unsuccessfully parsed: {stats["num_ligs_failed"]}\n'
                  f'Number of mols that failed sanitization: {stats["num_ligs_sanitize_failed"]}\n'
                  f'Number of mols skipped b/c not druglike: {stats["num_ligs_not_druglike"]}\n'
                  f'  ...of which were unparseable (bad SMILES syntax): '
                  f'{len(unparseable_undesired)}\n'
                  f'Number of malformed input rows: {stats["num_ligs_malformed_row"]}\n'
                  f'Number of SMARTS re-parse failures: {stats["num_smarts_parse_failed"]}\n'
                  f'Number of substructure comparison errors: {stats["num_comparison_errors"]}\n'
                  f'Number of empty fragment site groups: {stats["num_empty_site_groups"]}\n')
        if failed_smiles:
            log.write('\nFailed SMILES:\n')
            for line_num, lig_resname, reason, smi in failed_smiles:
                log.write(f'  line {line_num} ({lig_resname}) [{reason}]: {smi}\n')
        if unparseable_undesired:
            log.write('\nUnparseable rows carrying an undesired element (counted as '
                      'not druglike; typically metal complexes written with the CCD '
                      'dative-bond "|" notation, which RDKit rejects):\n')
            for line_num, lig_resname, smi in unparseable_undesired:
                log.write(f'  line {line_num} ({lig_resname}): {smi}\n')
        if comparison_errors:
            log.write('\nSubstructure comparison errors (count, example pair):\n')
            for reason, (count, (smi, existing_smi)) in sorted(
                    comparison_errors.items(), key=lambda kv: kv[1][0], reverse=True):
                log.write(f'  {count}x {reason}\n'
                          f'      e.g. {smi} vs {existing_smi}\n')

def record_frag(elements, frag_dict, smiles, lig_resname, smarts_cache,
                comparison_errors, alias_cache):

    '''
    `elements` is the fragment's element string from `fragment_elements` (computed
    by the caller, which may be a worker process that has already discarded the Mol).

    Update frag_dict to store these substructs and the molecule from which they are from.
    Dict Key: Alphabetically-sorted concatenation of element symbols in the substruct
              (elements_str; used to group structural isomers together)
    Dict Value: Subdicts where: 
        Subdict Key: annotated SMARTS key of fragments with the same composition as
                     the Dict Key
        Subdict Value: list of CCD resnames of ligands that contain that fragment.

    The Dict key is essentially a representation of its formula, and is used because when
    each new molecule frag is being added, it can be expensive to check whether the frag is
    already represented in the dict, as opposed to checking only its isomers.

    `smarts_cache` maps fragment key -> the isotope-encoded query Mol that key
    comparison actually consumes (or None if it wouldn't parse), so each key is
    parsed once per run rather than once per pair in this O(bucket^2) scan. It is
    passed in rather than stored on frag_dict, which gets pickled and must stay
    Mol-free. Keys that failed to parse are the None entries; main counts them once
    each after the run rather than once per comparison.

    `alias_cache` maps (elements, fragment key) -> the already-recorded key that key
    was found equivalent to. The scan below returns on the first match in insertion
    order and keys are only ever added to a bucket, never removed or reordered, so
    that answer is fixed once computed. Without it a degenerate key is never itself
    added as a bucket key, so every later ligand carrying it rescans the whole
    bucket -- the dominant cost of this module at ~5.7k fragments. One consequence:
    comparison errors raised on the pairs preceding a match are tallied once per key
    rather than once per occurrence of it.

    `comparison_errors` accumulates {exception message: [count, (smiles, existing_smiles)]}
    across the run, for report_stats to write out. Kept out of stdout: this loop runs
    once per bucket entry per fragment per ligand, so printing each occurrence would
    emit millions of lines.

    Note: All-carbon fragments are excluded upstream by Frags.get_fragments() before
          reaching this function.
    '''

    # Determine whether to add fragment to the dict
    if elements not in frag_dict.keys(): # automatically a new entry
        frag_dict[elements] = {smiles: [lig_resname]}
        return

    if smiles in frag_dict[elements].keys(): # Is the smiles already recorded?
        if lig_resname not in frag_dict[elements][smiles]:
            frag_dict[elements][smiles].append(lig_resname)
        return

    # Recorded under a degenerate key an earlier occurrence already resolved.
    alias = alias_cache.get((elements, smiles))
    if alias is not None:
        if lig_resname not in frag_dict[elements][alias]:
            frag_dict[elements][alias].append(lig_resname)
        return

    # Even if the smiles isn't already recorded, it still might be represented b/c the
    # smiles could be degenerate, even w/ canonical. Use HasSubstructMatch() to be sure.
    if smiles not in smarts_cache:
        smarts_cache[smiles] = fragment_key_query_mol(smiles)
    sub_mol = smarts_cache[smiles]

    # An unparseable key can never be equivalent to anything, so skip the scan.
    for existing_smiles in ([] if sub_mol is None else frag_dict[elements].keys()):
        if existing_smiles not in smarts_cache:
            smarts_cache[existing_smiles] = fragment_key_query_mol(existing_smiles)
        existing_sub_mol = smarts_cache[existing_smiles]
        if existing_sub_mol is None:
            continue
        try:
            # Key against key, not key against `substruct`. The submol has no
            # RingInfo (PathToSubmol output is never sanitized), so every
            # annotated key raised here and the bare `except` below recorded the
            # duplicate as a distinct fragment -- 190 of 1029 selected fragments
            # were redundant that way. It also has its rings cut, so an `r5` key
            # would not match its own fragment even with RingInfo supplied.
            if fragment_query_mols_equivalent(existing_sub_mol, sub_mol): # duplicate,
                                                       # but still add bc diff lig name
                alias_cache[(elements, smiles)] = existing_smiles
                if lig_resname not in frag_dict[elements][existing_smiles]:
                    frag_dict[elements][existing_smiles].append(lig_resname)
                return
        except Exception as e:
            # Can't reliably compare. Treating as non-duplicate may cause false positives
            # (duplicate recorded as distinct), but doesn't lose data. Tally by message
            # so report_stats can write every reason with its count and an example.
            entry = comparison_errors.setdefault(str(e), [0, (smiles, existing_smiles)])
            entry[0] += 1
            continue

    # Passed all checks and deemed nonredundant.
    frag_dict[elements][smiles] = [lig_resname]

UNDESIRED_ELEMENTS = {
    # metals
    'Be', 'Ba', 'Bi', 'Pt', 'Ru', 'Ir', 'Fe', 'Zn', 'Mg', 'Cu', 'Sn', 'Zr', 'Pb', 'Ga',
    'Hg', 'Pd', 'Si', 'Mo', 'W', 'Se', 'Mn', 'Ti', 'Y', 'V', 'Ni', 'Rh', 'Te', 'Au', 'Ag',
    'Co', 'Sb', 'Re', 'As', 'Cd', 'Hf', 'Na', 'Li', 'Ca', 'Os', 'Cr', 'In', 'Al', 'U', 'K',
    # lanthanides
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    # noble gases
    'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn',
}
# Boron is deliberately absent: boronic acids and boronate esters are real covalent
# warheads.

def undesired_elements(mol):
    '''Checks parsed atom symbols rather than regexing the raw SMILES string: a raw-string
    match (e.g. for "In" or "Sn") can false-positive on ordinary organic SMILES where an
    unbracketed atom happens to be immediately followed by another that completes one of
    these 2-letter codes (e.g. "In1cccc1" is iodine bonded to an aromatic ring N, not
    indium).'''
    return any(a.GetSymbol() in UNDESIRED_ELEMENTS for a in mol.GetAtoms())

def undesired_elements_in_dative_smiles(smiles):
    '''Element check for rows RDKit refuses to parse because of the CCD's '|' dative-bond
    notation. Drops the '|' marks and re-parses, so RDKit still does the element parsing
    and undesired_elements() still does the deciding -- removing bond marks cannot change
    the atom set, which is all this needs. Anything that remains unparseable returns False
    and stays a reported parse failure, so garbage never gets silently reclassified.'''
    mol = Chem.MolFromSmiles(smiles.replace('|', ''), sanitize=False)
    return mol is not None and undesired_elements(mol)

if __name__ == "__main__":
    main()
