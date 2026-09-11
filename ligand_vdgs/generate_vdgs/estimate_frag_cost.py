"""Estimate per-fragment mining cost by sampling the parent database.

One SMARTS pass over a random sample, scaled by the sampling fraction, matched against
every fragment per structure (one DB read, not one per fragment). Yields
structure count (sizes the job) and CG occurrence count (decides whether the
fragment is worth mining).

Usage:
    python ligand_vdgs/generate_vdgs/estimate_frag_cost.py \\
        --pdb-dir <path/to/parent_database/> --output resources/frag_cost_estimate.tsv
"""
import math
import os
import re
import sys
import glob
import random
import argparse
import multiprocessing as mp
import concurrent.futures
# Imported, not attribute-accessed: on 3.10 concurrent.futures exposes only public
# class names via __getattr__, so `concurrent.futures.process.X` raises
# AttributeError exactly when a worker has died.
from concurrent.futures.process import BrokenProcessPool
import numpy as np
from ligand_vdgs.functions.db_identity import identity_of
from ligand_vdgs.functions.utils import file_sha256
from ligand_vdgs.functions import parent_db

HERE = os.path.dirname(os.path.abspath(__file__))


def add_vdg_miner_paths():
    # The vdG-miner submodule is not a package, so it still needs sys.path.
    path = os.path.join(HERE, "..", "..", "external", "vdG-miner", "vdg_miner", "vdg")
    if path not in sys.path:
        sys.path.append(path)


# --min-instances thresholds the OCCURRENCE count, which is overdispersed (one
# structure carries many sites/NCS copies): observed SD is 2.6-9.5x the naive
# binomial estimate over five seeds. At n=3000, ~41% of the selected-fragment set
# is seed-dependent (resolve_include_fragments/--include covers this); raise
# --sample-size if the boundary matters more than runtime.
DEFAULT_SAMPLE_SIZE = 3000

# Above this share of the sample lost to unreadable or erroring structures, the
# extrapolation is reported as a failure rather than scaled up from what is left.
MAX_LOST_SAMPLE_FRACTION = 0.1


def sample_pdb_paths(pdb_dir, sample_size, seed=0):
    """A deterministic random sample of the parent DB, or all of it if smaller."""
    paths = [path for _stem, path in parent_db.iter_structures(pdb_dir)]
    if not paths:
        raise ValueError(f"No PDBs found under {pdb_dir}")
    if sample_size is None or sample_size >= len(paths):
        return paths, len(paths)
    return random.Random(seed).sample(paths, sample_size), len(paths)


_WORKER_FRAGMENTS = None
_WORKER_RING_CONSTRAINTS = None
_WORKER_PREFILTER = None

# Element symbols in the order cg.ring_size_constraints scans them: two-letter
# first, so 'Cl' is not read as 'C' followed by a ring-closure digit.
_ORGANIC_SUBSET = ('Cl', 'Br', 'B', 'C', 'N', 'O', 'P', 'S', 'F', 'I',
                   'b', 'c', 'n', 'o', 'p', 's')

# Atom classes counted by the prefilter; every key atom and every OBMol atom lands
# in exactly one, making a per-class count a necessary condition for a match.
_CLS_AROMATIC, _CLS_ACYCLIC, _CLS_RING, _CLS_ANY = 0, 1, 2, 3
_N_CLS = 4


def _parse_pattern_atoms(smarts):
    """[(element, is_aromatic, primitives)] in pattern-atom order, or None if the
    key uses grammar Frags doesn't emit -- disables the prefilter for that pattern
    rather than risking a wrong skip."""
    atoms, i, n = [], 0, len(smarts)
    while i < n:
        char = smarts[i]
        if char == '[':
            depth, j = 1, i + 1
            while j < n and depth:
                if smarts[j] == '[':
                    depth += 1
                elif smarts[j] == ']':
                    depth -= 1
                j += 1
            primitives = smarts[i + 1:j - 1]
            # `[c,n]` (disjunction) requires no specific element, so no per-element
            # count is valid; Frags doesn't emit these but a future source might.
            if ',' in primitives:
                return None
            match = re.match(r'([A-Z][a-z]?|[a-z])', primitives)
            if not match:
                return None
            symbol = match.group(1)
            atoms.append((symbol.capitalize(), symbol.islower(), primitives))
            i = j
            continue
        if char == '*':
            return None
        for symbol in _ORGANIC_SUBSET:
            if smarts.startswith(symbol, i):
                atoms.append((symbol.capitalize(), symbol.islower(), ''))
                i += len(symbol)
                break
        else:
            # Bond, branch, ring-closure digit, %nn, or dot: not an atom.
            i += 2 if char == '%' else 1
    return atoms


def _build_prefilter(fragments, ring_constraints):
    """(requirement matrix, element column index, ring-size column index).

    94% of pattern-vs-ligand calls in a real run cannot match (no sulfur, no
    5-ring, too few aromatic carbons), so each row holds the atom count a key
    needs per (element, class) and per ring size; a ligand below any of them is
    skipped without calling OpenBabel. Necessary, not sufficient -- surviving
    matches are unchanged. Unclassifiable keys get an all-zero row (no prefilter).
    """
    parsed = [_parse_pattern_atoms(f) for f in fragments]
    symbols = sorted({e for atoms in parsed if atoms for e, _, _ in atoms})
    sym_col = {s: i * _N_CLS for i, s in enumerate(symbols)}
    ring_offset = _N_CLS * len(symbols)
    sizes = sorted({size for constraints in ring_constraints
                    for size in constraints if size})
    ring_col = {size: ring_offset + i for i, size in enumerate(sizes)}
    req = np.zeros((len(fragments), ring_offset + len(sizes)), dtype=np.int16)
    for row, (atoms, constraints) in enumerate(zip(parsed, ring_constraints)):
        if atoms is None or len(atoms) != len(constraints):
            continue
        for (element, aromatic, primitives), ring_size in zip(atoms, constraints):
            base = sym_col[element]
            req[row, base + _CLS_ANY] += 1
            if aromatic:
                req[row, base + _CLS_AROMATIC] += 1
            elif ring_size is not None:
                req[row, base + _CLS_RING] += 1
                req[row, ring_col[ring_size]] += 1
            elif '!R' in primitives:
                req[row, base + _CLS_ACYCLIC] += 1
            # Otherwise the atom constrains only its element, already counted.
    return req, sym_col, ring_col


def _prefilter_candidates(mol):
    """Indices of fragments whose atom counts this molecule could supply."""
    from openbabel import openbabel as ob
    req, sym_col, ring_col = _WORKER_PREFILTER
    counts = np.zeros(req.shape[1], dtype=np.int16)
    for atom in ob.OBMolAtomIter(mol):
        base = sym_col.get(ob.GetSymbol(atom.GetAtomicNum()))
        if base is None:
            continue  # No fragment needs this element.
        counts[base + _CLS_ANY] += 1
        if atom.IsAromatic():
            counts[base + _CLS_AROMATIC] += 1
        elif not atom.IsInRing():
            counts[base + _CLS_ACYCLIC] += 1
        else:
            counts[base + _CLS_RING] += 1
            column = ring_col.get(atom.MemberOfRingSize())
            if column is not None:
                counts[column] += 1
    return np.flatnonzero((req <= counts).all(axis=1))


def _init_worker(fragments):
    """Compile the patterns once per worker, not once per structure."""
    global _WORKER_FRAGMENTS, _WORKER_RING_CONSTRAINTS, _WORKER_PREFILTER
    add_vdg_miner_paths()
    from cg import (check_ring_constraints, compile_smarts_patterns,
                    ring_size_constraints)
    _WORKER_FRAGMENTS = compile_smarts_patterns(fragments)
    # Applies the same r<n> re-filter find_cg_matches uses; without it
    # ring-constrained fragments are over-counted and tiers drift from reality.
    _WORKER_RING_CONSTRAINTS = [ring_size_constraints(f) for f in fragments]
    for smarts, pattern, constraints in zip(fragments, _WORKER_FRAGMENTS,
                                            _WORKER_RING_CONSTRAINTS):
        check_ring_constraints(smarts, pattern, constraints)
    _WORKER_PREFILTER = _build_prefilter(fragments, _WORKER_RING_CONSTRAINTS)


def _count_one(pdb_path):
    """({fragment index: CG occurrences}, unparsable blocks, file was unreadable).

    Occurrences, not presence: a fragment hit at 3 sites (or 3 NCS copies) is 3x
    the cost. `GetUMapList` dedupes only the pattern's own symmetry, so a symmetric
    CG counts once per site. Not cg.count_matching_structures, which stops at the
    first hit per pattern.
    """
    from cg import (match_satisfies_ring_sizes, read_ligand_blocks,
                    suppress_stdout_stderr)
    from ligand_vdgs.functions.ligand_perception import perceive_ligand_instance
    ligands = read_ligand_blocks(pdb_path)
    if ligands is None:
        # Unreadable after retries -- distinct from "no ligands", which would stay
        # in the denominator and bias every count downward (the direction that
        # gets jobs killed).
        return {}, 0, True
    if not ligands:
        return {}, 0, False
    counts = {}
    read_failures = 0
    with suppress_stdout_stderr():
        # Uses the miner's perception (deletes hydrogens; without it every
        # `D<n>`-annotated key counts zero). None means unreadable -> tallied as a
        # read failure, not zero occurrences, so the fragment isn't undercounted.
        for key, block in ligands.items():
            perceived = perceive_ligand_instance(block, key[-1])
            if perceived is None:
                read_failures += 1
                continue
            mol = perceived.obmol
            for i in _prefilter_candidates(mol):
                pattern = _WORKER_FRAGMENTS[i]
                if pattern.Match(mol):
                    constraints = _WORKER_RING_CONSTRAINTS[i]
                    n = sum(1 for m in pattern.GetUMapList()
                            if match_satisfies_ring_sizes(mol, m, constraints))
                    if n:
                        counts[int(i)] = counts.get(int(i), 0) + n
    return counts, read_failures, False


def _consume_result(pdb_path, get_result, tally, errored):
    """Tally one structure's result, recording a failure instead of raising.

    One bad structure shouldn't abort a multi-minute pass. Failures go to
    `errored` so the caller can drop them from the extrapolation's denominator
    (keeping them would bias every count downward). Shared by both branches so the
    untested single-process path can't diverge from production.

    BrokenProcessPool means the pool is dead, not this structure -- propagates.
    """
    try:
        tally(get_result())
    except BrokenProcessPool:
        raise
    except Exception as exc:
        errored.append((pdb_path, repr(exc)))


def estimate_fragment_counts(fragments, pdb_dir, sample_size=DEFAULT_SAMPLE_SIZE,
                             num_procs=1, seed=0):
    """({fragment: structures containing it}, {fragment: CG occurrences}).

    Same pass, extrapolated to the parent databse. Structure count sizes the job (tiers
    are calibrated against it); occurrence count decides whether the fragment is
    worth mining. A fragment missed entirely by the sample gets 0 in both.
    """
    # Cleared first so an early return below doesn't leave a stale scale in place.
    _LAST_SAMPLE_SCALE[0] = None
    fragments = list(fragments)
    if not fragments:
        return {}, {}
    # Guarded here too (not just parse_args): other scripts import this directly.
    if sample_size is not None and sample_size < 1:
        raise ValueError(f"sample_size must be a positive integer, got {sample_size}")
    if num_procs < 1:
        raise ValueError(f"num_procs must be at least 1, got {num_procs}")
    sample, total = sample_pdb_paths(pdb_dir, sample_size, seed=seed)
    struct_hits = [0] * len(fragments)
    occurrences = [0] * len(fragments)
    add_vdg_miner_paths()
    # Compiled once in the parent to surface a bad fragment clearly; the same
    # failure inside a worker initializer surfaces as BrokenProcessPool instead.
    from cg import check_ring_constraints, compile_smarts_patterns, ring_size_constraints
    patterns = compile_smarts_patterns(fragments)
    for smarts, pattern in zip(fragments, patterns):
        check_ring_constraints(smarts, pattern, ring_size_constraints(smarts))

    failures = [0]
    unreadable = [0]
    errored = []

    def _tally(result):
        matched, n_failed, was_unreadable = result
        failures[0] += n_failed
        unreadable[0] += bool(was_unreadable)
        for i, n in matched.items():
            struct_hits[i] += 1
            occurrences[i] += n

    if num_procs > 1:
        ctx = mp.get_context("spawn")
        # ProcessPoolExecutor, not mp.Pool: mp.Pool silently respawns a worker the
        # OOM killer takes and never reports the lost task (bpo-22393); the
        # executor marks the pool broken and raises instead.
        pool = concurrent.futures.ProcessPoolExecutor(
            max_workers=num_procs, mp_context=ctx,
            initializer=_init_worker, initargs=(fragments,))
        try:
            futures = {pool.submit(_count_one, p): p for p in sample}
            for fut in concurrent.futures.as_completed(futures):
                _consume_result(futures[fut], fut.result, _tally, errored)
        except BrokenProcessPool as exc:
            raise RuntimeError(
                f'A counting worker died without returning a result ({exc}). The '
                f'usual cause is the OOM killer; re-run with fewer --num-procs or '
                f'a larger -l mem_free (which is PER SLOT under -pe smp).') from exc
        finally:
            # cancel_futures: without it, __exit__ drains every queued structure
            # before a systematic failure surfaces.
            pool.shutdown(wait=True, cancel_futures=True)
    else:
        _init_worker(fragments)
        for pdb_path in sample:
            _consume_result(pdb_path, lambda p=pdb_path: _count_one(p),
                            _tally, errored)
    # Loud: these depress every count uniformly, so a large number means
    # selection and tiers are built on fewer ligands than the sample size suggests.
    if failures[0]:
        print(f'[WARNING] OpenBabel failed to read {failures[0]} ligand block(s) '
              f'in the {len(sample)}-structure sample; their fragments are '
              f'undercounted.', file=sys.stderr)
    # Must also leave the denominator, or counts get scaled against structures
    # nothing was counted from.
    if errored:
        print(f'[WARNING] {len(errored)} of {len(sample)} sampled structures raised '
              f'while being counted and were excluded; e.g. '
              f'{errored[0][0]}: {errored[0][1]}', file=sys.stderr)
    read_ok = len(sample) - unreadable[0] - len(errored)
    if unreadable[0]:
        print(f'[WARNING] {unreadable[0]} of {len(sample)} sampled structures were '
              f'unreadable; extrapolating from the {read_ok} that were read.',
              file=sys.stderr)
    if read_ok == 0:
        raise RuntimeError(
            f'none of the {len(sample)} sampled structures under {pdb_dir} could be '
            'read; refusing to emit an all-zero estimate.')
    # A handful of failures is normal; losing a large share means the estimate no
    # longer samples what it claims to.
    lost = len(sample) - read_ok
    if lost > MAX_LOST_SAMPLE_FRACTION * len(sample):
        raise RuntimeError(
            f'{lost} of {len(sample)} sampled structures under {pdb_dir} could not be '
            f'counted ({100.0 * lost / len(sample):.0f}%, limit '
            f'{100 * MAX_LOST_SAMPLE_FRACTION:.0f}%); refusing to emit an estimate '
            f'extrapolated from the remainder.')
    scale = total / read_ok
    _LAST_SAMPLE_SCALE[0] = scale
    return ({frag: int(round(n * scale)) for frag, n in zip(fragments, struct_hits)},
            {frag: int(round(n * scale)) for frag, n in zip(fragments, occurrences)})


# Set by estimate_fragment_counts; main() records it in the TSV header (counts are
# already scaled).
_LAST_SAMPLE_SCALE = [None]


def sampling_upper_bound(est_count, scale, z=2.0):
    """Upper confidence bound on an extrapolated count, for sizing a job.

    TSV counts are point estimates from one draw; across seeds, 15.5% of selected
    fragments changed resource tier, some into a queue too short to finish (killed
    at h_rt, symptom: no 'Job completed.' line). Loss is asymmetric
    (under-requesting kills the job; over-requesting only costs queue latency), so
    size on the upper end: `est_count / scale` recovers the raw sampled count
    (relative SE ~1/sqrt(raw), Poisson), inflated by `z` of those. Raw 0 gets no
    bound -- make_sge_scripts_for_frags' --include tier floor handles those.
    """
    if not est_count or not scale or scale <= 0:
        return est_count
    raw = est_count / scale
    if raw <= 0:
        return est_count
    return int(round(est_count * (1.0 + z / math.sqrt(raw))))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--frags-dict', default='resources/database_frags_dict.pkl',
                        help="Path to database_frags_dict.pkl.")
    parser.add_argument('--max-size', default=5, type=int,
                        help="Max fragment heavy-atom count. Default: 5.")
    parser.add_argument('--pdb-dir', required=True, help="Parent database.")
    parser.add_argument('--sample-size', default=DEFAULT_SAMPLE_SIZE, type=int,
                        help=f"Structures to sample. Default: {DEFAULT_SAMPLE_SIZE}.")
    parser.add_argument('--num-procs', default=10, type=int, help="Worker processes.")
    parser.add_argument('--seed', default=0, type=int, help="Sampling seed.")
    parser.add_argument('--output', default='resources/frag_cost_estimate.tsv',
                        help="Three-column TSV: fragment, estimated structure "
                             "count, estimated CG occurrences.")
    args = parser.parse_args()
    # Checked early: read late (after DB glob, SMARTS compile), a non-positive
    # sample divides by zero and num_procs=0 raises from inside the pool.
    if args.sample_size < 1:
        parser.error("--sample-size must be a positive number of structures.")
    if args.num_procs < 1:
        parser.error("--num-procs must be at least 1.")
    return args


def read_estimate_header(path):
    """The `# key\tvalue` provenance lines main() writes, as a dict.

    Empty for an older TSV without them -- callers must treat a missing key as
    "unknown", not "mismatched".
    """
    header = {}
    with open(path) as handle:
        for line in handle:
            if not line.startswith('#'):
                break
            parts = line[1:].strip().split('\t')
            if len(parts) == 2:
                header[parts[0]] = parts[1]
    return header


def read_estimate_tsv(path):
    """({fragment: structure count}, {fragment: CG occurrences}) from this CLI's TSV."""
    structures, occurrences = {}, {}
    with open(path) as handle:
        for line in handle:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            try:
                fragment, structure_count, occurrence_count = line.split('\t')
                structures[fragment] = int(structure_count)
                occurrences[fragment] = int(occurrence_count)
            except ValueError as exc:
                raise ValueError(
                    f'{path}: malformed line {line!r}. Two-column files predate the '
                    'occurrence count and must be regenerated.') from exc
    return structures, occurrences


def main():
    import pickle as pkl
    from ligand_vdgs.generate_vdgs.extract_fragment_smiles import select_fragments

    args = parse_args()
    with open(args.frags_dict, 'rb') as handle:
        frags_dict = pkl.load(handle)
    # Every candidate, not a pre-thresholded subset -- selection is applied
    # against this estimate, so it must exist first.
    fragments = select_fragments(frags_dict, 0, args.max_size)
    structures, occurrences = estimate_fragment_counts(
        fragments, args.pdb_dir, sample_size=args.sample_size,
        num_procs=args.num_procs, seed=args.seed)

    parent = os.path.dirname(args.output)
    if parent:
        os.makedirs(parent, exist_ok=True)
    # tmp + rename: a killed run would otherwise leave a truncated TSV that
    # read_estimate_tsv parses cleanly but silently covers only part of the set.
    tmp_output = args.output + '.tmp'
    with open(tmp_output, 'w') as handle:
        handle.write(f'# fragment\tstructures\tCG occurrences, estimated from a '
                     f'{args.sample_size}-structure sample of {args.pdb_dir}\n')
        # Pins the TSV to the dict it was computed from -- a membership check
        # alone would miss a representative silently demoted to an alias.
        handle.write(f'# frags_dict\t{os.path.abspath(args.frags_dict)}\n')
        handle.write(f'# frags_dict_sha256\t{file_sha256(args.frags_dict)}\n')
        handle.write(f'# max_size\t{args.max_size}\n')
        handle.write(f'# pdb_dir\t{os.path.abspath(args.pdb_dir)}\n')
        # Identity (not path) survives copying the DB between machines, so
        # that's what the consumer compares (functions/db_identity).
        handle.write(f'# pdb_db_identity\t{identity_of(args.pdb_dir)["sha256"]}\n')
        handle.write(f'# sample_size\t{args.sample_size}\n')
        handle.write(f'# seed\t{args.seed}\n')
        scale = _LAST_SAMPLE_SCALE[0]
        # None only when fragments is empty (returns before sampling).
        handle.write(f'# sample_scale\t{"n/a" if scale is None else f"{scale:.6f}"}\n')
        for fragment in fragments:
            handle.write(f'{fragment}\t{structures[fragment]}\t'
                         f'{occurrences[fragment]}\n')
    os.replace(tmp_output, args.output)
    print(f'Wrote estimates for {len(fragments)} fragments to {args.output}.')


if __name__ == '__main__':
    main()
