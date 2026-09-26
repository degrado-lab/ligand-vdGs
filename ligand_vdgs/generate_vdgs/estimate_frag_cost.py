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
import sys
import random
import argparse
import multiprocessing as mp
import concurrent.futures
# Imported, not attribute-accessed: on 3.10 concurrent.futures exposes only public
# class names via __getattr__, so `concurrent.futures.process.X` raises
# AttributeError exactly when a worker has died.
from concurrent.futures.process import BrokenProcessPool
from ligand_vdgs.functions.db_identity import identity_of
from ligand_vdgs.functions.utils import file_sha256
from ligand_vdgs.functions import parent_db
from ligand_vdgs.functions.fragment_prefilter import build_prefilter, prefilter_candidates

HERE = os.path.dirname(os.path.abspath(__file__))


def add_vdg_miner_paths():
    # The vdG-miner submodule is not a package, so it still needs sys.path.
    path = os.path.join(HERE, "..", "..", "external", "vdG-miner", "vdg_miner", "vdg")
    if path not in sys.path:
        sys.path.append(path)


# This estimates mining cost, not fragment membership. Structure counts size jobs;
# occurrence counts determine whether a fragment is worth mining.
DEFAULT_SAMPLE_SIZE = 3000

# Above this share of the sample lost to unreadable or erroring structures, the
# extrapolation is reported as a failure rather than scaled up from what is left.
MAX_LOST_SAMPLE_FRACTION = 0.1

# Set by estimate_fragment_counts; main() records it in the TSV header (counts are
# already scaled).
_LAST_SAMPLE_SCALE = [None]


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
    _WORKER_PREFILTER = build_prefilter(fragments, _WORKER_RING_CONSTRAINTS)


def _count_one(pdb_path):
    """Return occurrence counts and read/perception failure statistics."""
    from cg import (match_satisfies_ring_sizes, read_ligand_blocks,
                    suppress_stdout_stderr)
    from ligand_vdgs.functions.ligand_perception import (
        perceive_ligand_instance, phantom_atom_indices)
    ligands = read_ligand_blocks(pdb_path)
    if ligands is None:
        return {}, 0, True, False
    if not ligands:
        return {}, 0, False, False
    counts = {}
    read_failures = 0
    with suppress_stdout_stderr():
        # Uses the miner's perception (deletes hydrogens; without it every
        # `D<n>`-annotated key counts zero). None means unreadable; the miner skips
        # it too, so it contributes zero but is tallied for the warning.
        for key, block in ligands.items():
            perceived = perceive_ligand_instance(block, key[-1])
            if perceived is None:
                read_failures += 1
                continue
            mol = perceived.obmol
            phantoms = phantom_atom_indices(mol)
            for i in prefilter_candidates(mol, _WORKER_PREFILTER):
                pattern = _WORKER_FRAGMENTS[i]
                if pattern.Match(mol):
                    constraints = _WORKER_RING_CONSTRAINTS[i]
                    n = sum(1 for m in pattern.GetUMapList()
                            if match_satisfies_ring_sizes(mol, m, constraints)
                            and phantoms.isdisjoint(m))
                    if n:
                        counts[int(i)] = counts.get(int(i), 0) + n
    return counts, read_failures, False, read_failures == len(ligands)


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
    # Also validates the SMARTS/ring constraints early, before spinning up a pool.
    _init_worker(fragments)

    stats = {'failures': 0, 'unreadable': 0, 'all_failed': 0}
    errored = []

    def _tally(result):
        matched, n_failed, was_unreadable, was_all_failed = result
        stats['failures'] += n_failed
        stats['unreadable'] += bool(was_unreadable)
        stats['all_failed'] += bool(was_all_failed)
        for i, n in matched.items():
            struct_hits[i] += 1
            occurrences[i] += n

    if num_procs > 1:
        ctx = mp.get_context("spawn")
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
            pool.shutdown(wait=True, cancel_futures=True)
    else:
        for pdb_path in sample:
            _consume_result(pdb_path, lambda p=pdb_path: _count_one(p),
                            _tally, errored)
    if stats['failures']:
        print(f"[WARNING] {stats['failures']} ligand block(s) in the {len(sample)}-structure "
              f'sample failed perception; the miner skips them too, so they count as '
              f'zero.', file=sys.stderr)
    if errored:
        print(f'[WARNING] {len(errored)} of {len(sample)} sampled structures raised '
              f'while being counted and were excluded; e.g. '
              f'{errored[0][0]}: {errored[0][1]}', file=sys.stderr)
    # Perception failure is deterministic and the miner skips the same ligands, so a
    # structure whose every ligand failed is a real zero and stays in the denominator.
    read_ok = len(sample) - stats['unreadable'] - len(errored)
    if stats['unreadable']:
        print(f"[WARNING] {stats['unreadable']} of {len(sample)} sampled structures were "
              f'unreadable; extrapolating from the {read_ok} that were read.',
              file=sys.stderr)
    if stats['all_failed']:
        print(f"[WARNING] {stats['all_failed']} of {len(sample)} sampled structures had "
              f'every ligand fail perception; counted as zero, since the miner '
              f'skips them too.', file=sys.stderr)
    if read_ok == 0:
        raise RuntimeError(
            f'none of the {len(sample)} sampled structures under {pdb_dir} could be '
            'read; refusing to emit an all-zero estimate.')
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

def sampling_upper_bound(est_count, scale, z=2.0):
    """Upper confidence bound on an extrapolated count, for sizing a job.

    TSV counts are point estimates from one draw; loss from under-sizing is
    asymmetric with over-sizing (a too-small tier gets killed at h_rt; a too-large
    one only costs queue latency), so size on the upper end: `est_count / scale`
    recovers the raw sampled count (relative SE ~1/sqrt(raw), Poisson), inflated by
    `z` of those.

    `scale == 1` (a --census run with no read failures) is a no-op case: `est_count`
    is then an exact count, not a sample draw, so there is no sampling error to
    widen against.
    """
    if not est_count or not scale or scale <= 0 or scale == 1.0:
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
    parser.add_argument('--census', action='store_true',
                        help="Use every structure under --pdb-dir instead of "
                             "--sample-size, so counting a growing mirror never "
                             "silently degrades into a sample.")
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


def main():
    from ligand_vdgs.functions import ligand_structure
    from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (load_frags_dict,
                                                                   select_fragments)

    args = parse_args()
    frags_dict, _support, _meta = load_frags_dict(args.frags_dict)
    fragments = select_fragments(frags_dict, 0, args.max_size)
    frags_dict_sha256 = file_sha256(args.frags_dict)
    pdb_db_identity = identity_of(args.pdb_dir)['sha256']
    structures, occurrences = estimate_fragment_counts(
        fragments, args.pdb_dir, sample_size=None if args.census else args.sample_size,
        num_procs=args.num_procs, seed=args.seed)

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    tmp_output = args.output + '.tmp'
    with open(tmp_output, 'w') as handle:
        scale = _LAST_SAMPLE_SCALE[0]
        headers = [
            f'# fragment\tstructures\tCG occurrences, estimated from '
            f'{"a census of every structure" if args.census else f"a {args.sample_size}-structure sample"} '
            f'in {args.pdb_dir}',
            f'# frags_dict\t{os.path.abspath(args.frags_dict)}',
            f'# frags_dict_sha256\t{frags_dict_sha256}',
            f'# key_schema\t{ligand_structure.KEY_SCHEMA}',
            f'# max_size\t{args.max_size}',
            f'# pdb_dir\t{os.path.abspath(args.pdb_dir)}',
            f'# pdb_db_identity\t{pdb_db_identity}',
            f'# sample_size\t{"census" if args.census else args.sample_size}',
            f'# seed\t{args.seed}',
            f'# sample_scale\t{"n/a" if scale is None else f"{scale:.6f}"}',]
        handle.write('\n'.join(headers) + '\n')
        for fragment in fragments:
            handle.write(f'{fragment}\t{structures[fragment]}\t'
                         f'{occurrences[fragment]}\n')
    os.replace(tmp_output, args.output)
    print(f'Wrote estimates for {len(fragments)} fragments to {args.output}.')


if __name__ == '__main__':
    main()
