"""Build the parent-database ligand roster used by ``fragment_database_ligs``.

Each roster type is ``(CCD code, observed canonical heavy-atom names)``.  A
biounit stores each type at most once, so NCS copies count as one observation.
OpenBabel-fallback and unreadable instances are excluded and counted in stats.
"""

import argparse
import concurrent.futures
import multiprocessing as mp
import os
import pickle as pkl
import time
from collections import defaultdict
from concurrent.futures.process import BrokenProcessPool

from ligand_vdgs.functions import ccd_templates, db_identity, parent_db
from ligand_vdgs.generate_vdgs.estimate_frag_cost import (
    MAX_LOST_SAMPLE_FRACTION, add_vdg_miner_paths)

def parse_args():
    parser = argparse.ArgumentParser(
        description='Record observed atom names for every ligand in a parent PDB database.')
    parser.add_argument('--pdb-dir', required=True,
                        help='Parent PDB database (prepwizard_BioLiP2_repaired).')
    parser.add_argument('--out', default='resources/ligand_roster.pkl')
    parser.add_argument('--logfile', default='logs/build_ligand_roster.log')
    parser.add_argument('--num-procs', type=int, default=1)
    parser.add_argument('--limit', type=int, default=None,
                        help='Testing only: stop after this many structures.')
    parser.add_argument('--force', action='store_true', help='Overwrite --out.')
    return parser.parse_args()

def pdb_paths(pdb_dir):
    return [path for _, path in parent_db.iter_structures(pdb_dir)]

def _init_worker():
    add_vdg_miner_paths()

def _roster_one(pdb_path):
    """Return ``(stem, instances, fallback_count, unreadable_count, bad_file)``."""
    from cg import read_ligand_blocks, suppress_stdout_stderr
    from ligand_vdgs.functions import parent_db
    from ligand_vdgs.functions.ligand_perception import (
        PERCEPTION_CCD_TEMPLATE, perceive_ligand_instance)

    stem = parent_db.stem_of(pdb_path)
    ligands = read_ligand_blocks(pdb_path)
    if ligands is None:
        return stem, [], 0, 0, True

    instances, fallback, unreadable = [], 0, 0
    with suppress_stdout_stderr():
        for key, block in ligands.items():
            perceived = perceive_ligand_instance(block, key[-1])
            if perceived is None:
                unreadable += 1
            elif perceived.provenance != PERCEPTION_CCD_TEMPLATE:
                fallback += 1
            else:
                names = _canonical_names(key[-1], block)
                if names is None:
                    unreadable += 1
                else:
                    instances.append((key[-1], names))
    return stem, instances, fallback, unreadable, False

def _canonical_names(resname, block):
    """Return canonical observed heavy-atom names, or ``None`` for an unknown name."""
    from ligand_vdgs.functions import ccd_templates

    try:
        template = ccd_templates.get_template(resname)
    except Exception:
        return None
    if template is None:
        return None

    names = set()
    for line in block.splitlines():
        if not line.startswith('HETATM') or line[76:78].strip() in ('H', 'D'):
            continue
        entry = ccd_templates.index_by_name(template).get(line[12:16].strip())
        if entry is None:
            return None
        if entry.element not in ('H', 'D'):
            names.add(entry.name)
    return tuple(sorted(names))

def build_roster(pdb_dir, num_procs=1, limit=None, progress_every=2000, log=print):
    """Build and return ``(types, biounits, stats, errors)``."""
    if num_procs < 1:
        raise ValueError(f'num_procs must be at least 1, got {num_procs}')
    paths = pdb_paths(pdb_dir)[:limit] if limit is not None else pdb_paths(pdb_dir)
    if not paths:
        raise SystemExit(f'[ERROR]: no PDB files found under {pdb_dir}.')

    types, type_index, biounits = [], {}, {}
    stats = defaultdict(int, num_structures=len(paths))
    errors, started = [], time.time()

    def tally(result):
        stem, instances, fallback, unreadable, bad_file = result
        stats['num_ob_fallback_instances'] += fallback
        stats['num_unreadable_instances'] += unreadable
        if bad_file:
            stats['num_unreadable_files'] += 1
            return
        stats['num_instances'] += len(instances)
        ids = set()
        for entry in instances:
            idx = type_index.get(entry)
            if idx is None:
                idx = type_index[entry] = len(types)
                types.append(entry)
            ids.add(idx)
        if ids:
            biounits[stem] = sorted(ids)
            stats['num_biounits_with_ligands'] += 1

    def progress(done):
        if progress_every and done % progress_every == 0:
            log(f'  {done}/{len(paths)} structures, {len(types)} distinct types, '
                f'{done / max(time.time() - started, 1e-9):.1f}/s')

    _init_worker()
    if num_procs == 1:
        for done, path in enumerate(paths, 1):
            try:
                tally(_roster_one(path))
            except Exception as exc:
                errors.append((path, repr(exc)))
            progress(done)
    else:
        pool = concurrent.futures.ProcessPoolExecutor(
            max_workers=num_procs, mp_context=mp.get_context('spawn'), initializer=_init_worker)
        try:
            futures = {pool.submit(_roster_one, path): path for path in paths}
            for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
                try:
                    tally(future.result())
                except BrokenProcessPool:
                    raise
                except Exception as exc:
                    errors.append((futures[future], repr(exc)))
                progress(done)
        except BrokenProcessPool as exc:
            raise RuntimeError(
                f'A roster worker died without returning a result ({exc}). Re-run with '
                f'fewer --num-procs or a larger -l mem_free.') from exc
        finally:
            pool.shutdown(wait=True, cancel_futures=True)

    stats.update(num_errored_structures=len(errors), num_types=len(types),
                 elapsed_s=round(time.time() - started, 1))
    return types, biounits, dict(stats), errors

def write_roster(out_path, payload):
    os.makedirs(os.path.dirname(out_path) or '.', exist_ok=True)
    tmp = out_path + '.tmp'
    with open(tmp, 'wb') as handle:
        pkl.dump(payload, handle, protocol=pkl.HIGHEST_PROTOCOL)
    os.replace(tmp, out_path)

def report(logfile, stats, errors, identity, ccd_identity, pdb_dir, limit):
    os.makedirs(os.path.dirname(logfile) or '.', exist_ok=True)
    with open(logfile, 'w') as handle:
        handle.write(f'{"=" * 25} Ligand roster {"=" * 25}\n')
        handle.write(f'pdb_dir:     {pdb_dir}\n' f'db_identity: {identity}\n'
                     f'ccd_identity: {ccd_identity}\n')
        if limit is not None:
            handle.write(f'*** PARTIAL: --limit {limit} was set. ***\n')
        for key in sorted(stats):
            handle.write(f'{key}: {stats[key]}\n')
        if errors:
            handle.write(f'\nStructures that raised ({len(errors)}):\n')
            for path, error in errors[:200]:
                handle.write(f'  {path}: {error}\n')
            if len(errors) > 200:
                handle.write(f'  ... and {len(errors) - 200} more\n')

def check_loss_threshold(stats, pdb_dir, logfile):
    """Refuse a roster if unreadable structures exceed the allowed fraction."""
    lost = stats.get('num_unreadable_files', 0) + stats['num_errored_structures']
    fraction = lost / stats['num_structures'] if stats['num_structures'] else 0.0
    if fraction > MAX_LOST_SAMPLE_FRACTION:
        raise SystemExit(
            f'[ERROR] {lost} of {stats["num_structures"]} structures under {pdb_dir} '
            f'were unreadable or raised ({100 * fraction:.0f}%, limit '
            f'{100 * MAX_LOST_SAMPLE_FRACTION:.0f}%). '
            f'See {logfile}.')

def main():
    args = parse_args()
    if os.path.exists(args.out) and not args.force:
        raise SystemExit(f'ERROR: {args.out} already exists. Pass --force to overwrite it.')
    if not os.path.isdir(args.pdb_dir):
        raise SystemExit(f'ERROR: --pdb-dir {args.pdb_dir} is not a directory.')

    identity = db_identity.identity_of(args.pdb_dir)
    from ligand_vdgs.functions import ligand_perception
    ligand_perception.require_template_store()
    ccd_identity = ccd_templates.store_identity()
    types, biounits, stats, errors = build_roster(
        args.pdb_dir, args.num_procs, args.limit)
    report(args.logfile, stats, errors, identity, ccd_identity, args.pdb_dir,
           args.limit)
    current_ccd_identity = ccd_templates.store_identity()
    if current_ccd_identity != ccd_identity:
        raise SystemExit(
            f'[ERROR]: CCD template store changed while building the roster; started with '
            f'{ccd_identity}, now {current_ccd_identity}.')
    if not biounits:
        raise SystemExit(
            f'[ERROR]: {stats["num_structures"]} structures yielded no templated ligand '
            f'instance; refusing to write an empty roster. See {args.logfile}.')
    check_loss_threshold(stats, args.pdb_dir, args.logfile)
    write_roster(args.out, {
        'db_identity': identity, 'pdb_dir': os.path.abspath(args.pdb_dir),
        'ccd_identity': ccd_identity,
        'partial': args.limit is not None, 'types': types, 'biounits': biounits,
        'stats': stats})
    print(f'Wrote {args.out}: {len(biounits)} biounits, {len(types)} distinct types. '
          f'Log: {args.logfile}', flush=True)

if __name__ == '__main__':
    main()
