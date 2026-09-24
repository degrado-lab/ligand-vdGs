# Generate per-fragment SGE job scripts for running the vdG generation pipeline.

import os
import re
import json
import fcntl
import argparse
from ligand_vdgs.functions import utils
from ligand_vdgs.functions.utils import _int_or_none, file_sha256
from ligand_vdgs.functions.db_identity import identity_of
from ligand_vdgs.functions.Frags import check_vdg_job_status
from ligand_vdgs.functions.vdg_npz_utils import load_fragment_aliases
from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
    load_frags_dict, alias_kind, fragment_dict_keys, prepare_fragments_cached,
    select_fragments, write_fragment_aliases)
from ligand_vdgs.generate_vdgs.estimate_frag_cost import (
    DEFAULT_SAMPLE_SIZE, estimate_fragment_counts, sampling_upper_bound, _LAST_SAMPLE_SCALE)

TIER_1 = (20000, 10, '48:00:00')
RESOURCE_TIERS = (TIER_1, (float('inf'), 20, '48:00:00'))
MEM_FREE_PER_SLOT = {10: '1G', 20: '2G'}  # mem_free is PER SLOT under -pe smp; fixed, not user-configurable.
DEFAULT_MEM_FREE = '2G'
SCRATCH = '5G'
WRAPPER_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vdg_generation_wrapper.py')
PROVENANCE_FILENAME = 'library_provenance.json'

def _h_rt_to_hours(h_rt):
    hours, minutes, seconds = (int(part) for part in h_rt.split(':'))
    return hours + minutes / 60 + seconds / 3600

def _h_rt(value):
    match = re.fullmatch(r'(\d+):(\d{1,2}):(\d{1,2})', value)
    if not match or int(match.group(2)) > 59 or int(match.group(3)) > 59:
        raise argparse.ArgumentTypeError(f"expected HH:MM:SS (e.g. 36:00:00), got {value!r}")
    return value

def resources_for(est_occurrences, fixed_num_procs=None, fixed_h_rt=None):
    for upper, slots, h_rt in RESOURCE_TIERS:
        if est_occurrences < upper:
            break
    slots = int(fixed_num_procs) if fixed_num_procs is not None else slots
    h_rt = fixed_h_rt if fixed_h_rt is not None else h_rt
    return slots, h_rt

def provenance_path(vdg_lib_dir):
    return os.path.join(vdg_lib_dir, PROVENANCE_FILENAME)

def write_provenance(vdg_lib_dir, args, key_schema=None):
    record = {'frags_dict': os.path.abspath(args.frags_dict), 'frags_dict_sha256': file_sha256(args.frags_dict),
              'max_size': args.max_size, 'min_support': args.min_support,
              'parent_pdb_dir': os.path.abspath(args.pdb_dir), 'parent_db_identity': identity_of(args.pdb_dir)['sha256']}
    if key_schema is not None:
        record['key_schema'] = key_schema
    os.makedirs(vdg_lib_dir, exist_ok=True)
    final, tmp = provenance_path(vdg_lib_dir), provenance_path(vdg_lib_dir) + '.tmp'
    with open(tmp, 'w') as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write('\n')
    os.replace(tmp, final)

def check_library_vocabulary(vdg_lib_dir, key_schema):
    path = provenance_path(vdg_lib_dir)
    if not os.path.isfile(path):
        return
    with open(path) as handle:
        recorded = json.load(handle).get('key_schema')
    if recorded != key_schema:
        raise SystemExit(f'[ERROR] {vdg_lib_dir}: fragment-vocabulary mismatch (recorded {recorded!r}, '
                          f'this run emits {key_schema!r}). Use a new --vdg-lib-dir, or rebuild.')

def check_recorded_inputs(recorded, args, source, top_up, max_size_is_error):
    problems = []
    recorded_identity = recorded.get('db_identity')
    if recorded_identity is None:
        problems.append('missing parent-database identity')
    elif recorded_identity != identity_of(args.pdb_dir)['sha256']:
        problems.append(f'db identity: recorded {recorded_identity[:12]}..., now '
                         f'{identity_of(args.pdb_dir)["sha256"][:12]}...')

    recorded_pdb_dir = recorded.get('pdb_dir')
    if recorded_pdb_dir is not None and os.path.abspath(recorded_pdb_dir) != os.path.abspath(args.pdb_dir):
        print(f'[WARNING] {source}: --pdb-dir was {recorded_pdb_dir}, now {os.path.abspath(args.pdb_dir)} '
              f'(informational; contents are what is checked).')

    recorded_sha, current_sha = recorded.get('frags_dict_sha256'), file_sha256(args.frags_dict)
    if recorded_sha is None:
        problems.append('missing frags_dict_sha256')
    elif recorded_sha != current_sha:
        detail = f'sha256: recorded {recorded_sha[:12]}..., now {current_sha[:12]}...'
        if top_up:
            print(f'[WARNING] {source}: {detail}. Top-up grows the dict by design -- regenerate to be sure.')
        else:
            problems.append(detail)

    recorded_max_size = recorded.get('max_size')
    if recorded_max_size is None:
        problems.append('missing max_size')
    elif int(recorded_max_size) != args.max_size:
        detail = f'--max-size: recorded {recorded_max_size}, now {args.max_size}'
        (problems.append if max_size_is_error else lambda d: print(f'[WARNING] {source}: {d}.'))(detail)

    if problems:
        raise SystemExit(f'[ERROR] {source} mismatches this run ({"; ".join(problems)}). '
                          'Regenerate the record, or restore the recorded inputs.')

def check_frags_dict_identity(frags_meta, args):
    recorded = frags_meta.get('db_identity')
    if not recorded:
        raise SystemExit(f'[ERROR] {args.frags_dict} records no parent-database identity; regenerate it.')
    if recorded.get('sha256') != identity_of(args.pdb_dir)['sha256']:
        raise SystemExit(f'[ERROR] {args.frags_dict} was built from a different --pdb-dir than this run.')

def check_provenance(vdg_lib_dir, args):
    path = provenance_path(vdg_lib_dir)
    if not os.path.isfile(path):
        raise SystemExit(f'[ERROR] {path} is missing; rebuild the library before adding fragments.')
    with open(path) as handle:
        record = json.load(handle)
    check_recorded_inputs(
        {'pdb_dir': record.get('parent_pdb_dir'), 'db_identity': record.get('parent_db_identity'),
         'frags_dict_sha256': record.get('frags_dict_sha256'), 'max_size': record.get('max_size')},
        args, f"{path} (this library's provenance)", top_up=args.include_only, max_size_is_error=True)

def parse_args():
    parser = argparse.ArgumentParser(description="Create SGE submission scripts for vdG generation.")
    parser.add_argument('--frags-dict', default='resources/database_frags_dict.pkl', help="Path to database_frags_dict.pkl.")
    parser.add_argument('--template', default='resources/frag_sge_template.sh', help="Path to SGE job script template.")
    parser.add_argument('--include', nargs='*', default=[], metavar='SMARTS', help="Fragment SMARTS to build regardless of estimated count.")
    parser.add_argument('--resume', action='store_true',
                         help="Emit scripts only for unfinished fragments (per Frags.check_vdg_job_status), and allow "
                              "--sge-out-dir to be non-empty. Use --include-only to ADD fragments instead.")
    parser.add_argument('--include-only', action='store_true',
                         help="Top-up mode: build scripts for --include fragments only, skipping the count threshold.")
    parser.add_argument('--min-support', default=None, type=int,
                         help="Min distinct parent biounits containing a fragment (DR-5), from the fragment dict. Not "
                              "the same unit as the old --min-instances; required for a full/resumed build, unused "
                              "under --include-only.")
    parser.add_argument('--max-size', default=5, type=int, help="Max fragment heavy-atom count. Default: 5.")
    parser.add_argument('--sge-out-dir', default='ligand_vdgs/generate_vdgs/frag_submit_scripts/', help="Output directory for SGE scripts.")
    parser.add_argument('--vdg-lib-dir', default='/wynton/home/degradolab/skt/docking/frag_lib', help="Output directory for vdG library.")
    parser.add_argument('--log-dir', default='/wynton/home/degradolab/skt/docking/frag_sge_logs', help="SGE log directory.")
    parser.add_argument('--pdb-dir', default='/wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2_repaired/', help="Path to parent PDB database.")
    parser.add_argument('--max-num-clus', default=None, type=_int_or_none, help="Max vdGs to cluster per subset. Default: no limit.")
    parser.add_argument('--h-rt', default=None, type=_h_rt, help="Fixed h_rt for every fragment. Default: tiered alongside slot count.")
    parser.add_argument('--num-procs', default=None, type=_int_or_none, help="Fixed slot count for every fragment. Default: tiered by estimated cost.")
    parser.add_argument('--sample-size', default=DEFAULT_SAMPLE_SIZE, type=int, help=f"Structures sampled for the inline cost estimate. Default: {DEFAULT_SAMPLE_SIZE}.")
    parser.add_argument('--estimate-procs', default=10, type=int, help="Worker processes for the inline cost estimate. Default: 10.")
    parser.add_argument('--subset-sizes', nargs='+', type=int, default=[1, 2], choices=[1, 2], help="vdG subset sizes to build. Default: 1 2.")
    args = parser.parse_args()
    if args.include_only and not args.include:
        parser.error("--include-only needs at least one --include fragment; with none it would generate an empty fleet.")
    return args

def main():
    args = parse_args()

    for path, flag, check, exc in [
            (args.frags_dict, '--frags-dict', os.path.isfile, FileNotFoundError),
            (args.template, '--template', os.path.isfile, FileNotFoundError),
            (args.log_dir, '--log-dir', os.path.isdir, NotADirectoryError),
            (args.pdb_dir, '--pdb-dir', os.path.isdir, NotADirectoryError)]:
        if not check(path):
            raise exc(f"{flag} path does not exist: {path}")
    if not os.path.isfile(WRAPPER_PATH):
        raise SystemExit(f'[ERROR] wrapper not found at {WRAPPER_PATH}.')
    if args.min_support is None and not args.include_only:
        raise SystemExit('[ERROR] --min-support is required for a full or resumed build (min distinct parent '
                          'biounits containing a fragment, DR-5). No default.')

    os.makedirs(args.sge_out_dir, exist_ok=True)
    if not args.resume and [f for f in os.listdir(args.sge_out_dir) if not f.startswith('.')]:
        raise FileExistsError(f"{args.sge_out_dir} already has files; refusing to overwrite. "
                               "Use --resume to write only unfinished scripts.")
    if args.include_only or args.resume:
        check_provenance(args.vdg_lib_dir, args)

    frags_dict, support_pooled, frags_meta = load_frags_dict(args.frags_dict)
    check_frags_dict_identity(frags_meta, args)
    check_library_vocabulary(args.vdg_lib_dir, frags_meta.get('key_schema'))
    prepared = prepare_fragments_cached(frags_dict, args.max_size, args.frags_dict)
    candidates = select_fragments(frags_dict, 0, args.max_size, prepared=prepared)

    print(f'Counting {len(candidates)} candidate fragments over {args.sample_size} sampled structures...')
    _, occurrences = estimate_fragment_counts(candidates, args.pdb_dir, sample_size=args.sample_size,
                                                num_procs=args.estimate_procs)
    sample_scale = _LAST_SAMPLE_SCALE[0]

    smiles_to_run, aliases = select_fragments(frags_dict, args.min_support, args.max_size, support=support_pooled,
                                               return_aliases=True, prepared=prepared, include=args.include,
                                               include_only=args.include_only)

    skipped_finished, partial_dirs = [], []
    if args.resume:
        smiles_to_run, skipped_finished, partial_dirs = partition_by_completion(smiles_to_run, args.vdg_lib_dir)
        if not smiles_to_run:
            print(f'Resume: all {len(skipped_finished)} selected fragment(s) already have '
                  f"'Job completed.' in their log. Nothing to submit.")
            return

    if sample_scale is None:
        print('[WARNING] cost estimate has no sample_scale; tiers use the point estimate and carry '
              '~15% seed-to-seed tier instability. Regenerate to fix.')

    tier_counts, written_scripts, submission_rows = {}, [], []
    with open(args.template) as f:
        template_lines = f.readlines()

    for smiles in smiles_to_run:
        tier_count = sampling_upper_bound(occurrences[smiles], sample_scale)
        slots, h_rt = resources_for(tier_count, fixed_num_procs=args.num_procs, fixed_h_rt=args.h_rt)
        per_frag = {'$WRAPPER': WRAPPER_PATH, '$LOG_DIR': args.log_dir, '$PDB_DIR': args.pdb_dir,
                    '$OUTPUT_DIR': args.vdg_lib_dir, '$MAX_NUM_CLUS': str(args.max_num_clus),
                    '$SUBSET_SIZES': ' '.join(str(s) for s in sorted(set(args.subset_sizes))),
                    '$SCRATCH': SCRATCH, '$NUM_PROCS': str(slots), '$RUN_TIME': h_rt,
                    '$MEM_FREE': MEM_FREE_PER_SLOT.get(slots, DEFAULT_MEM_FREE)}
        tier_counts[(slots, h_rt)] = tier_counts.get((slots, h_rt), 0) + 1
        output_script(template_lines, smiles, args.sge_out_dir, per_frag)
        script_path = os.path.join(args.sge_out_dir, utils.smiles_to_filename(smiles) + '.sh')
        written_scripts.append(script_path)
        submission_rows.append((script_path, smiles, slots, h_rt, tier_count))

    alias_path = os.path.join(args.vdg_lib_dir, 'fragment_aliases.tsv')
    dict_keys = fragment_dict_keys(frags_dict)
    if args.include_only or args.resume:
        os.makedirs(args.vdg_lib_dir, exist_ok=True)
        with open(alias_path + '.lock', 'w') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            merged = dict(load_fragment_aliases(args.vdg_lib_dir))
            merged.update(aliases)
            write_fragment_aliases(alias_path, merged, dict_keys)
        aliases = merged
    else:
        write_fragment_aliases(alias_path, aliases, dict_keys)
    if not args.resume and not args.include_only:
        write_provenance(args.vdg_lib_dir, args, frags_meta.get('key_schema'))

    manifest = os.path.join(args.sge_out_dir, 'submission_order.tsv')
    submission_rows.sort(key=lambda row: (-row[2], -row[4], row[0]))
    with open(manifest, 'w') as handle:
        handle.write('# Submit THESE scripts, in this order. Do not use `qsub *.sh`: a resumed or topped-up '
                      'directory still holds scripts for fragments that already finished, and the glob silently '
                      'redoes them.\n')
        handle.write('order\tscript\tfragment\tslots\th_rt\ttier_count\n')
        for order, (path, smiles, slots, h_rt, tier_count) in enumerate(submission_rows):
            handle.write(f'{order}\t{path}\t{smiles}\t{slots}\t{h_rt}\t{tier_count}\n')
    print(f'Submission order written to {manifest}. Submit exactly these, in order:\n'
          f'  python ligand_vdgs/generate_vdgs/submit_frag_jobs_by_size.py {args.sge_out_dir}')

    if args.resume:
        print(f'Resume: {len(skipped_finished)} fragment(s) already finished and were skipped; wrote '
              f'{len(written_scripts)} script(s) to {args.sge_out_dir}.')
        if partial_dirs:
            print(f'[WARNING] {len(partial_dirs)} fragment(s) have a leftover output dir from a killed job; '
                  f'the wrapper refuses a non-empty output dir, so those jobs will die immediately. Remove them:')
            for smiles in partial_dirs:
                print(f'  {os.path.join(args.vdg_lib_dir, utils.smiles_to_filename(smiles))}')
    elif args.include_only:
        print(f'Top-up: created scripts for {len(smiles_to_run)} requested fragment(s) in {args.sge_out_dir} '
              f'(count threshold not applied).')
    else:
        print(f'Created scripts for {len(smiles_to_run)} of {len(candidates)} candidate fragments '
              f'(>= {args.min_support} parent biounits'
              f'{f", plus {len(args.include)} requested" if args.include else ""}) in {args.sge_out_dir}.')

    print('Resources requested: ' + ', '.join(f'{n} fragment(s) at -pe smp {slots}, h_rt {h_rt}'
                                                for (slots, h_rt), n in sorted(tier_counts.items())))
    slot_hours = {key: n * key[0] * _h_rt_to_hours(key[1]) for key, n in tier_counts.items()}
    print('Worst-case slot-hours (ceiling, not a forecast): ' +
          ', '.join(f'{slot_hours[key]:,.0f} at -pe smp {key[0]}/{key[1]}' for key in sorted(tier_counts)) +
          f'; TOTAL {sum(slot_hours.values()):,.0f} slot-hours over {sum(tier_counts.values())} job(s), '
          f'peak {sum(k[0] * n for k, n in tier_counts.items()):,} slots if every job ran at once.')
    if aliases:
        print(f'Collapsed {len(aliases)} protonation variant(s); wrote {alias_path}.')
        promoted = sorted(r for r in set(aliases.values()) if alias_kind(r, dict_keys) == 'promoted')
        if promoted:
            print(f'{len(promoted)} representative(s) are promoted charge-stripped keys absent from the fragment '
                  f'dict (library dirs named for a SMARTS no CCD ligand is drawn with). See kind=promoted in '
                  f'{alias_path} and scripts/lookup_fragment_key.py: {promoted}')

def partition_by_completion(smiles_to_run, vdg_lib_dir):
    unfinished, finished, partial = [], [], []
    for smiles in smiles_to_run:
        label = utils.smiles_to_filename(smiles)
        if check_vdg_job_status(label, vdg_lib_dir):
            finished.append(smiles)
            continue
        unfinished.append(smiles)
        if os.path.isdir(os.path.join(vdg_lib_dir, label)):
            partial.append(smiles)
    return unfinished, finished, partial

def output_script(template_lines, smiles, sge_out_dir, replace):
    local_replace = dict(replace, **{'$SMILES': f'"{smiles}"', '$CG': f'"{utils.smiles_to_filename(smiles)}"',
                                      '$JOB_NAME': utils.smiles_to_job_name(smiles)})
    with open(os.path.join(sge_out_dir, utils.smiles_to_filename(smiles) + '.sh'), 'w') as f:
        start_copy = False
        for line in template_lines:
            start_copy = start_copy or line.startswith('#!/bin/bash')
            if start_copy:
                f.write(re.sub('|'.join(re.escape(key) for key in local_replace),
                                lambda m: local_replace[m.group(0)], line))

if __name__ == '__main__':
    main()
