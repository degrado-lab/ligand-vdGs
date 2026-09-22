# Generate per-fragment SGE job scripts for running the vdG generation pipeline.

import os
import re
import json
import shlex
import fcntl
import getpass
import argparse
import pickle as pkl
import subprocess
import xml.etree.ElementTree as ET
from ligand_vdgs.functions import utils
from ligand_vdgs.functions.utils import _int_or_none, file_sha256
from ligand_vdgs.functions.db_identity import identity_of
from ligand_vdgs.functions.Frags import check_vdg_job_status
from ligand_vdgs.functions.vdg_npz_utils import load_fragment_aliases
from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
    load_frags_dict,
    alias_kind, fragment_dict_keys, prepare_fragments_cached,
    select_fragments, write_fragment_aliases)
from ligand_vdgs.generate_vdgs.estimate_frag_cost import (
    DEFAULT_SAMPLE_SIZE, estimate_fragment_counts, read_estimate_tsv,
    read_estimate_header, sampling_upper_bound, _LAST_SAMPLE_SCALE)

TIER_1 = (20000, 10, '48:00:00')
RESOURCE_TIERS = (
    TIER_1,
    (float('inf'), 20, '48:00:00'),
)
# mem_free is PER SLOT under -pe smp. 
MEM_FREE_PER_SLOT = {10: '1G', 20: '2G'}
TIER_KEY_COLUMN = 'occurrences'
TOP_TIER_LOWER_BOUND = RESOURCE_TIERS[-2][0]
TOP_TIER_H_RT = RESOURCE_TIERS[-1][2]
SCRATCH = '5G'

def _h_rt_to_hours(h_rt):
    hours, minutes, seconds = (int(part) for part in h_rt.split(':'))
    return hours + minutes / 60 + seconds / 3600

def _h_rt(value):
    match = re.fullmatch(r'(\d+):(\d{1,2}):(\d{1,2})', value)
    if not match or int(match.group(2)) > 59 or int(match.group(3)) > 59:
        raise argparse.ArgumentTypeError(
            f"expected HH:MM:SS (e.g. 36:00:00), got {value!r}")
    return value

def resources_for(est_occurrences, max_h_rt, fixed_num_procs=None, fixed_h_rt=None):
    tiers = RESOURCE_TIERS
    for upper, slots, h_rt in tiers:
        if est_occurrences < upper:
            break
    if fixed_num_procs is not None:
        slots = int(fixed_num_procs)
    if fixed_h_rt is not None:
        h_rt = fixed_h_rt
    if _h_rt_to_hours(h_rt) > _h_rt_to_hours(max_h_rt):
        h_rt = max_h_rt
    return slots, h_rt

WRAPPER_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            'vdg_generation_wrapper.py')

PROVENANCE_FILENAME = 'library_provenance.json'

def provenance_path(vdg_lib_dir):
    return os.path.join(vdg_lib_dir, PROVENANCE_FILENAME)

def _selection_inputs(args):
    return {'frags_dict': os.path.abspath(args.frags_dict),
            'frags_dict_sha256': file_sha256(args.frags_dict),
            'max_size': args.max_size}

def write_provenance(vdg_lib_dir, args, key_schema=None):
    record = _selection_inputs(args)
    record['min_support'] = args.min_support
    if key_schema is not None:
        record['key_schema'] = key_schema
    record['parent_pdb_dir'] = os.path.abspath(args.pdb_dir)
    record['parent_db_identity'] = identity_of(args.pdb_dir)['sha256']
    os.makedirs(vdg_lib_dir, exist_ok=True)
    final = provenance_path(vdg_lib_dir)
    tmp = final + '.tmp'
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
    if recorded == key_schema:
        return
    raise SystemExit(
        f'[ERROR] {vdg_lib_dir} is an existing library with a different fragment '
        f'vocabulary: provenance records {recorded!r}, this run emits {key_schema!r} '
        f'(None means no key_schema at all, i.e. pre-annotated-key). Point '
        f'--vdg-lib-dir at a new directory, or rebuild from scratch.')

def check_recorded_inputs(recorded, args, source, top_up, max_size_is_error):
    problems = []

    recorded_identity = recorded.get('db_identity')
    if recorded_identity is None:
        problems.append('missing parent-database identity')
    else:
        current = identity_of(args.pdb_dir)
        if recorded_identity != current['sha256']:
            problems.append(
                f'db identity: recorded {recorded_identity[:12]}..., '
                f'{os.path.abspath(args.pdb_dir)} is {current["sha256"][:12]}...')

    recorded_pdb_dir = recorded.get('pdb_dir')
    if (recorded_pdb_dir is not None
            and os.path.abspath(recorded_pdb_dir) != os.path.abspath(args.pdb_dir)):
        print(f'[WARNING] {source} used --pdb-dir {recorded_pdb_dir}, this run uses '
              f'{os.path.abspath(args.pdb_dir)} (informational; contents are what '
              f'is checked).')

    recorded_sha = recorded.get('frags_dict_sha256')
    current_sha = file_sha256(args.frags_dict)
    dict_name = os.path.basename(args.frags_dict)
    if recorded_sha is None:
        problems.append('missing frags_dict_sha256')
    elif recorded_sha != current_sha:
        detail = f'{dict_name} sha256: recorded {recorded_sha[:12]}..., now {current_sha[:12]}...'
        if top_up:
            print(f'[WARNING] {source}: {detail}. Top-up grows the dict by design, but '
                  f'grouping may have changed -- regenerate to be sure.')
        else:
            problems.append(detail)

    recorded_max_size = recorded.get('max_size')
    if recorded_max_size is None:
        problems.append('missing max_size')
    elif int(recorded_max_size) != args.max_size:
        detail = f'--max-size: recorded {recorded_max_size}, now {args.max_size}'
        if max_size_is_error:
            problems.append(detail)
        else:
            print(f'[WARNING] {source}: {detail}.')

    if problems:
        raise SystemExit(
            f'[ERROR] {source} does not match this run ({"; ".join(problems)}). '
            'Recorded inputs decide which fragments are built and which key a '
            'SMARTS resolves to. Regenerate the record, or restore the recorded inputs.')

def check_frags_dict_identity(frags_meta, args):
    recorded = frags_meta.get('db_identity')
    if not recorded:
        raise SystemExit(f'[ERROR] {args.frags_dict} records no parent-database identity; '
                         'regenerate it with the current fragment-dictionary writer.')
    current = identity_of(args.pdb_dir)
    if recorded.get('sha256') != current['sha256']:
        raise SystemExit(
            f'[ERROR] {args.frags_dict} was built from a different parent database than '
            '--pdb-dir describes; support thresholds must come from the same database.')

def check_provenance(vdg_lib_dir, args):
    path = provenance_path(vdg_lib_dir)
    if not os.path.isfile(path):
        raise SystemExit(f'[ERROR] {path} is missing; rebuild the library with the '
                         'current writer before adding fragments.')
    with open(path) as handle:
        record = json.load(handle)
    recorded = {'pdb_dir': record.get('parent_pdb_dir'),
                'db_identity': record.get('parent_db_identity'),
                'frags_dict_sha256': record.get('frags_dict_sha256'),
                'max_size': record.get('max_size')}
    check_recorded_inputs(recorded, args, f'{path} (this library\'s provenance)',
                          top_up=args.include_only, max_size_is_error=True)

def check_estimate_header(header, args):
    recorded = dict(header)
    recorded['db_identity'] = header.get('pdb_db_identity')
    check_recorded_inputs(recorded, args, args.frag_cost_estimate,
                          top_up=args.include_only, max_size_is_error=False)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Create SGE submission scripts for vdG generation.")
    parser.add_argument('--frags-dict', default='resources/database_frags_dict.pkl',
                        help="Path to database_frags_dict.pkl.")
    parser.add_argument('--template', default='resources/frag_sge_template.sh',
                        help="Path to SGE job script template.")
    parser.add_argument('--include', nargs='*', default=[], metavar='SMARTS',
                        help="Fragment SMARTS to build regardless of estimated count.")
    parser.add_argument('--include-file', default=None,
                        help="File of --include SMARTS, one per line ('#'/blank ignored).")
    parser.add_argument('--resume', action='store_true',
                        help="Emit scripts only for unfinished fragments (per "
                             "Frags.check_vdg_job_status), and allow --sge-out-dir to be "
                             "non-empty. Use --include-only to ADD fragments instead.")
    parser.add_argument('--clear-partial', action='store_true',
                        help="With --resume, delete each fragment's leftover output dir "
                             "before rerunning. Off by default: it deletes inside the "
                             "library.")
    parser.add_argument('--include-only', action='store_true',
                        help="Top-up mode: build scripts for --include fragments only, "
                             "skipping the count threshold.")
    parser.add_argument('--min-support', default=None, type=int,
                        help="Min distinct parent biounits containing a fragment (DR-5), "
                             "from the fragment dict. Not the same unit as the old "
                             "--min-instances; required for a full/resumed build, unused "
                             "under --include-only.")
    parser.add_argument('--max-size', default=5, type=int,
                        help="Max fragment heavy-atom count. Default: 5.")
    parser.add_argument('--sge-out-dir',
                        default='ligand_vdgs/generate_vdgs/frag_submit_scripts/',
                        help="Output directory for SGE scripts.")
    parser.add_argument('--vdg-lib-dir',
                        default='/wynton/home/degradolab/skt/docking/frag_lib',
                        help="Output directory for vdG library.")
    parser.add_argument('--log-dir',
                        default='/wynton/home/degradolab/skt/docking/frag_sge_logs',
                        help="SGE log directory.")
    parser.add_argument('--pdb-dir',
                        default='/wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2_repaired/',
                        help="Path to parent PDB database.")
    parser.add_argument('--max-num-clus', default=None, type=_int_or_none,
                        help="Max vdGs to cluster per subset. Default: no limit.")
    parser.add_argument('--h-rt', default=None, type=_h_rt,
                        help="Fixed h_rt for every fragment, still clamped by --max-h-rt. "
                             "Default: tiered alongside slot count.")
    parser.add_argument('--num-procs', default=None, type=_int_or_none,
                        help="Fixed slot count for every fragment. Default: tiered by "
                             "estimated cost.")
    parser.add_argument('--max-h-rt', required=True, type=_h_rt,
                        help="Hard ceiling on -l h_rt, as HH:MM:SS; every tier is clamped "
                             "to this.")
    parser.add_argument('--frag-cost-estimate', default=None,
                        help="TSV from estimate_frag_cost.py. Omit to compute inline from "
                             "--pdb-dir (default).")
    parser.add_argument('--sample-size', default=DEFAULT_SAMPLE_SIZE, type=int,
                        help="Structures sampled for inline cost estimate. "
                             f"Default: {DEFAULT_SAMPLE_SIZE}.")
    parser.add_argument('--estimate-procs', default=10, type=int,
                        help="Worker processes for the inline cost estimate. Default: 10.")
    parser.add_argument('--mem-free', default='2G',
                        help="SGE -l mem_free, PER SLOT; fallback when --num-procs isn't "
                             f"in {sorted(MEM_FREE_PER_SLOT)} (those use "
                             f"{MEM_FREE_PER_SLOT}). Default: 2G.")
    parser.add_argument('--subset-sizes', nargs='+', type=int, default=[1, 2],
                        choices=[1, 2],
                        help="vdG subset sizes to build. Default: 1 2.")
    args = parser.parse_args()
    if args.include_file:
        with open(args.include_file) as handle:
            args.include = list(args.include) + [
                line.strip() for line in handle
                if line.strip() and not line.startswith('#')]
    if args.include_only and not args.include:
        parser.error("--include-only needs at least one --include/--include-file "
                     "fragment; with none it would generate an empty fleet.")
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

    replace = {'$WRAPPER':        WRAPPER_PATH,
               '$LOG_DIR':       args.log_dir,
               '$PDB_DIR':       args.pdb_dir,
               '$OUTPUT_DIR':    args.vdg_lib_dir,
               '$MAX_NUM_CLUS':  str(args.max_num_clus),
               '$SUBSET_SIZES':  ' '.join(str(s) for s in sorted(set(args.subset_sizes))),
               '$SCRATCH':       SCRATCH}

    if args.clear_partial and not args.resume:
        raise SystemExit('[ERROR] --clear-partial only means anything with --resume.')

    if args.min_support is None and not args.include_only:
        raise SystemExit(
            '[ERROR] --min-support is required for a full or resumed build: the min '
            'distinct parent biounits (DR-5) containing a fragment, from the frontier '
            'table. No default -- --min-instances counted a different, non-convertible unit.')

    if not os.path.exists(args.sge_out_dir):
        os.makedirs(args.sge_out_dir)
    visible = [f for f in os.listdir(args.sge_out_dir) if not f.startswith('.')]
    if visible and not args.resume:
        raise FileExistsError(f"{args.sge_out_dir} already has files; refusing to "
                              "overwrite. Use --resume to write only unfinished scripts.")

    if args.include_only or args.resume:
        check_provenance(args.vdg_lib_dir, args)

    frags_dict, support_pooled, frags_meta = load_frags_dict(args.frags_dict)
    check_frags_dict_identity(frags_meta, args)
    check_library_vocabulary(args.vdg_lib_dir, frags_meta.get('key_schema'))
    prepared = prepare_fragments_cached(frags_dict, args.max_size, args.frags_dict)
    candidates = select_fragments(frags_dict, 0, args.max_size, prepared=prepared)
    if args.frag_cost_estimate:
        est, occurrences = read_estimate_tsv(args.frag_cost_estimate)
        _hdr = read_estimate_header(args.frag_cost_estimate)
        check_estimate_header(_hdr, args)
        _scale = _hdr.get('sample_scale')
        sample_scale = float(_scale) if _scale not in (None, 'n/a') else None
        missing = [s for s in candidates if s not in est]
        if missing:
            raise ValueError(
                f"--frag-cost-estimate is missing {len(missing)} candidate(s), e.g. "
                f"{missing[:3]}. Regenerate against the same --max-size.")
    else:
        print(f'Counting {len(candidates)} candidate fragments over '
              f'{args.sample_size} sampled structures...')
        est, occurrences = estimate_fragment_counts(
            candidates, args.pdb_dir, sample_size=args.sample_size,
            num_procs=args.estimate_procs)
        sample_scale = _LAST_SAMPLE_SCALE[0]

    smiles_to_run, aliases = select_fragments(
        frags_dict, args.min_support, args.max_size, support=support_pooled,
        return_aliases=True, prepared=prepared,
        include=args.include, include_only=args.include_only)

    skipped_finished, partial_dirs = [], []
    if args.resume:
        smiles_to_run, skipped_finished, partial_dirs = partition_by_completion(
            smiles_to_run, args.vdg_lib_dir)
        if not smiles_to_run:
            print(f'Resume: all {len(skipped_finished)} selected fragment(s) already '
                  f"have 'Job completed.' in their log. Nothing to submit.")
            return
        if args.clear_partial and partial_dirs:
            check_no_live_partial(partial_dirs, active_sge_job_names())

    if sample_scale is None:
        print('[WARNING] cost estimate has no sample_scale; tiers use the point '
              'estimate and carry ~15% seed-to-seed tier instability. Regenerate to fix.')

    tier_counts = {}
    clamped_top_tier = []
    with open(args.template, 'r') as f:
        template_lines = f.readlines()
    partial_set = set(partial_dirs)
    written_scripts = []
    submission_rows = []

    for smiles in smiles_to_run:
        tier_count = sampling_upper_bound(occurrences[smiles], sample_scale)
        slots, h_rt = resources_for(tier_count, args.max_h_rt,
                                    fixed_num_procs=args.num_procs,
                                    fixed_h_rt=args.h_rt)
        if (tier_count >= TOP_TIER_LOWER_BOUND
                and _h_rt_to_hours(h_rt) < _h_rt_to_hours(TOP_TIER_H_RT)):
            clamped_top_tier.append(smiles)
        per_frag = dict(replace, **{
            '$NUM_PROCS': str(slots), '$RUN_TIME': h_rt,
            '$MEM_FREE': MEM_FREE_PER_SLOT.get(slots, args.mem_free)})
        per_frag['$PRE_RUN'] = (
            clear_partial_snippet(args.vdg_lib_dir, utils.smiles_to_filename(smiles))
            if (args.clear_partial and smiles in partial_set) else '')
        tier_counts[(slots, h_rt)] = tier_counts.get((slots, h_rt), 0) + 1
        output_script(template_lines, smiles, args.sge_out_dir, per_frag)
        script_path = os.path.join(args.sge_out_dir,
                                   utils.smiles_to_filename(smiles) + '.sh')
        written_scripts.append(script_path)
        submission_rows.append((script_path, smiles, slots, h_rt))

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
    with open(manifest, 'w') as handle:
        handle.write('# Submit THESE scripts, in this order. Do not use `qsub *.sh`: a '
                     'resumed or topped-up directory still holds scripts for fragments '
                     'that already finished, and the glob silently redoes them.\n')
        handle.write('order\tscript\tfragment\tslots\th_rt\n')
        for order, (path, smiles, slots, h_rt) in enumerate(submission_rows):
            handle.write(f'{order}\t{path}\t{smiles}\t{slots}\t{h_rt}\n')
    print(f'Submission order written to {manifest}. Submit exactly these, in order:\n'
          f"  awk -F'\\t' 'NR>2 {{print $2}}' {manifest} | xargs -n1 qsub")

    if args.resume:
        print(f'Resume: {len(skipped_finished)} fragment(s) already finished and were '
              f'skipped; wrote {len(written_scripts)} script(s) to {args.sge_out_dir}.')
        if partial_dirs:
            if args.clear_partial:
                print(f'[WARNING] {len(partial_dirs)} fragment(s) have leftover output '
                      f'dirs; --clear-partial will DELETE each before rerunning. e.g. '
                      f'{[utils.smiles_to_filename(s) for s in partial_dirs[:3]]}')
            else:
                print(f'[WARNING] {len(partial_dirs)} fragment(s) have a leftover output '
                      f'dir from a killed job; the wrapper refuses a non-empty output '
                      f'dir, so those jobs will die immediately. Remove them, or rerun '
                      f'with --clear-partial:')
                for smiles in partial_dirs:
                    print(f'  {os.path.join(args.vdg_lib_dir, utils.smiles_to_filename(smiles))}')
    elif args.include_only:
        print(f'Top-up: created scripts for {len(smiles_to_run)} requested '
              f'fragment(s) in {args.sge_out_dir} (count threshold not applied).')
    else:
        print(f'Created scripts for {len(smiles_to_run)} of {len(candidates)} candidate '
              f'fragments (>= {args.min_support} parent biounits'
              f'{f", plus {len(args.include)} requested" if args.include else ""}) '
              f'in {args.sge_out_dir}.')
    print('Resources requested: ' + ', '.join(
        f'{n} fragment(s) at -pe smp {slots}, h_rt {h_rt}'
        for (slots, h_rt), n in sorted(tier_counts.items())))
    slot_hours = {key: n * key[0] * _h_rt_to_hours(key[1])
                  for key, n in tier_counts.items()}
    print('Worst-case slot-hours (ceiling, not a forecast): ' + ', '.join(
        f'{slot_hours[key]:,.0f} at -pe smp {key[0]}/{key[1]}'
        for key in sorted(tier_counts)) +
        f'; TOTAL {sum(slot_hours.values()):,.0f} slot-hours over '
        f'{sum(tier_counts.values())} job(s), peak {sum(k[0] * n for k, n in tier_counts.items()):,} '
        f'slots if every job ran at once.')
    if clamped_top_tier:
        print(f'[WARNING] {len(clamped_top_tier)} fragment(s) at >= {TOP_TIER_LOWER_BOUND} '
              f'estimated occurrences request less than the top tier\'s {TOP_TIER_H_RT}, '
              f'which is itself unbounded.')
    if aliases:
        print(f'Collapsed {len(aliases)} protonation variant(s); wrote {alias_path}.')
        promoted = sorted(r for r in set(aliases.values())
                          if alias_kind(r, dict_keys) == 'promoted')
        if promoted:
            print(f'{len(promoted)} representative(s) are promoted charge-stripped keys '
                  f'absent from the fragment dict (library dirs named for a SMARTS no '
                  f'CCD ligand is drawn with). See kind=promoted in {alias_path} and '
                  f'scripts/lookup_fragment_key.py: {promoted}')

def active_sge_job_names():
    try:
        return {el.text for el in ET.fromstring(subprocess.run(
            ['qstat', '-xml', '-u', getpass.getuser()], capture_output=True, text=True,
            check=True, timeout=30).stdout).iter('JB_name') if el.text}
    except (OSError, subprocess.SubprocessError) as exc:
        raise SystemExit(f'[ERROR] --clear-partial could not query qstat for live jobs '
            f'before deleting partial output ({exc}). Fix or drop --clear-partial.')

def check_no_live_partial(partial_dirs, active):
    live = [s for s in partial_dirs if utils.smiles_to_job_name(s) in active]
    if live:
        raise SystemExit(
            f'[ERROR] --clear-partial would delete output of {len(live)} fragment(s) '
            f'among {len(active)} job(s) qstat shows live for this user. Wait for them '
            f'to finish, or drop --clear-partial.')

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

def clear_partial_snippet(vdg_lib_dir, label):
    target = os.path.join(os.path.abspath(vdg_lib_dir), label)
    return (
        f'# --clear-partial: this fragment did not finish and its directory would\n'
        f'# make the wrapper refuse to start. Remove only that directory, and only\n'
        f'# if it looks like this fragment\'s own output.\n'
        f'RESUME_DIR={shlex.quote(target)}\n'
        f'if [ -n "$RESUME_DIR" ] && [ -d "$RESUME_DIR" ] && '
        f'[ -e "$RESUME_DIR/{label}_log" ]; then\n'
        f'    echo "resume: clearing partial output $RESUME_DIR"\n'
        f'    rm -rf -- "$RESUME_DIR"\n'
        f'fi\n')

def output_script(template_lines, smiles, sge_out_dir, replace):
    local_replace = dict(replace)

    script_name = os.path.join(sge_out_dir, utils.smiles_to_filename(smiles) + '.sh')
    local_replace['$SMILES'] = f'"{smiles}"'
    local_replace['$CG'] = f'"{utils.smiles_to_filename(smiles)}"'
    local_replace['$JOB_NAME'] = utils.smiles_to_job_name(smiles)
    local_replace.setdefault('$PRE_RUN', '')

    lines = template_lines

    with open(script_name, 'w') as f:
        start_copy = False
        pattern = '|'.join(re.escape(key) for key in local_replace.keys())
        for line in lines:
            if line.startswith('#!/bin/bash'):
                start_copy = True
            if start_copy:
                line = re.sub(pattern, lambda m: local_replace[m.group(0)], line)
                f.write(line)

if __name__ == '__main__':
    main()
