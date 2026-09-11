'''
Generate per-fragment SGE job scripts for running the vdG generation pipeline.

For each fragment in database_frags_dict.pkl that meets the counts and size thresholds,
this script writes one SGE shell script that calls vdg_generation_wrapper.py for that
fragment. The scripts are written to --sge-out-dir and can be submitted with:
    qsub <script>.sh

The generated scripts use `#$ -cwd` plus the relative path
ligand_vdgs/generate_vdgs/vdg_generation_wrapper.py, so `qsub` must be run from the
repo root.

All defaults are set for the Wynton HPC cluster. If you are running on a different
cluster or with a custom database, override the relevant flags (see usage below).
Paths that typically need to change for a custom setup:
    --pdb-dir      directory of prepared, protonated PDB files (your parent database)
    --probe-dir    directory of Probe output files for those PDBs
    --vdg-lib-dir  where to write the vdG library output
    --log-dir      where SGE should write job logs


Usage:
    python ligand_vdgs/generate_vdgs/make_sge_scripts_for_frags.py \\
        --pdb-dir   <path/to/your/pdb_database/> \\
        --probe-dir <path/to/your/probe_output/> \\
        --vdg-lib-dir <path/to/output/vdg_library/> \\
        --log-dir   <path/to/sge_logs/> \\
        --max-h-rt  <HH:MM:SS>
'''

import os
import re
import json
import fcntl
import argparse
import pickle as pkl
from ligand_vdgs.functions import utils
from ligand_vdgs.functions.utils import _int_or_none, file_sha256
from ligand_vdgs.functions.vdg_npz_utils import load_fragment_aliases
from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
    alias_kind, fragment_dict_keys, prepare_fragments, resolve_include_fragments,
    select_fragments, write_fragment_aliases)
from ligand_vdgs.generate_vdgs.estimate_frag_cost import (
    DEFAULT_SAMPLE_SIZE, estimate_fragment_counts, read_estimate_tsv,
    read_estimate_header, sampling_upper_bound, _LAST_SAMPLE_SCALE)


# Slot/wall-time tiers, keyed by the number of structures a fragment
# occurs in (the first return of estimate_frag_cost.estimate_fragment_counts).
# Validated against a
# real 317-fragment build: this predictor reaches log-log correlation 0.63 with
# wall time where the fragment's CCD ligand count reaches only 0.28, and the
# bins are monotone with a hard ceiling in the cheap ones --
#
#   est. structures     n     median      p90        max     (-m CAPPED build)
#   0 - 1000          170     0.04 h     0.11 h     1.2 h
#   1000 - 5000        92     0.18 h     0.68 h     2.3 h
#   5000+              55     1.48 h    12.93 h    24.3 h
#
# so the small tier is safe rather than merely typical. Handing every fragment
# 20 slots makes ~75% of the fleet queue for an allocation it cannot use.
#
# THAT BUILD WAS CAPPED (-m), so every number above is a LOWER BOUND on wall
# time, not a bound: -m truncates clustering, the one phase that scales worst,
# while leaving mining and streaming faithful. It is evidence for the ranking
# the predictor produces and for the cheap tiers (whose jobs finish before the
# cap can bite); it is not evidence for the top tier's h_rt. The uncapped
# profile runs (run_uncapped_frags.sh, 20 slots, subset sizes 1 2) say how far
# off it is:
#
#   fragment                              est structs  autos    uncapped wall
#   cnnnn                                        219      1           0.06 h
#   [C;!R][C;!R](=[O;!R])[O;!R]               15,635      2           0.62 h
#   [C;!R][C;!R][C;!R][O;!R]                  22,020      1           0.99 h
#   [C;r5][C;r5]([C;r5])[O;!R]                21,451      2           4.59 h
#   cccnc                                     35,380      1           8.37 h
#   [O;!R]=[P;!R]([O;!R])([O;!R])[O;!R]       22,151     24    >19.4 h, UNFINISHED
#
# The phosphate is the case the top tier exists for -- structure count in the
# middle of the set, 24 CG automorphisms -- and it was still in
# clus_and_deduplicate_vdgs at 19.4 h when this was written, so the top tier's
# worst case has never been measured to completion and no h_rt here is known to
# cover it. Note also that automorphism count, which the tiers do not use, moves
# wall time harder than structure count does at the top end (phosphate vs. cccnc).
# That is why the last tier asks for Wynton's maximum rather than a number
# derived from the table: for an unbounded worst case the only defensible
# request is all of it. Clamping it with --max-h-rt is a real risk, not a
# formality -- see resources_for.
#
# The first tier exists to reach Wynton's short queue, which admits jobs
# requesting h_rt <= 30 min and gives them a much larger slice of the cluster --
# worth having when a whole build has to land inside a deadline. It is bounded by
# the same 317-fragment build, re-estimated with the current predictor:
#
#   est. structures     n     median     p99       max      (wall, on 20 slots)
#   0 - 750           159     2.5 min   8.2 min   8.8 min
#   750 - 1000         18     4.6 min             19.9 min
#
# so the tier stops just under a real cliff -- the first job past 10 min sits at
# est 853. Sizing the slots is the safety margin, not the threshold: the maxima
# are flat across the whole tier (7.0 min already at est < 50), so tightening the
# cutoff buys almost nothing, while slots scale the wall time directly. The
# pipeline is pool-parallel with only ~1 min of fixed overhead (the floor at
# est < 10), so wall time on S slots is at worst t_20 * 20/S; 16 slots bounds the
# tier at ~11 min against a 30 min limit, a ~2.7x margin on a deliberately
# pessimistic model. That is why the cheapest tier asks for *more* slots than the
# next one: the slots buy wall-clock compression to fit the window, not
# throughput. Do not trade them back for a wider threshold -- being killed at the
# limit costs a whole resubmit, while a fragment left out of the tier still runs.
#
# A short-tier job that overruns is killed at 30 min and leaves a partial
# fragment directory, which is recoverable but not free: it is invisible except
# as a missing 'Job completed.' in <cg_label>/<cg_label>_log (what
# Frags.check_vdg_job_status reads). After a build, resubmit any fragment whose
# log lacks that line. Pass --no-short-queue to skip the tier entirely.
#
# (est_structures_upper_bound, slots, h_rt)
SHORT_QUEUE_TIER = (750, 16, '0:29:00')
RESOURCE_TIERS = (
    SHORT_QUEUE_TIER,
    (1000, 4, '6:00:00'),
    (5000, 8, '24:00:00'),
    (float('inf'), 20, '336:00:00'),   # 2 weeks is the Wynton maximum
)
TOP_TIER_LOWER_BOUND = RESOURCE_TIERS[-2][0]   # est. structures entering the top tier
TOP_TIER_H_RT = RESOURCE_TIERS[-1][2]
# Longest uncapped profile run that actually finished (cccnc, 35,380 structures,
# 20 slots). Not a bound on the tier -- the phosphate run passed 19.4 h without
# finishing -- just the largest number anyone has measured end to end.
LONGEST_FINISHED_UNCAPPED_H = 8.4


def _h_rt_to_hours(h_rt):
    hours, minutes, seconds = (int(part) for part in h_rt.split(':'))
    return hours + minutes / 60 + seconds / 3600


def _h_rt(value):
    """argparse type for an SGE wall-clock limit. Rejects bare hours ('36') here
    rather than in _h_rt_to_hours, where it would surface as an unpack error after
    the estimate pass has already run."""
    match = re.fullmatch(r'(\d+):(\d{1,2}):(\d{1,2})', value)
    if not match or int(match.group(2)) > 59 or int(match.group(3)) > 59:
        raise argparse.ArgumentTypeError(
            f"expected HH:MM:SS (e.g. 36:00:00), got {value!r}")
    return value


def resources_for(est_structures, max_h_rt, fixed_num_procs=None, fixed_h_rt=None,
                  short_queue=True):
    """(slots, h_rt) for a fragment estimated to occur in `est_structures` structures.

    `max_h_rt` is a hard ceiling, not a default: a request above what the queue
    will actually run (a maintenance window, say) is not scheduled at all, so
    clamping is better than a job that never starts. A clamped job may be killed
    mid-run, which is recoverable; an unscheduled one is not.

    "Recoverable" holds only where a rerun can finish. For the top tier it may
    not: resubmitting under the same ceiling reproduces the same kill, and that
    tier's cost is unmeasured (see the phosphate row above). main() reports how
    many fragments were clamped out of it.
    """
    tiers = RESOURCE_TIERS if short_queue else RESOURCE_TIERS[1:]
    for upper, slots, h_rt in tiers:
        if est_structures < upper:
            break
    if fixed_num_procs is not None:
        slots = int(fixed_num_procs)
    if fixed_h_rt is not None:
        h_rt = fixed_h_rt
    if _h_rt_to_hours(h_rt) > _h_rt_to_hours(max_h_rt):
        h_rt = max_h_rt
    return slots, h_rt


PROVENANCE_FILENAME = 'library_provenance.json'


def provenance_path(vdg_lib_dir):
    return os.path.join(vdg_lib_dir, PROVENANCE_FILENAME)


def _selection_inputs(args):
    """The inputs that decide which key a fragment request resolves to.

    Only these belong here. --min-instances is deliberately absent: a top-up
    skips the threshold entirely, so a different value is not a conflict. The
    fragment dict and --max-size are different -- prepare_fragments derives the
    representative/alias mapping from them, so a change can silently resolve the
    same SMARTS to a different key and write a directory the rest of the library
    does not match.
    """
    return {'frags_dict': os.path.abspath(args.frags_dict),
            'frags_dict_sha256': file_sha256(args.frags_dict),
            'max_size': args.max_size}


def write_provenance(vdg_lib_dir, args):
    record = _selection_inputs(args)
    record['min_instances'] = args.min_instances      # recorded, not compared
    # Also recorded, not compared: the parent db is machine-local, so a top-up run
    # elsewhere legitimately differs. Buckets store this per file too; here it is
    # one editable place saying what the library was built against, for a copy
    # whose reader must point $PARENT_PDBS_DIR somewhere local.
    record['parent_pdb_dir'] = os.path.abspath(args.pdb_dir)
    os.makedirs(vdg_lib_dir, exist_ok=True)
    # Written via a temp file: a crash mid-dump would otherwise leave truncated
    # JSON that makes every later check_provenance raise on json.load, which
    # reads as a corrupted library rather than an interrupted write.
    final = provenance_path(vdg_lib_dir)
    tmp = final + '.tmp'
    with open(tmp, 'w') as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write('\n')
    os.replace(tmp, final)


def check_provenance(vdg_lib_dir, args):
    """Refuse a top-up whose fragment resolution would not match the library's.

    A missing file means the library predates this check, which cannot be
    verified either way -- warn rather than block, so an older library stays
    usable.
    """
    path = provenance_path(vdg_lib_dir)
    if not os.path.isfile(path):
        print(f'[WARNING] {path} is missing, so the fragment dict this library '
              f'was built from cannot be verified against the one being used '
              f'now. If they differ, a requested SMARTS can resolve to a key '
              f'the rest of the library does not use.')
        return
    with open(path) as handle:
        recorded = json.load(handle)
    current = _selection_inputs(args)
    differing = {k: (recorded.get(k), v) for k, v in current.items()
                 if recorded.get(k) != v}
    # A moved-but-identical dict is not a conflict; the hash is what matters.
    differing.pop('frags_dict', None)
    if differing:
        detail = '; '.join(f'{k}: library has {was!r}, this run has {now!r}'
                           for k, (was, now) in sorted(differing.items()))
        raise SystemExit(
            f'ERROR: --include-only cannot add to {vdg_lib_dir}: the inputs that '
            f'decide how a fragment SMARTS resolves have changed since it was '
            f'built ({detail}). The same request can now resolve to a different '
            f'fragment key than the rest of the library uses. Either rebuild the '
            f'library, or restore the recorded inputs.')


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create SGE submission scripts for vdG generation.")
    parser.add_argument('--frags-dict', default='resources/database_frags_dict.pkl',
                        help="Path to database_frags_dict.pkl.")
    parser.add_argument('--template', default='resources/frag_sge_template.sh',
                        help="Path to SGE job script template.")
    parser.add_argument('--include', nargs='*', default=[], metavar='SMARTS',
                        help="Fragment SMARTS to build regardless of their "
                             "estimated count. Matched up to equivalent atom "
                             "ordering, and a charged variant resolves to its "
                             "representative. Use for fragments you need whether "
                             "or not the sampling estimate happened to clear the "
                             "threshold.")
    parser.add_argument('--include-file', default=None,
                        help="File of fragment SMARTS to --include, one per line "
                             "('#' comments and blank lines ignored).")
    parser.add_argument('--include-only', action='store_true',
                        help="Build scripts for the --include fragments alone, "
                             "skipping the count threshold. This is the top-up "
                             "mode: it adds fragments to an existing library "
                             "without rebuilding what is already there.")
    parser.add_argument('--min-instances', default=250, type=int,
                        help="Min candidate vdG sites for a fragment to be built: CG "
                             "occurrences, i.e. SMARTS matches summed over every ligand "
                             "copy in the parent PDB db. Default: 250.")
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
                        default='/wynton/group/degradolab/skt/docking/databases/prepwizard_BioLiP2/',
                        help="Path to parent PDB database.")
    parser.add_argument('--probe-dir',
                        default='/wynton/group/degradolab/skt/docking/databases/probe_output/',
                        help="Path to Probe output directory.")
    parser.add_argument('--max-num-clus', default=None, type=_int_or_none,
                        help="Max number of vdGs to cluster per subset. Default: no limit.")
    parser.add_argument('--h-rt', default=None, type=_h_rt,
                        help="Fixed SGE wall-clock limit for every fragment, still "
                             "clamped by --max-h-rt. Default: unset, meaning it is "
                             "tiered alongside the slot count.")
    parser.add_argument('--num-procs', default=None, type=_int_or_none,
                        help="Fixed slot count for every fragment. Default: unset, "
                             "meaning slots are tiered per fragment by estimated cost.")
    parser.add_argument('--no-short-queue', dest='short_queue', action='store_false',
                        help="Do not put the cheapest fragments in the "
                             f"{SHORT_QUEUE_TIER[2]} tier that qualifies for Wynton's "
                             "short queue. Default: the tier is on, which trades a small "
                             "risk of a killed, resubmittable job for much shorter "
                             "queueing on roughly half the fleet.")
    parser.add_argument('--max-h-rt', required=True, type=_h_rt,
                        help="Hard ceiling on -l h_rt, as HH:MM:SS. Required: a request "
                             "longer than the queue will run (a maintenance window, say) "
                             "is never scheduled, so every tier is clamped to this.")
    parser.add_argument('--frag-cost-estimate', default=None,
                        help="TSV from estimate_frag_cost.py, covering every candidate "
                             "fragment. Omit to compute it inline from --pdb-dir, which "
                             "is the default: ~1 min for 5713 fragments over a "
                             "3000-structure sample at 16 procs, and it keeps a build "
                             "dependent on nothing but the parent db and the fragment dict.")
    parser.add_argument('--sample-size', default=DEFAULT_SAMPLE_SIZE, type=int,
                        help="Structures sampled when estimating cost inline. "
                             f"Default: {DEFAULT_SAMPLE_SIZE}.")
    parser.add_argument('--estimate-procs', default=10, type=int,
                        help="Worker processes for the inline cost estimate. Default: 10.")
    parser.add_argument('--mem-free', default='4G',
                        help="SGE -l mem_free, PER SLOT under -pe smp. Default: 4G. "
                             "Size this from a finished job's qacct maxvmem divided by "
                             "its slot count, NOT from the profile's rss_peak_child_mb: "
                             "that counter is the largest single child's peak, not the "
                             "sum over concurrent children. Measured worst case over six "
                             "uncapped 20-slot runs is 24.4G total = 1.22G/slot.")
    parser.add_argument('--scratch', default='20G',
                        help="SGE -l scratch (not per slot). Default: 20G. Set this from "
                             "a profiled run's scratch_peak_total_mb; measured worst case "
                             "is 5.4G, and scratch stays under 0.152 MB/structure across "
                             "the profiled set. Oversizing costs scheduling latency, since "
                             "it narrows the set of eligible nodes.")
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

    for path, flag in [(args.frags_dict, '--frags-dict'), (args.template, '--template')]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"{flag} path does not exist: {path}")
    # These are interpolated into every generated script (--log-dir also into
    # `#$ -o`), so a missing one fails at qsub time across all ~850 jobs rather
    # than here. --pdb-dir additionally drives fragment *selection* via the
    # inline cost estimate, so a wrong-but-existing parent db silently changes
    # which fragments get built -- validate its existence at minimum.
    for path, flag in [(args.log_dir, '--log-dir'), (args.pdb_dir, '--pdb-dir'),
                       (args.probe_dir, '--probe-dir')]:
        if not os.path.isdir(path):
            raise NotADirectoryError(f"{flag} directory does not exist: {path}")

    replace = {'$LOG_DIR':       args.log_dir,
               '$PDB_DIR':       args.pdb_dir,
               '$PROBE_DIR':     args.probe_dir,
               '$OUTPUT_DIR':    args.vdg_lib_dir,
               '$MAX_NUM_CLUS':  str(args.max_num_clus),
               '$SUBSET_SIZES':  ' '.join(str(s) for s in sorted(set(args.subset_sizes))),
               '$MEM_FREE':      args.mem_free,
               '$SCRATCH':       args.scratch}

    if not os.path.exists(args.sge_out_dir):
        os.makedirs(args.sge_out_dir)
    visible = [f for f in os.listdir(args.sge_out_dir) if not f.startswith('.')]
    if visible:
        raise FileExistsError(f"Output directory {args.sge_out_dir} already has "
                              "files. Terminating to prevent overwriting.")

    # Before the expensive dict walk, and before any script is written: a
    # mismatch here invalidates every fragment key this run would produce.
    if args.include_only:
        check_provenance(args.vdg_lib_dir, args)

    with open(args.frags_dict, 'rb') as f:
        frags_dict = pkl.load(f)
    # Shared with the scheduler-agnostic path so the two cannot select
    # different fragment sets. Returns a sorted list, so runs are reproducible.
    # One estimate pass serves both jobs: which fragments are worth building, and
    # how many slots each needs. Candidates are everything that survives the size
    # and solvent filters, since the threshold is applied against the estimate.
    prepared = prepare_fragments(frags_dict, args.max_size)
    candidates = select_fragments(frags_dict, 0, args.max_size, prepared=prepared)
    if args.frag_cost_estimate:
        est, occurrences = read_estimate_tsv(args.frag_cost_estimate)
        # The membership check below only catches keys that vanished. A dict whose
        # *grouping* changed keeps every key present while silently reassigning
        # which of them is a representative, so compare the recorded identity too.
        _hdr = read_estimate_header(args.frag_cost_estimate)
        _recorded = _hdr.get('frags_dict_sha256')
        if _recorded is None:
            print(f'[WARNING] {args.frag_cost_estimate} records no frags_dict hash '
                  f'(written before that was added), so it cannot be verified '
                  f'against --frags-dict. Regenerate it if the dict has changed.')
        elif _recorded != file_sha256(args.frags_dict):
            print(f'[WARNING] {args.frag_cost_estimate} was computed from a '
                  f'different {os.path.basename(args.frags_dict)} (sha256 '
                  f'{_recorded[:12]}... vs {file_sha256(args.frags_dict)[:12]}...). '
                  f'Counts still cover every current candidate or the check below '
                  f'would fail, but fragment grouping may have changed since. '
                  f'Regenerate with estimate_frag_cost.py to be sure.')
        _hdr_max_size = _hdr.get('max_size')
        if _hdr_max_size is not None and int(_hdr_max_size) != args.max_size:
            print(f'[WARNING] {args.frag_cost_estimate} was computed at --max-size '
                  f'{_hdr_max_size}, this run uses {args.max_size}.')
        _scale = _hdr.get('sample_scale')
        sample_scale = float(_scale) if _scale is not None else None
        missing = [s for s in candidates if s not in est]
        if missing:
            raise ValueError(
                f"--frag-cost-estimate is missing {len(missing)} candidate fragment(s), "
                f"e.g. {missing[:3]}. Regenerate it against the same --max-size.")
    else:
        print(f'Counting {len(candidates)} candidate fragments over '
              f'{args.sample_size} sampled structures...')
        est, occurrences = estimate_fragment_counts(
            candidates, args.pdb_dir, sample_size=args.sample_size,
            num_procs=args.estimate_procs)
        sample_scale = _LAST_SAMPLE_SCALE[0]

    # Two quantities from one pass, and they are not interchangeable: selection asks
    # how many vdG sites a fragment offers (occurrences), while the resource tiers
    # were calibrated against how many structures the job must read (est).
    smiles_to_run, aliases = select_fragments(
        frags_dict, args.min_instances, args.max_size, return_aliases=True,
        instance_counts=occurrences, prepared=prepared,
        include=args.include, include_only=args.include_only)

    # --include bypasses the occurrence threshold, so an included fragment is
    # typically one the sample barely saw -- est 0 puts it in the short-queue
    # tier, where an expensive fragment (many automorphisms, which the tiers do
    # not model) is killed at h_rt on every resubmit. Deny it the short tier.
    included_reps = set(resolve_include_fragments(args.include, prepared)[0].values())

    # Tier on the upper end of the sampling interval, not the point estimate --
    # see sampling_upper_bound. Absent a recorded sample_scale (a TSV written
    # before the header carried it) there is nothing to widen, so the point
    # estimate is used and the run says so rather than silently under-sizing.
    if sample_scale is None:
        print('[WARNING] the cost estimate records no sample_scale, so resource '
              'tiers use the point estimate and carry the full seed-to-seed tier '
              'instability (~15% of fragments). Regenerate it to remove this.')

    # Each scheduler passes only the fragment SMARTS. The wrapper/core derive
    # the exact automorphisms identically at execution time.
    tier_counts = {}
    clamped_top_tier = []
    with open(args.template, 'r') as f:  # read once; ~1 NFS read, not one per fragment
        template_lines = f.readlines()
    for smiles in smiles_to_run:
        tier_count = sampling_upper_bound(est[smiles], sample_scale)
        slots, h_rt = resources_for(tier_count, args.max_h_rt,
                                    fixed_num_procs=args.num_procs,
                                    fixed_h_rt=args.h_rt,
                                    short_queue=(args.short_queue and
                                                 smiles not in included_reps))
        # The top tier is the one whose cost was never measured to completion,
        # so a ceiling that cuts into it is worth saying out loud rather than
        # applying silently.
        if (est[smiles] >= TOP_TIER_LOWER_BOUND
                and _h_rt_to_hours(h_rt) < _h_rt_to_hours(TOP_TIER_H_RT)):
            clamped_top_tier.append(smiles)
        per_frag = dict(replace, **{'$NUM_PROCS': str(slots), '$RUN_TIME': h_rt})
        tier_counts[(slots, h_rt)] = tier_counts.get((slots, h_rt), 0) + 1
        output_script(template_lines, smiles, args.sge_out_dir, per_frag)

    # Written into the library root, not next to the scripts: a consumer holding
    # a charged fragment name resolves it against the library it is reading.
    alias_path = os.path.join(args.vdg_lib_dir, 'fragment_aliases.tsv')
    dict_keys = fragment_dict_keys(frags_dict)
    if args.include_only:
        # Top-up: this run knows only about the fragments it was asked for, so
        # overwriting would drop every alias the original build recorded and make
        # those charged variants unresolvable against a library that still holds
        # their vdGs.
        #
        # Locked because this is a read-modify-write: two concurrent top-ups
        # would otherwise each read the pre-existing file and the second writer
        # would drop the first's rows. The lock file is separate from the target
        # so the lock survives write_fragment_aliases replacing it.
        os.makedirs(args.vdg_lib_dir, exist_ok=True)
        with open(alias_path + '.lock', 'w') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            merged = dict(load_fragment_aliases(args.vdg_lib_dir))
            merged.update(aliases)
            write_fragment_aliases(alias_path, merged, dict_keys)
        aliases = merged
    else:
        write_fragment_aliases(alias_path, aliases, dict_keys)
        write_provenance(args.vdg_lib_dir, args)

    if args.include_only:
        print(f'Top-up: created scripts for {len(smiles_to_run)} requested '
              f'fragment(s) in {args.sge_out_dir} (count threshold not applied).')
    else:
        print(f'Created scripts for {len(smiles_to_run)} of {len(candidates)} candidate '
              f'fragments (>= {args.min_instances} CG occurrences'
              f'{f", plus {len(args.include)} requested" if args.include else ""}) '
              f'in {args.sge_out_dir}.')
    print('Resources requested: ' + ', '.join(
        f'{n} fragment(s) at -pe smp {slots}, h_rt {h_rt}'
        for (slots, h_rt), n in sorted(tier_counts.items())))
    short = sum(n for (_, h_rt), n in tier_counts.items()
                if _h_rt_to_hours(h_rt) <= 0.5)
    if short:
        print(f'{short} of these qualify for the short queue (h_rt <= 30 min). If one '
              'overruns it is killed with a partial fragment directory: after the build, '
              "resubmit any fragment whose log lacks 'Job completed.'")
    if clamped_top_tier:
        print(f'[WARNING] {len(clamped_top_tier)} fragment(s) at >= {TOP_TIER_LOWER_BOUND} '
              f'estimated structures request less than the top tier\'s {TOP_TIER_H_RT}. '
              f'Nothing bounds that tier: the longest uncapped profile run that finished '
              f'took {LONGEST_FINISHED_UNCAPPED_H} h, but the phosphate '
              f'(24 CG automorphisms) passed 19.4 h without finishing. A fragment that '
              f'cannot finish under this ceiling is killed the same way on every resubmit, '
              f"so check for 'Job completed.' and rerun those under a longer one. "
              f'e.g. {clamped_top_tier[:3]}')
    if aliases:
        print(f'Collapsed {len(aliases)} protonation variant(s); wrote {alias_path}.')
        promoted = sorted(r for r in set(aliases.values())
                          if alias_kind(r, dict_keys) == 'promoted')
        if promoted:
            print(f'{len(promoted)} representative(s) are promoted charge-stripped keys '
                  f'absent from the fragment dict; their library directories are named '
                  f'for a SMARTS no CCD ligand is drawn with. See kind=promoted in '
                  f'{alias_path} and scripts/lookup_fragment_key.py: {promoted}')


def output_script(template_lines, smiles, sge_out_dir, replace):
    # Per-fragment copy prevents placeholder values leaking between scripts.
    local_replace = dict(replace)

    script_name = os.path.join(sge_out_dir, utils.smiles_to_filename(smiles) + '.sh')
    local_replace['$SMILES'] = f'"{smiles}"'
    local_replace['$CG'] = f'"{utils.smiles_to_filename(smiles)}"'
    local_replace['$JOB_NAME'] = utils.smiles_to_job_name(smiles)  # SGE job names: # truncates directives

    lines = template_lines

    with open(script_name, 'w') as f:
        start_copy = False
        pattern = '|'.join(re.escape(key) for key in local_replace.keys())
        for line in lines:
            # skip header; wait until #!/bin/bash
            if line.startswith('#!/bin/bash'):
                start_copy = True
            if start_copy:
                # replace all placeholders in a single pass to avoid order-dependent bugs
                line = re.sub(pattern, lambda m: local_replace[m.group(0)], line)
                f.write(line)


if __name__ == '__main__':
    main()
