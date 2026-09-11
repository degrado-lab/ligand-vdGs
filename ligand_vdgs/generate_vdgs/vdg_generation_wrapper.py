import sys
import os
import time
import subprocess
import argparse
import shutil
import shlex
import uuid
import concurrent.futures
from ligand_vdgs.generate_vdgs.clus_and_deduplicate_vdgs import EXIT_NO_VDGS
from ligand_vdgs.functions.compute_profile import ComputeProfile
from ligand_vdgs.functions import sasa
from ligand_vdgs.functions.utils import (convert_time_elapsed, identify_mol_automorphisms,
                   mol_from_fragment, set_up_outdir, _int_or_none, smiles_to_filename,
                   aromatic_h_constrained_atoms)

HERE = os.path.dirname(os.path.abspath(__file__))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--smarts', type=str, required=True, help="SMARTS pattern.")
    parser.add_argument('-c', '--cg', type=str, help="CG label; defaults to SMARTS.")
    parser.add_argument('-p', "--pdb-dir", type=str, required=True, help="PDB database directory.")
    parser.add_argument('-o', "--out-dir", type=str, required=True, help="Output directory.")
    parser.add_argument('-m', "--max-num-vdgs-to-clus", default=None, type=_int_or_none,
                        help="Max PDBs to cluster per AA composition.")
    parser.add_argument('--num-procs', type=int, default=10, help="Number of processes.")
    # On by default. Measured cost is ~1 s per fragment run (three os.walk passes
    # over the scratch tree; 0.34 s for the largest shard count observed, 33,464
    # files) against runs of 212 s to 30,134 s, and every other profile call
    # short-circuits when disabled. The information is not reconstructable after
    # the fact, and when this was opt-in the generated SGE scripts never passed
    # the flag, so a whole production fleet could finish having recorded nothing.
    parser.add_argument('--no-profile-compute', dest='profile_compute',
        action='store_false', default=True,
        help="Skip the <cg>_compute_profile.json sidecar of per-phase timings and "
             "counters. Profiling is on by default and costs ~1 s per run.")
    parser.add_argument('--subset-sizes', nargs='+', type=int, default=[1, 2], choices=[1, 2],
                        help="vdG subset sizes (default: 1 2).")
    return parser.parse_args()


def main():
    args = parse_args()
    smarts_raw = args.smarts.strip('"')
    pdb_dir = args.pdb_dir
    out_dir = args.out_dir

    smarts_mol = mol_from_fragment(smarts_raw)
    if smarts_mol is None:
        raise ValueError(f"Could not parse --smarts value as SMARTS: {smarts_raw!r}")
    # Reported here for the log only. The group clustering actually uses is derived
    # again from --cg-smarts in clus_and_deduplicate_vdgs.py, which additionally
    # checks it against the elements of every mined CG before any of it is used.
    num_automorphisms = len(identify_mol_automorphisms(smarts_mol))
    max_num_to_clus = args.max_num_vdgs_to_clus
    num_procs = args.num_procs
    subset_sizes = sorted(set(args.subset_sizes))

    main_script_start = time.time()
    profile = ComputeProfile(enabled=args.profile_compute, smarts=smarts_raw,
                             num_procs=num_procs, subset_sizes=subset_sizes,
                             num_automorphisms=num_automorphisms)

    # The -c label doubles as a path component: smarts_to_cgs.py writes
    # f'{cg}_matches.pkl' from whatever -c it is given, so a raw SMILES
    # containing "/" (e.g. C/C=C/O, which is what -c defaults to when omitted)
    # would point at a nonexistent nested subdir. Encode it once here and pass
    # the *encoded* label to every subprocess, so both sides agree. Safe as a
    # chemistry key too: the only raw -c consumer is the cg_atoms[cg] lookup in
    # clus_and_deduplicate_vdgs.py, whose keys ('ccn', 'gn', 'coo') contain no
    # / or \ and so survive the encoding unchanged.
    cg = smiles_to_filename(args.cg.strip('"')) if args.cg else smiles_to_filename(smarts_raw)
    out_dir = os.path.join(out_dir, cg)

    smarts_to_cg_script = os.path.join(HERE, "..", "..", "external", "vdG-miner",
                                        "vdg_miner", "programs", "smarts_to_cgs.py")
    gen_environments_script = os.path.join(HERE, "..", "..", "external", "vdG-miner",
                                            "vdg_miner", "programs", "generate_environments.py")

    # Checked before set_up_outdir claims the output directory, not after. A typo in
    # --pdb-dir or a missing submodule otherwise surfaces as a CalledProcessError from
    # smarts_to_cgs.py, and the obvious response -- fix the flag and rerun -- then
    # fails a second time because the directory is no longer empty.
    for path, flag in [(pdb_dir, '--pdb-dir')]:
        if not os.path.isdir(path):
            raise NotADirectoryError(f"{flag} directory does not exist: {path}")
    for path, what in [(smarts_to_cg_script, 'smarts_to_cgs.py'),
                       (gen_environments_script, 'generate_environments.py')]:
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"Missing vdG-miner script {what} at {path}. The external/vdG-miner "
                "submodule is probably not checked out (git submodule update --init).")

    set_up_outdir(out_dir)

    # Set logfile directly in cg outdir
    logfile = os.path.join(out_dir, f"{cg}_log")

    write_out_commandline_params(logfile, smarts_raw, cg, pdb_dir, out_dir,
                                 num_automorphisms, max_num_to_clus=max_num_to_clus,
                                 num_procs=num_procs, subset_sizes=subset_sizes)

    # Warn rather than fail: the query is legal and someone may want exactly one
    # tautomer. It is a mistake often enough to be worth saying out loud.
    h_pinned = aromatic_h_constrained_atoms(smarts_raw)
    if h_pinned:
        warning = (
            f"[WARNING] --smarts pins a hydrogen count on an aromatic atom "
            f"({', '.join(h_pinned)}). Automorphism detection normalizes H "
            f"across an aromatic component, but SMARTS matching does not, so "
            f"this mines only the tautomer as modeled and will undercount "
            f"badly. Drop the H constraint (e.g. 'n' for '[nH]') unless the "
            f"single tautomer is the point.\n")
        print(warning, file=sys.stderr)
        with open(logfile, 'a') as _log:
            _log.write(warning)

    tmp_root = None
    for tmp_candidate in [os.environ.get("TMPDIR"), "/scratch", "/tmp"]:
        if tmp_candidate and os.path.isdir(tmp_candidate):
            tmp_root = tmp_candidate
            break

    if tmp_root is None:
        raise RuntimeError("No temporary directory available (TMPDIR, /scratch, /tmp not found).")

    out_leaf = os.path.basename(os.path.normpath(out_dir))
    environments_out_root = os.path.join(
        tmp_root, f"vdg_environments_{cg}_{out_leaf}_{uuid.uuid4().hex[:8]}")
    os.makedirs(environments_out_root, exist_ok=True)

    try:
        # ----- Run smarts_to_cg.py -----
        smarts_to_cg_cmd = (
            f'python {shlex.quote(smarts_to_cg_script)} '
            f'-s {shlex.quote(smarts_raw)} -c {shlex.quote(cg)} '
            f'-p {shlex.quote(pdb_dir)} -o {shlex.quote(out_dir)} '
            f'-l {shlex.quote(logfile)} -n {num_procs}')

        with profile.phase('smarts_to_cgs'):
            subprocess.run(smarts_to_cg_cmd, shell=True, check=True)

        # Check if match output exists and is non-empty (indicates successful matches).
        # smarts_to_cgs.py (external) writes this file as f'{cg}_matches.pkl' from
        # the same -c value passed above, which is why both use the encoded label.
        match_pkl = os.path.join(out_dir, f'{cg}_matches.pkl')
        if not os.path.exists(match_pkl) or os.path.getsize(match_pkl) == 0:
            # A zero-match run still writes its profile sidecar and completion
            # marker before returning, or it is indistinguishable from a crash.
            print('No ligands contain the specified SMARTS pattern.')
            _finish(profile, logfile, main_script_start, out_dir, cg)
            return

        # ----- Mine contact-defined environments (no unused fingerprints) -----
        environments_cmd = (
            f'python {shlex.quote(gen_environments_script)} '
            f'-c {shlex.quote(cg)} -m {shlex.quote(match_pkl)} '
            f'-p {shlex.quote(pdb_dir)} -o {shlex.quote(environments_out_root)}')
        gen_environments_start = time.time()

        with profile.phase('mine_environments'):
            with concurrent.futures.ProcessPoolExecutor(max_workers=num_procs) as executor:
                futures = [executor.submit(run_gen_environments, i, num_procs, environments_cmd)
                           for i in range(num_procs)]
                try:
                    for future in concurrent.futures.as_completed(futures):
                        future.result()
                except Exception:
                    # One dead shard means the run is over, so the queued ones are
                    # pure waste -- and the `finally` below deletes
                    # environments_out_root, which any shard still running is
                    # writing into. Drop what has not started and wait out what has.
                    executor.shutdown(wait=True, cancel_futures=True)
                    raise

        gen_environments_elapsed = time.time() - gen_environments_start
        hours, minutes, seconds = convert_time_elapsed(gen_environments_elapsed)

        # One walk for both the shard count and the scratch footprint.
        num_environment_shards = 0
        environments_bytes = 0
        for root, _dirs, files in os.walk(environments_out_root):
            for fname in files:
                if fname.endswith('.jsonl'):
                    num_environment_shards += 1
                try:
                    environments_bytes += os.path.getsize(os.path.join(root, fname))
                except OSError:
                    pass

        with open(logfile, 'a') as f:
            f.write(f"Completed generate_environments.py in {hours} h, {minutes} mins, and {seconds} secs.\n")
            f.write(f"\t{num_environment_shards} environment shards generated.\n")
        profile.set('environment_shards', num_environment_shards)
        profile.record_peak_rss('mine.')
        # The environment shards are not deleted until the run ends, so they sit
        # on scratch alongside everything clustering writes there. Sizing an SGE
        # -l scratch request means adding this to the clustering section's
        # scratch_peak_mb, not taking the larger of the two.
        profile.record_size_bytes('scratch_environments_mb', environments_bytes)

        # ----- Run clus_and_deduplicate_vdgs.py -----
        with open(logfile, 'a') as f:
            f.write(f"\n----- Starting clus_and_deduplicate_vdgs.py -----\n")

        clus_script = os.path.join(HERE, "clus_and_deduplicate_vdgs.py")
        subset_sizes_arg = ' '.join(str(size) for size in subset_sizes)
        deduplicate_cmd = (
            f'python {shlex.quote(clus_script)} '
            f'-c {shlex.quote(cg)} -v {shlex.quote(out_dir)} '
            f'--cg-smarts {shlex.quote(smarts_raw)} '
            f'-P {shlex.quote(pdb_dir)} -E {shlex.quote(environments_out_root)} '
            f'--cg-match-dict-pkl {shlex.quote(match_pkl)} '
            f'-l {shlex.quote(logfile)} '
            + (f'-m {max_num_to_clus} ' if max_num_to_clus is not None else '')
            + f'--num-procs {num_procs} '
            f'--subset-sizes {subset_sizes_arg}')
        clustering_profile_path = os.path.join(out_dir, f'{cg}_clustering_profile.json')
        if profile.enabled:
            deduplicate_cmd += f' --profile-json {shlex.quote(clustering_profile_path)}'
        with profile.phase('cluster_and_deduplicate'):
            # Not check=True: EXIT_NO_VDGS is a warning, not a crash. The job
            # still exits 0 (nothing malfunctioned on this node), but we return
            # before _finish so the log never gets 'Job completed.' and the
            # post-build sweep lists the fragment as incomplete.
            clus_result = subprocess.run(deduplicate_cmd, shell=True)
        if clus_result.returncode == EXIT_NO_VDGS:
            # Clustering still wrote its profile sidecar before exiting. _finish is
            # skipped on this path, so nothing would ever fold or remove it; drop it
            # rather than leave a half-pair (clustering profile, no compute profile)
            # in the fragment directory.
            if os.path.exists(clustering_profile_path):
                os.remove(clustering_profile_path)
            # Printed to stderr, not written to the fragment log: the log is what
            # check_vdg_job_status substring-searches for the completion marker, and
            # nothing that lands in it may contain that phrase.
            print(f'[WARNING] {cg}: streamed no vdG records; see {logfile}. '
                  f'Leaving this fragment marked incomplete.', file=sys.stderr)
            return
        if clus_result.returncode != 0:
            raise subprocess.CalledProcessError(
                clus_result.returncode, deduplicate_cmd)
        if profile.enabled:
            # Clustering runs in its own process, so its numbers arrive through a
            # file; fold them in and leave one sidecar rather than two.
            profile.section('clustering', ComputeProfile.load(clustering_profile_path))
            if os.path.exists(clustering_profile_path):
                os.remove(clustering_profile_path)

        # Clean up nr_vdgs tree (remove empty dirs)
        delete_empty_dirs(os.path.join(out_dir, 'nr_vdgs'))
    finally:
        # Removed on every exit path (success, early return, or exception),
        # not just the success path -- this is a $TMPDIR/scratch tree, not out_dir.
        shutil.rmtree(environments_out_root, ignore_errors=True)

    _finish(profile, logfile, main_script_start, out_dir, cg)


def _finish(profile, logfile, main_script_start, out_dir, cg):
    """Write the completion marker and the profile sidecar.

    Shared by the normal exit and the zero-match early return: a run that
    matched nothing is a legitimate outcome, and without these it looks exactly
    like a job that died partway.
    """
    main_script_elapsed = time.time() - main_script_start
    hours, minutes, seconds = convert_time_elapsed(main_script_elapsed)

    with open(logfile, 'a') as _log:
        _log.write(f'{"="*79}\nJob completed.\nTotal job time: {hours} h, {minutes} mins, and {seconds} secs.\n')

    # One pointer line, not the numbers: the point of the sidecar is that a
    # production log stays readable.
    profile.set('total_wall_s', main_script_elapsed)
    profile.record_peak_rss()
    if profile.enabled:
        clustering_counters = profile.sections.get('clustering', {}).get('counters', {})
        profile.set('scratch_peak_total_mb',
                    round(profile.counters.get('scratch_environments_mb', 0.0)
                          + clustering_counters.get('scratch_peak_mb', 0.0), 1))
    profile_path = profile.write(os.path.join(out_dir, f'{cg}_compute_profile.json'))
    if profile_path:
        with open(logfile, 'a') as _log:
            _log.write(f'Compute profile: {profile_path}\n')


def delete_empty_dirs(_dir, keep_root=True):
    """Prune empty subdirectories. `_dir` itself is kept by default: a fragment
    that yields no vdGs must still leave nr_vdgs/ behind, or downstream readers
    (materialize_vdg_pdbs.py) fail with a bare 'no such directory' instead of
    seeing an empty library."""
    root_abs = os.path.abspath(_dir)
    for root, dirs, files in os.walk(_dir, topdown=False):
        if keep_root and os.path.abspath(root) == root_abs:
            continue
        if not dirs and not files:
            try:
                os.rmdir(root)
            except OSError:
                pass

def run_gen_environments(job_index, num_procs, environments_cmd):
    environments_cmd = f'{environments_cmd} -j {job_index} -n {num_procs}'
    subprocess.run(environments_cmd, shell=True, check=True)

def write_out_commandline_params(logfile, smarts, cg, pdb_dir, out_dir,
                                 num_automorphisms, max_num_to_clus, num_procs,
                                 subset_sizes):
    with open(logfile, 'w') as _log:
        _log.write(f'SMARTS: {smarts}\nCG: {cg}\n'
                   f'Number of exact CG automorphisms: {num_automorphisms}\n'
                   f'Max num vdgs to cluster: {max_num_to_clus}\nParent PDB dir: {pdb_dir}\n'
                   f'Contact-area threshold: {sasa.MIN_CONTACT_AREA} A^2 on '
                   f'buried+shared ({sasa.N_SPHERE_POINTS} SASA points)\n'
                   f'Output dir: {out_dir}\nNumber of processes: {num_procs}\n'
                   f'Subset sizes: {subset_sizes}\n')

if __name__ == '__main__':
    main()
