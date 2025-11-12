import os
import time
import subprocess
import argparse 
import shutil
import multiprocessing

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--smarts', type=str, required=True, 
                        help="SMARTS pattern.")
    parser.add_argument('-c', '--cg', type=str, 
                        help="The common name for the chemical group. Defaults "
                        "to the SMARTS pattern.")
    parser.add_argument('-p', "--pdb-dir", type=str, required=True,
                        help="Path to your PDB parent database with the same "
                        "directory structure as a mirror of the RCSB PDB; see "
                        "docs/database_generation_guide.md")
    parser.add_argument('-b', "--probe-dir", 
                        help="Path to your directory containing gzipped Probe "
                        "files; see docs/database_generation_guide.md")
    parser.add_argument('-o', "--out-dir", type=str, required=True,
                        help="Output directory path.")
    parser.add_argument('-k', "--keep-clustered-pdbs", action='store_true', 
                        help="Keep the `clusters/flankseq_and_bb/` dir. to keep track of "
                        "which clusters each nr_pdb came from.")
    parser.add_argument('-t', "--trial-run", type=int, 
                        help="Number of PDBs to process in a trial run used to "
                        "determine if smarts_to_cgs.py can run to completion "
                        "without errors.")
    parser.add_argument('--symmetry-classes', nargs='+', type=str,
                        help='Integers representing the symmetry classes of the CG '
                        'atoms on which clustering is to be performed.')
    parser.add_argument('-w', "--align-cg-weight", type=float, default=0.99, 
                        help="Fraction of weights to assign to CG atoms (collectively) "
                        "when superposing output vdGs. Not weights for clustering. "
                        "Example: 0.5 means 1/2 of weight is assigned to CG atoms and "
                        "the remaining 1/2 goes to the vdM backbone atoms. Defaults to 0.99.")
    parser.add_argument('-x', "--print-flankbb", action='store_true', 
                        help="Include flanking bb residues when writing out PDB.")
    parser.add_argument('-m', "--max-num-vdgs-to-clus", default=2500, type=int, 
                        help="Maximum number of PDBs w/ vdgs of the same AA compositions "
                        "to cluster.")
    parser.add_argument('--num-procs', type=int, default=10,
                        help="Number of processes to use for multiprocessing.")
    return parser.parse_args()


def main():
    args = parse_args()
    smarts = f'"{args.smarts}"'
    cg = args.cg 
    pdb_dir = args.pdb_dir
    probe_dir = args.probe_dir
    out_dir = args.out_dir
    keep_clustered_pdbs = args.keep_clustered_pdbs
    trial_run = args.trial_run
    align_cg_weight = args.align_cg_weight
    print_flankbb = args.print_flankbb
    symm_classes = args.symmetry_classes
    if symm_classes is not None:
        symm_classes = ' '.join(symm_classes)
    max_num_to_clus = args.max_num_vdgs_to_clus
    num_procs = args.num_procs

    if not cg: 
        cg = smarts

    # Set up outdir
    cg = cg.rstrip('"').lstrip('"')
    out_dir = os.path.join(out_dir, cg)
    if trial_run: 
        out_dir = out_dir.rstrip('/') + '_trial' 
    out_dir = set_up_outdir(out_dir)

    # Set logfile. If logfile exists, set a new logfile name.
    logdir = os.path.join(out_dir, 'logs')
    logfile = os.path.join(logdir, f'{cg}_log')
    if os.path.exists(logfile):
        logfile = logfile + '_' + str(time.time())
    write_out_commandline_params(logfile, smarts, cg, pdb_dir, probe_dir, out_dir, 
                                 symm_classes, logdir, max_num_to_clus, align_cg_weight, 
                                 num_procs)

    # Run smarts_to_cg.py
    if trial_run:
        smarts_to_cg_cmd = f'python external/vdG-miner/vdg_miner/programs/smarts_to_cgs.py -s {smarts} -c "{cg}" -p {pdb_dir} -o "{out_dir}" -l "{logfile}" -t {trial_run}'
    else:
        smarts_to_cg_cmd = f'python external/vdG-miner/vdg_miner/programs/smarts_to_cgs.py -s {smarts} -c "{cg}" -p {pdb_dir} -o "{out_dir}" -l "{logfile}"'
    
    subprocess.run(smarts_to_cg_cmd, shell=True, check=True)

    # Run generate_fingerprints.py
    match_pkl = os.path.join(out_dir, f'{cg}_matches.pkl') # output from smarts_to_cg.py
    fingerprints_cmd = f'python external/vdG-miner/vdg_miner/programs/generate_fingerprints.py -c "{cg}" -l "{logfile}" -m "{match_pkl}" -p {pdb_dir} -b {probe_dir} -o "{out_dir}"'

    # Pass tuple of arguments to starmap
    args_for_fingerprint = [(i, num_procs, fingerprints_cmd) for i in range(num_procs)]

    with multiprocessing.Pool(processes=num_procs, maxtasksperchild=1) as pool:
        pool.starmap(run_gen_fingerprints, args_for_fingerprint)

    # Run fingerprints_to_pdbs.py
    fingerprints = os.path.join(out_dir, 'fingerprints') # output from generate_fingerprints.py
    to_pdbs_cmd_template = (
        f'python external/vdG-miner/vdg_miner/programs/fingerprints_to_pdbs.py '
        f'-c "{cg}" -m "{match_pkl}" -l "{logfile}" -f "{fingerprints}" '
        f'-p {pdb_dir} -o "{out_dir}" -s -e')


    args_for_to_pdbs = [(i, num_procs, to_pdbs_cmd_template) for i in range(num_procs)]
    with multiprocessing.Pool(processes=num_procs, maxtasksperchild=1) as pool:
        pool.starmap(run_fingerprints_to_pdbs, args_for_to_pdbs)

    written_vdg_pdbs = os.listdir(os.path.join(out_dir, 'vdg_pdbs'))
    final_num_vdg_pdbs = len([f for f in written_vdg_pdbs  if f.endswith('.pdb.gz') 
        and os.path.isfile(os.path.join(out_dir, 'vdg_pdbs', f))])
    with open(logfile, 'a') as f:
        f.write(f"\nFinal total: {final_num_vdg_pdbs} vdg pdb files written out.\n")

    # Run clus_and_deduplicate_vdgs.py
    if symm_classes is not None:
        deduplicate_template = f'python ligand_vdgs/generate_vdgs/clus_and_deduplicate_vdgs.py -c "{cg}" -v "{out_dir}" -s {symm_classes} -l "{logfile}" -w {align_cg_weight} -m {max_num_to_clus} --num-procs {num_procs}'
    else:
        deduplicate_template = f'python ligand_vdgs/generate_vdgs/clus_and_deduplicate_vdgs.py -c "{cg}" -v "{out_dir}" -l "{logfile}" -w {align_cg_weight} -m {max_num_to_clus} --num-procs {num_procs}'

    # Add optional flags, if requested 
    if keep_clustered_pdbs:
        deduplicate_template += ' -k'
    if print_flankbb:
        deduplicate_template += ' -x'

    # Create a Pool of workers and run the deduplication commands concurrently
    num_vdms_list = [1, 2]
    with multiprocessing.Pool(processes=len(num_vdms_list), maxtasksperchild=1) as pool:
        pool.starmap(run_deduplicate, [(num_vdms, deduplicate_template) for num_vdms 
                                       in num_vdms_list])

    # Clean up the final state of clusters dir. Each subset's tempdir, flankseq, and 
    # flankbb dirs were cleaned up along the way, but the highest level of these dirs 
    # need to be deleted too.
    clean_up_clusdirs(out_dir, 'temp', logfile) 
    clean_up_clusdirs(out_dir, 'cgvdmbb', logfile) 
    delete_empty_dirs(os.path.join(out_dir, 'nr_vdgs'))
    if not keep_clustered_pdbs:
        clean_up_clusdirs(out_dir, 'flankseq_and_bb', logfile) 
        delete_empty_dirs(os.path.join(out_dir, 'clusters'))
    
    with open(logfile, 'a') as _log:
        _log.write(f'='*79 + '\n')
        _log.write(f'Job completed.\n')

def delete_empty_dirs(_dir):
    for root, dirs, files in os.walk(_dir):
        if not dirs and not files:
            os.system(f'rmdir "{root}"')

def clean_up_clusdirs(out_dir, clus_level, logfile):
    direc = os.path.join(out_dir, 'clusters', clus_level)
    with open(logfile, 'a') as _log:
        if not os.path.exists(direc):
            _log.write(f"\t{direc} does not exist.\n")
        else:
            # Check if there are files. 
            # They should have been deleted in the last set of clus_and_deduplicate_vdgs.py
            for root, dirs, files in os.walk(direc):
                if files:
                    _log.write(f"\t{direc} contains these files and will be deleted: \n")
                    for file in files:
                        _log.write(f"\t\t{file}\n")
            shutil.rmtree(direc)

def set_up_outdir(out_dir):
    # Set up output directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else:
        # Check to see if output directory is empty
        if os.listdir(out_dir):
            dir_exists_msg = (
                  f'\tThe output directory {out_dir} is not empty. Please remove its '
                  'contents or specify a different output directory path to avoid '
                  'accidental overwriting.\n')
            print('\n', dir_exists_msg)
            raise ValueError(dir_exists_msg)
    return out_dir

def run_deduplicate(num_vdms, deduplicate_template):
    deduplicate_cmd = f'{deduplicate_template} -n {num_vdms}'
    subprocess.run(deduplicate_cmd, shell=True, check=True)

def run_gen_fingerprints(job_index, num_procs, fingerprints_cmd):
    fingerprints_cmd = f'{fingerprints_cmd} -j {job_index} -n {num_procs}'
    subprocess.run(fingerprints_cmd, shell=True, check=True)

def run_fingerprints_to_pdbs(job_index, num_procs, cmd_template):
    cmd = f'{cmd_template} -j {job_index} -n {num_procs}'
    subprocess.run(cmd, shell=True, check=True)

def write_out_commandline_params(logfile, smarts, cg, pdb_dir, probe_dir, out_dir, 
                                 symm_classes, logdir, max_num_to_clus, align_cg_weight, 
                                 num_procs):
    if not os.path.exists(logdir):
        os.mkdir(logdir)
    #print(f'\nLogdir: {logdir}\n')
    with open(logfile, 'w') as _log:
        _log.write(f'SMARTS: {smarts} \n')
        _log.write(f'CG: {cg} \n')
        _log.write(f'Symmetry classes: {symm_classes} \n')
        _log.write(f'Weight to align CGs: {align_cg_weight} \n')
        _log.write(f'Max num vdgs to cluster: {max_num_to_clus} \n')
        _log.write(f'Parent PDB dir: {pdb_dir} \n')
        _log.write(f'Probe dir: {probe_dir} \n')
        _log.write(f'Output dir: {out_dir} \n')
        _log.write(f'Number of processes: {num_procs} \n')

if __name__ == '__main__':
    main()
