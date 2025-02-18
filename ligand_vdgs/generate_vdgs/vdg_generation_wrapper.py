import os
import sys
import time
import subprocess
import argparse 
import shutil

def parse_args():
    parser = argparse.ArgumentParser(
        description="Determine CGs matching a SMARTS pattern."
    )
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
    parser.add_argument('-t', "--trial-run", type=int, 
                        help="Number of PDBs to process in a trial run used to "
                        "determine if smarts_to_cgs.py can run to completion "
                        "without errors.")
    parser.add_argument('--symmetry-classes', nargs='+', type=str,
                        help='Integers representing the symmetry classes of the CG '
                        'atoms on which clustering is to be performed.')
    return parser.parse_args()


def main():
    args = parse_args()
    smarts = f'"{args.smarts}"'
    cg = args.cg 
    pdb_dir = args.pdb_dir
    probe_dir = args.probe_dir
    out_dir = args.out_dir
    trial_run = args.trial_run
    symm_classes = args.symmetry_classes
    if symm_classes is not None:
        symm_classes = ' '.join(symm_classes)

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
                                 symm_classes, logdir)

    # Run smarts_to_cg.py
    if trial_run:
        smarts_to_cg_cmd = f'python ligand_vdgs/programs/vdG-miner/vdg_miner/programs/smarts_to_cgs.py -s {smarts} -c "{cg}" -p {pdb_dir} -o "{out_dir}" -l "{logfile}" -t {trial_run}'
    else:
        smarts_to_cg_cmd = f'python ligand_vdgs/programs/vdG-miner/vdg_miner/programs/smarts_to_cgs.py -s {smarts} -c "{cg}" -p {pdb_dir} -o "{out_dir}" -l "{logfile}"'
    
    subprocess.run(smarts_to_cg_cmd, shell=True, check=True)

    # Run generate_fingerprints.py
    match_pkl = os.path.join(out_dir, f'{cg}_matches.pkl') # output from smarts_to_cg.py
    fingerprints_cmd = f'python ligand_vdgs/programs/vdG-miner/vdg_miner/programs/generate_fingerprints.py -c "{cg}" -l "{logfile}" -m "{match_pkl}" -p {pdb_dir} -b {probe_dir} -o "{out_dir}"'

    subprocess.run(fingerprints_cmd, shell=True, check=True)

    # Run fingerprints_to_pdbs.py
    fingerprints = os.path.join(out_dir, 'fingerprints') # output from generate_fingerprints.py
    to_pdbs_cmd = f'python ligand_vdgs/programs/vdG-miner/vdg_miner/programs/fingerprints_to_pdbs.py -c "{cg}" -m "{match_pkl}" -l "{logfile}" -f "{fingerprints}" -p {pdb_dir} -o "{out_dir}" -s -e'

    subprocess.run(to_pdbs_cmd, shell=True, check=True)

    # Run clus_and_deduplicate_vdgs.py
    if symm_classes is not None:
        deduplicate_template = f'python ligand_vdgs/generate_vdgs/clus_and_deduplicate_vdgs.py -c "{cg}" -v "{out_dir}" -s {symm_classes} -l "{logfile}" -n '
    else:
        deduplicate_template = f'python ligand_vdgs/generate_vdgs/clus_and_deduplicate_vdgs.py -c "{cg}" -v "{out_dir}" -l "{logfile}" -n '

    for num_vdms in [1,2,3,4]:
        deduplicate_cmd = deduplicate_template + str(num_vdms)
        subprocess.run(deduplicate_cmd, shell=True, check=True)
        
    # Clean up the final state of clusters dir. Each subset's tempdir, flankseq, and 
    # flankbb dirs were cleaned up along the way, but the highest level of these dirs 
    # need to be deleted too.
    clean_up_dirs(out_dir, 'temp', logfile) 
    clean_up_dirs(out_dir, 'flankseq', logfile) 
    clean_up_dirs(out_dir, 'flankbb', logfile) 
    
    with open(logfile, 'a') as _log:
        _log.write(f'='*79 + '\n')
        _log.write(f'Job completed.\n')

def clean_up_dirs(out_dir, clus_level, logfile):
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
                    return
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

def write_out_commandline_params(logfile, smarts, cg, pdb_dir, probe_dir, out_dir, 
                                 symm_classes, logdir):
    if not os.path.exists(logdir):
        os.mkdir(logdir)
    print(f'\nLogdir: {logdir}\n')
    with open(logfile, 'w') as _log:
        _log.write(f'SMARTS: {smarts} \n')
        _log.write(f'CG: {cg} \n')
        _log.write(f'Symmetry classes: {symm_classes} \n')
        _log.write(f'Parent PDB dir: {pdb_dir} \n')
        _log.write(f'Probe dir: {probe_dir} \n')
        _log.write(f'Output dir: {out_dir} \n')

if __name__ == '__main__':
    main()
