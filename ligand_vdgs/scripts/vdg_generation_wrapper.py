import os
import sys
import time
import subprocess
import argparse 

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
    return parser.parse_args()


def main():
    args = parse_args()
    smarts = f'"{args.smarts}"'
    cg = args.cg 
    pdb_dir = args.pdb_dir
    probe_dir = args.probe_dir
    out_dir = args.out_dir
    trial_run = args.trial_run

    if not cg: 
        cg = smarts

    # Set logfile. If logfile exists, set a new logfile name.
    logfile = os.path.join(out_dir, 'logs', f'{cg}_log')
    if os.path.exists(logfile):
        logfile = logfile + '_' + time.time()

    # Run smarts_to_cg.py
    if trial_run: 
        smarts_to_cg_cmd = f'python ligand_vdgs/programs/vdG-miner/vdg_miner/programs/smarts_to_cgs.py -s {smarts} -c {cg} -p {pdb_dir} -o {out_dir} -l {logfile} -t {trial_run}'
    else:
        smarts_to_cg_cmd = f'python ligand_vdgs/programs/vdG-miner/vdg_miner/programs/smarts_to_cgs.py -s {smarts} -c {cg} -p {pdb_dir} -o {out_dir} -l {logfile}'
    
    subprocess.run(smarts_to_cg_cmd, shell=True, check=True)

    # Run generate_fingerprints.py
    match_pkl = os.path.join(out_dir, cg, f'{cg}_matches.pkl') # output from smarts_to_cg.py
    fingerprints_cmd = f'python ligand_vdgs/programs/vdG-miner/vdg_miner/programs/generate_fingerprints.py -c {cg} -l {logfile} -m {match_pkl} -p {pdb_dir} -b {probe_dir} -o {out_dir}'

    subprocess.run(fingerprints_cmd, shell=True, check=True)

    # Run fingerprints_to_pdbs.py
    fingerprints = os.path.join(out_dir, cg, 'fingerprints') # output from generate_fingerprints.py
    to_pdbs_cmd = f'python ligand_vdgs/programs/vdG-miner/vdg_miner/programs/fingerprints_to_pdbs.py -c {cg} -m {match_pkl} -l {logfile} -f {fingerprints} -p {pdb_dir} -o {out_dir} -s -e'

    subprocess.run(to_pdbs_cmd, shell=True, check=True)

if __name__ == '__main__':
    main()