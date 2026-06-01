'''
Generate per-fragment SGE job scripts for running the vdG generation pipeline.

For each fragment in database_frags_dict.pkl that meets the counts and size thresholds,
this script writes one SGE shell script that calls vdg_generation_wrapper.py for that
fragment. The scripts are written to --sge-out-dir and can be submitted with:
    qsub <script>.sh

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
        --log-dir   <path/to/sge_logs/>
'''

import os
import sys
import argparse
import pickle as pkl
from rdkit import Chem
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import utils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create SGE submission scripts for vdG generation.")
    parser.add_argument('--frags-dict', default='resources/database_frags_dict.pkl',
                        help="Path to database_frags_dict.pkl.")
    parser.add_argument('--template', default='resources/frag_sge_template.sh',
                        help="Path to SGE job script template.")
    parser.add_argument('--counts-threshold', default=80, type=int,
                        help="Min number of ligands containing a fragment to run vdG "
                             "generation. Default: 80.")
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
    parser.add_argument('--max-num-clus', default='5000',
                        help="Max number of vdGs to cluster per subset. Default: 5000.")
    parser.add_argument('--h-rt', default='300:00:00',
                        help="SGE wall-clock time limit. Default: 300:00:00.")
    parser.add_argument('--num-procs', default='20',
                        help="Number of parallel processes. Default: 20.")
    return parser.parse_args()


def main():
    args = parse_args()

    for path, flag in [(args.frags_dict, '--frags-dict'), (args.template, '--template')]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"{flag} path does not exist: {path}")

    replace = {'$LOG_DIR':      args.log_dir,
               '$PDB_DIR':      args.pdb_dir,
               '$PROBE_DIR':    args.probe_dir,
               '$OUTPUT_DIR':   args.vdg_lib_dir,
               '$MAX_NUM_CLUS': args.max_num_clus,
               '$RUN_TIME':     args.h_rt,
               '$NUM_PROCS':    args.num_procs}

    if not os.path.exists(args.sge_out_dir):
        os.makedirs(args.sge_out_dir)
    if os.listdir(args.sge_out_dir):
        raise FileExistsError(f"Output directory {args.sge_out_dir} already has "
                              "files. Terminating to prevent overwriting.")

    smiles_to_run = set()
    with open(args.frags_dict, 'rb') as f:
        frags_dict = pkl.load(f)
    for _, smiles_dict in frags_dict.items():
        for smiles, lignames in smiles_dict.items():
            # apply counts threshold
            if len(set(lignames)) < args.counts_threshold:
                continue
            # apply size threshold
            mol = Chem.MolFromSmarts(smiles)
            if mol is None:
                continue
            if mol.GetNumHeavyAtoms() > args.max_size:
                continue
            smiles_to_run.add(smiles)

    # iterate over smiles and create scripts
    for smiles in smiles_to_run:
        # determine symmetry
        symm_classes = utils.define_symmetry(smiles)
        # output the job script for this frag
        output_script(args.template, smiles, symm_classes, args.sge_out_dir, replace)

    print(f'Created scripts for {len(smiles_to_run)} fragments in {args.sge_out_dir}.')


def output_script(template, smiles, symm_classes, sge_out_dir, replace):
    # per-fragment copy of the base mapping. `local_replace` is b/c symm_classes
    # definitions would persist onto the next frag.
    local_replace = dict(replace)

    script_name = os.path.join(sge_out_dir, utils.smiles_to_filename(smiles) + '.sh')
    local_replace['$SMILES'] = f'"{smiles}"'
    local_replace['$CG'] = f'"{smiles}"'
    local_replace['$JOB_NAME'] = utils.smiles_to_job_name(smiles)  # SGE job names: # truncates directives

    # always define how to handle --num-procs
    if symm_classes:
        local_replace['--num-procs'] = f'--symmetry-classes {symm_classes} --num-procs'
    else:
        local_replace['--num-procs'] = '--num-procs'

    with open(template, 'r') as f:
        lines = f.readlines()

    with open(script_name, 'w') as f:
        start_copy = False
        for line in lines:
            # skip header; wait until #!/bin/bash
            if line.startswith('#!/bin/bash'):
                start_copy = True
            if start_copy:
                # replace placeholders based on the `replace` dict
                for key, value in local_replace.items():
                    line = line.replace(key, value)
                f.write(line)


if __name__ == '__main__':
    main()
