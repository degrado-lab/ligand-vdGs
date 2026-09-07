'''
Create SGE submission scripts to run vdg_hit_finder.py.

Reads `models_to_run_file` (set as a module-level constant below, not a CLI
argument), a csv file where each line is a comma-separated triplet of
(num_procs, query_dir [dir containing model pdb files], and a SMILES string
for the ligand in those pdb files). All models in query_dir must have the same ligand.

Note: provide the ligand SMILES instead of deriving it from the structure files because 
rdkit sometimes parses them incorrectly.

Hit-finder options supplied to this generator are added to every generated job, e.g.:

    python ligand_vdgs/score_poses/make_sge_scripts_for_hit_finder.py \
        --rmsd-threshold 1.5 --contact-cutoff 3.8 --no-dedup
'''

import argparse
import csv
import os
import shlex

template = 'resources/vdg_hit_finder_sge_template.sh' 
out_dir_for_sge_scripts = 'ligand_vdgs/score_poses/hit_finder_submit_scripts/'
vdg_lib_dir = '/wynton/home/degradolab/skt/docking/frag_lib'
models_to_run_file = 'resources/validation_set.csv'
log_dir = '/wynton/home/degradolab/skt/logs/hit_finder_logs'
vdg_hit_outdir = '/wynton/home/degradolab/skt/docking/vdg_hits' 

# ------------------------------------------------------------------------------------------

base_replace = {'$LOG_DIR': log_dir,
                '$VDG_LIB_DIR': shlex.quote(vdg_lib_dir),}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Create SGE submission scripts for vdg_hit_finder.py."
    )
    parser.add_argument(
        "--rmsd-threshold",
        type=float,
        help="Pass --rmsd-threshold to every hit-finder job.",
    )
    parser.add_argument(
        "--ref-pdb",
        help="Pass --ref-pdb to every hit-finder job.",
    )
    parser.add_argument(
        "--print-bsr-selection",
        action="store_true",
        help="Pass --print-bsr-selection to every hit-finder job.",
    )
    parser.add_argument(
        "--contact-cutoff",
        type=float,
        help="Pass --contact-cutoff to every hit-finder job.",
    )
    parser.add_argument(
        "--no-dedup",
        action="store_true",
        help="Pass --no-dedup to every hit-finder job.",
    )
    parser.add_argument(
        "--min-shared-atoms",
        type=int,
        help="Pass --min-shared-atoms to every hit-finder job.",
    )
    return parser.parse_args(argv)


def optional_hit_finder_args(args):
    """Return shell-quoted optional CLI arguments for vdg_hit_finder.py."""
    command_args = []
    for flag, value in (
        ("--rmsd-threshold", args.rmsd_threshold),
        ("--ref-pdb", args.ref_pdb),
        ("--contact-cutoff", args.contact_cutoff),
        ("--min-shared-atoms", args.min_shared_atoms),
    ):
        if value is not None:
            command_args.extend((flag, str(value)))
    if args.print_bsr_selection:
        command_args.append("--print-bsr-selection")
    if args.no_dedup:
        command_args.append("--no-dedup")
    return shlex.join(command_args)


def main(argv=None):
    args = parse_args(argv)
    optional_args = optional_hit_finder_args(args)
    if not os.path.exists(out_dir_for_sge_scripts):
        os.makedirs(out_dir_for_sge_scripts)
    if os.listdir(out_dir_for_sge_scripts):
        raise FileExistsError(f"Output directory {out_dir_for_sge_scripts} already has "
                              "files. Terminating to prevent overwriting.")
    # The template is shared by every generated job script, so read it once
    # rather than reopening it for each row in models_to_run_file.
    with open(template, 'r') as f_tmpl:
        template_lines = f_tmpl.readlines()

    # Read models_to_run_file
    num_structs = 0
    with open(models_to_run_file, 'r', newline='') as f_models:
        reader = csv.reader(f_models)
        for line_number, row in enumerate(reader, start=1):
            # Ignore empty lines, and normalize fields before unpacking.  This
            # keeps whitespace-only lines from causing an unpacking error and
            # prevents whitespace in paths/SMILES from leaking into scripts.
            row = [field.strip() for field in row]
            if not row or not any(row):
                continue
            if len(row) != 3:
                raise ValueError(
                    f'{models_to_run_file}:{line_number}: expected 3 CSV fields '
                    f'(num_procs, query_dir, smiles), got {len(row)}'
                )
            num_structs += 1
            num_procs, query_dir, smiles = row
            if not num_procs.isdigit() or int(num_procs) <= 0:
                raise ValueError(
                    f'{models_to_run_file}:{line_number}: num_procs must be a '
                    f'positive integer, got {num_procs!r}'
                )
            # Keep substitutions local to this row.  Mutating a module-level dict
            # here leaks values between rows and between repeated calls to main().
            replace = {
                **base_replace,
                '$QUERY_DIR': shlex.quote(query_dir),
                '$OUTDIR': shlex.quote(os.path.join(
                    vdg_hit_outdir,
                    os.path.basename(os.path.normpath(query_dir)),
                )),
                '$SMILES': shlex.quote(smiles),
                '$JOB_NAME': f'_{os.path.basename(os.path.normpath(query_dir))}',
                '$NUM_PROCS': shlex.quote(num_procs),
                '$OPTIONAL_ARGS': optional_args,
            }
            name_structs = os.path.basename(os.path.normpath(query_dir))
            print(name_structs)
            job_script = os.path.join(out_dir_for_sge_scripts, name_structs + '.sh')
            # Create SGE script from template
            with open(job_script, 'w') as f_out:
                start_copy = False
                for line in template_lines:
                    # skip header. wait until #!/bin/bash
                    if line.startswith('#!/bin/bash'):
                        start_copy = True
                    if start_copy:
                        # replace placeholders based on the `replace` dict
                        for key, value in replace.items():
                            line = line.replace(key, value)
                        f_out.write(line)
    print(f'RUNNING {num_structs} STRUCTS')


if __name__ == "__main__":
    main()
