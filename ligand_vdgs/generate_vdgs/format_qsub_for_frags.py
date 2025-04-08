import os

outfile = 'run_frags.sh'
submit_scripts_dir = 'ligand_vdgs/generate_vdgs/lib_gen_frag_submit_scripts/'

# ------------------------------------------------------------------------------------------

with open(outfile, 'w') as out:
    for f in os.listdir(submit_scripts_dir):
        f_path = os.path.join(submit_scripts_dir, f)
        out.write(f'qsub "{f_path}"\n')

