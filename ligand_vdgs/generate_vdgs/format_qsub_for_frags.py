import os
import sys
from rdkit import Chem
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import utils

outfile = 'run_frags.sh'
submit_scripts_dir = 'ligand_vdgs/generate_vdgs/frag_submit_scripts/'

# ------------------------------------------------------------------------------------------

with open(outfile, 'w') as out:
    for f in os.listdir(submit_scripts_dir):
        smiles = f.rstrip('.sh')
        f_path = os.path.join(submit_scripts_dir, f)
        out.write(f'qsub "{f_path}"\n')

