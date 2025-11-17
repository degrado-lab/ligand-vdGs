import os
import sys
from rdkit import Chem
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))
import utils

outfile_no_symm = 'run_frags_no_symm.sh'
outfile_with_symm = 'run_frags_symm.sh'
submit_scripts_dir = 'ligand_vdgs/generate_vdgs/frag_submit_scripts/'

# ------------------------------------------------------------------------------------------

with open(outfile_no_symm, 'w') as out_no_symm, open(outfile_with_symm, 'w') as out_with_symm:
    for f in os.listdir(submit_scripts_dir):
        smiles = f.rstrip('.sh')
        # determine if smiles is symmetric
        mol = Chem.MolFromSmarts(smiles)
        symm_classes = utils.identify_mol_symmetry(mol)
        if len(symm_classes) == mol.GetNumAtoms():
            # no symmetry
            out = out_no_symm
        else:
            out = out_with_symm
        f_path = os.path.join(submit_scripts_dir, f)
        out.write(f'qsub "{f_path}"\n')

