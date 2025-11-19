import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'ligand_vdgs', 'functions'))
import utils

# -------- Example input -----------
smiles = 'C=NC=NC'

# ----------------------------------
symm_classes = utils.define_symmetry(smiles)
print(symm_classes)