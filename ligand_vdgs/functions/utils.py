import os
import shutil
from rdkit import Chem

def handle_existing_files(out_dir):
   os.makedirs(out_dir, exist_ok=True)
   if len(os.listdir(out_dir)) > 0:
      raise ValueError(f'The output dir {out_dir} is not empty. Remove files or define a new '
            'output dir name to prevent accidental overwriting.')

def set_up_outdir(outdir, overwrite=False):
    '''Create outdir if it does not exist, and require overwrite_existing if it does, 
    so that there are no stale files.'''

    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            raise ValueError(f'The filename you designated as the output directory, {outdir}, '
                             'already exists and is not a directory.')
        if overwrite: 
            print(f'\nWarning: overwriting existing output directory {outdir} because '
                  'overwrite_existing was set to True.\n')
            shutil.rmtree(outdir) # removes the dir and everything in it
            os.makedirs(outdir)
        else: # check if there are already files. ignore if there are empty dirs.
            if len(os.listdir(outdir)) > 0:
                raise ValueError(f'The output directory {outdir} is not empty. Remove files or '
                                'set overwrite_existing to True to prevent accidental overwriting.')
    
    else: # dir doesn't exist, so create it
        parent_dir = os.path.dirname(outdir)
        if parent_dir: # skip cases where parent_dir is "" for the current working dir.
            if not os.path.exists(parent_dir):
                os.makedirs(parent_dir) 
            os.makedirs(outdir)

def valid_database_subdir_format(input_dir):
    ''' Ensure that the pdb database dir has subdirs formatted similarly to the RCSB
    mirror format (see docs/database_generation.txt). '''
    
    def error_message():
        print('Database structure must be similar to the RCSB mirror format; '
              'see docs/database_generation.txt')
    
    # Input path must be a dir
    if not os.path.exists(input_dir):
        error_message()
        return False
    # Input dir must not be empty
    if not len(os.listdir(input_dir)) > 0: 
        error_message()
        return False
    # Ensure that there is at least one correctly formatted subdir
    for p in os.listdir(input_dir):
        subdir_path = os.path.join(input_dir, p)
        if os.path.isdir(subdir_path):
            # Subdirs should be the inner 2 characters of a 4-char pdb file name
            if len(p) == 2: 
                for pdb in os.listdir(subdir_path):
                    if pdb[1:3] == p:
                        return True
    # No valid subdir was found
    error_message()
    return False

vdw_radii = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.80, 'S': 1.80, 'Cl': 1.75, 
            'Br': 1.85,'I': 1.98,}

def min_clash_dist(elem1, elem2, hbond_tol=0.6, general_tol=0.4):
    # Params: tolerance (Å) for decreasing the minimum allowable distance for atoms 
    # that may participate in hbonds (hbond_tol) or non-hbond interactions (general_tol).
    # Returns the minimum distance at which atoms are considered to clash.

    can_hbond = {frozenset(['N', 'O']), frozenset(['O', 'O']), frozenset(['N', 'N'])}

    if elem1 not in vdw_radii:
        print('Undefined vdw radius for elemen:', elem1)
        vdw1 = max(vdw_radii.values())
    else:
        vdw1 = vdw_radii[elem1]

    if elem2 not in vdw_radii:
        print('Undefined vdw radius for element:', elem2)
        vdw2 = max(vdw_radii.values())
    else:
        vdw2 = vdw_radii[elem2]

    min_dist = vdw1 + vdw2

    # Apply tolerance
    if frozenset([elem1, elem2]) in can_hbond:
        min_dist -= hbond_tol 
    else:
        min_dist -= general_tol 

    return min_dist

def get_atom_coords(prody_obj, atom_name):
    atom = prody_obj.select(f'name {atom_name}')
    coords = atom.getCoords()[0]
    return coords

def smiles_equiv(existingfrag, sub_smiles):
    mol1 = Chem.MolFromSmarts(existingfrag)
    mol2 = Chem.MolFromSmarts(sub_smiles)
    is_equivalent = (mol1.HasSubstructMatch(mol2) and 
                 mol2.HasSubstructMatch(mol1) and 
                 mol1.GetNumAtoms() == mol2.GetNumAtoms())
    return is_equivalent