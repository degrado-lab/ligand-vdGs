'''
The `ligand-vdGs` and `vdG-miner` packages require the input PDB database (the "parent database")
to match the directory structure
used in the RCSB PDB mirror. This script formats your PDB database to conform to that structure.

Each PDB file will be organized into subdirectories based on the inner two characters of its name. 
For example, a PDB file named 
`1ABC.pdb` will be placed into a subdirectory named `AB`, resulting in the file path `AB/1ABC.pdb`.

Please note that your PDB files should be named with a 4-character identifier (XXXX.pdb).

Usage:
    >> cd $YOUR_LIGAND-VDGS_DIR
    >> python scripts/format_parent_database.py \
        --input_dir <input_dir> \
        --output_dir <output_dir>
'''

import os
import argparse
import shutil
from ligand_vdgs.functions import parent_db

pdb_input_dir = 'consolidated_BioLiP2'
pdb_output_dir = 'consolidated_BioLiP2_split'

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', 
                        required=True, 
                        help='Path to the PDB directory that needs to be reformatted.')
    parser.add_argument('--output_dir', 
                        required=True, 
                        help='Path to output the reformatted database.')
    return parser.parse_args()

def main():
    args = parse_args()
    pdb_input_dir = args.input_dir
    pdb_output_dir = args.output_dir
    # Set up output directory
    if os.path.exists(pdb_output_dir):
        if len(os.listdir(pdb_output_dir)) > 0:
            raise ValueError('The output directory is not empty. Process terminated to avoid overwriting existing files.')
    else:
        os.makedirs(pdb_output_dir)

    pdbs = os.listdir(pdb_input_dir)

    # Validate all file names up front so a bad name aborts before any file is touched
    for pdb in pdbs:
        if not pdb.endswith(parent_db.STRUCTURE_EXT) or len(parent_db.stem_of(pdb)) < 3:
            raise ValueError('PDB file names are expected to be <stem>.pdb with a stem of '
                             'at least 3 characters. Invalid file: {}'.format(
                                 os.path.join(pdb_input_dir, pdb)))

    # Copy the pdbs (source database is left intact)
    for pdb in pdbs:
        print(pdb)
        # Determine its subdir, based on the inner 2 characters in the pdb name
        subdir_code = parent_db.shard(parent_db.stem_of(pdb))
        subdir_path = os.path.join(pdb_output_dir, subdir_code)
        if not os.path.exists(subdir_path):
            os.mkdir(subdir_path)
        source_path = os.path.join(pdb_input_dir, pdb)

        # Copy the file
        shutil.copy2(source_path, subdir_path)

if __name__ == "__main__":
    main()
