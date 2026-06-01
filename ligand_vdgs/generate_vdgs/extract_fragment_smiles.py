'''
Extract qualifying fragment SMILES from database_frags_dict.pkl and write them
to a text file (one SMILES per line).

Fragments qualify if they appear in at least --counts-threshold database ligands
and contain no more than --max-size heavy atoms. These are the same thresholds
used by make_sge_scripts_for_frags.py.

The output file can be used to submit one vdg_generation_wrapper.py job per
fragment on any cluster scheduler (SLURM, PBS, etc.).

Usage:
    python ligand_vdgs/generate_vdgs/extract_fragment_smiles.py
    python ligand_vdgs/generate_vdgs/extract_fragment_smiles.py \
        --frags-dict <path/to/database_frags_dict.pkl> \
        --output     <path/to/fragment_smiles.txt>

Output:
    A text file with one fragment SMILES per line.
'''

import os
import sys
import argparse
import pickle as pkl
from rdkit import Chem
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'functions'))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract qualifying fragment SMILES from database_frags_dict.pkl.")
    parser.add_argument('--frags-dict', default='resources/database_frags_dict.pkl',
                        help="Path to database_frags_dict.pkl. "
                             "Default: resources/database_frags_dict.pkl.")
    parser.add_argument('--output', default='resources/fragment_smiles.txt',
                        help="Path to write the output SMILES file. "
                             "Default: resources/fragment_smiles.txt.")
    parser.add_argument('--counts-threshold', default=80, type=int,
                        help="Min number of ligands containing a fragment to qualify. "
                             "Default: 80.")
    parser.add_argument('--max-size', default=5, type=int,
                        help="Max fragment heavy-atom count. Default: 5.")
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.isfile(args.frags_dict):
        raise FileNotFoundError(f"--frags-dict path does not exist: {args.frags_dict}")

    with open(args.frags_dict, 'rb') as f:
        frags_dict = pkl.load(f)

    smiles_to_run = []
    for _, smiles_dict in frags_dict.items():
        for smiles, lignames in smiles_dict.items():
            if len(set(lignames)) < args.counts_threshold:
                continue
            mol = Chem.MolFromSmarts(smiles)
            if mol is None or mol.GetNumHeavyAtoms() > args.max_size:
                continue
            smiles_to_run.append(smiles)

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.output, 'w') as f:
        for smi in smiles_to_run:
            f.write(smi + '\n')

    print(f'Wrote {len(smiles_to_run)} fragment SMILES to {args.output}.')


if __name__ == '__main__':
    main()
