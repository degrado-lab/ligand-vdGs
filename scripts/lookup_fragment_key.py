'''Trace a fragment key to the frags-dict keys and library directory it covers.

Library directory names are representatives, not necessarily dict keys: a
charged variant lives under its neutral twin, and a group with no twin lives
under a *promoted* charge-stripped SMARTS that appears in no CCD drawing
(extract_fragment_smiles.group_protonation_variants). Looking such a key up in
the dict by string finds nothing. This resolves it the way the pipeline does --
exact, then equivalent atom order, then charge-loose -- and prints what it covers.

    python scripts/lookup_fragment_key.py 'c[N;!R](=[O;!R])[O;!R]'
    python scripts/lookup_fragment_key.py KEY --vdg-lib-dir ../frag_lib
'''
import argparse
import os
import pickle as pkl

from ligand_vdgs.functions import utils
from ligand_vdgs.functions.vdg_npz_utils import load_fragment_aliases
from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (load_frags_dict,
                                                                  resolve_fragment_key)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('key', help='Annotated fragment SMARTS (ring annotations included).')
    p.add_argument('--frags-dict', default='resources/database_frags_dict.pkl')
    p.add_argument('--vdg-lib-dir', help='Also report alias rows and the library '
                   'directory for the key.')
    p.add_argument('--ligands', action='store_true',
                   help='List the CCD ligand names under each covering dict key.')
    return p.parse_args()


def main():
    args = parse_args()
    # load_frags_dict, not a raw pkl.load: the dict is a wrapper carrying
    # key_schema, and a raw load would answer questions about a
    # vocabulary the running code cannot match.
    frags_dict, _support, _meta = load_frags_dict(args.frags_dict)
    ligands = {}
    for smiles_dict in frags_dict.values():
        for smiles, names in smiles_dict.items():
            ligands.setdefault(smiles, set()).update(names)

    found = resolve_fragment_key(args.key, ligands)
    if not found:
        print(f'{args.key}: matches no frags-dict key, not even charge-different or '
              'reordered. Check the spelling, ring annotations included.')
    else:
        print(f'{args.key}: covers {len(found)} frags-dict key(s), '
              f'{len(set().union(*(ligands[k] for k, _ in found)))} distinct ligand(s)')
        for cand, how in found:
            print(f'  {how:<14s} {cand}  ({len(ligands[cand])} ligands)')
            if args.ligands:
                print('      ' + ' '.join(sorted(ligands[cand])))

    if args.vdg_lib_dir:
        aliases = load_fragment_aliases(args.vdg_lib_dir)
        rep = aliases.get(args.key, args.key)
        if rep != args.key:
            print(f'library: {args.key} is an alias of {rep}')
        covered = sorted(a for a, r in aliases.items() if r == rep)
        if covered:
            print(f'library: {rep} covers alias(es) {covered}')
        frag_dir = os.path.join(args.vdg_lib_dir, utils.smiles_to_filename(rep))
        state = 'exists' if os.path.isdir(frag_dir) else 'does not exist'
        print(f'library: directory {frag_dir} {state}')


if __name__ == '__main__':
    main()
