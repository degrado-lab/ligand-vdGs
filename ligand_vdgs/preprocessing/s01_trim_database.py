'''
Stand-alone script.
Reduces size of the 50G consolidated_BioLiP2_split/ database by 
- filtering by b-factor
- crudely determining if a binding site is equivalent by collecting all the vdm
  resnums and resnames of a lig, and determining if that same set of vdm resnums,
  resnames, and lig resnums/resnames are repeated in that same pdb or a diff pdb.
  (there's some logic about taking a guess at whether the pdbs are within the same
  series or not - double check that and describe here.)

Usage:
cd $YOUR_LIGAND-VDGS_DIR
python ligand_vdgs/preprocessing/s01_trim_database.py

NOTE: large assemblies in BioLiP2 carry two-character chain IDs written across PDB
columns 21-22 (``ASPA4  66`` = chain A4). Every reader downstream sees only column
22, so chains ``A4`` and ``4`` merge residue by residue and vdG-miner drops the
structure whole. _chain_ids.remap_chain_ids rewrites them to single characters
before ProDy parses the file (parsing first would bake the merge in silently);
the mapping is recorded as REMARK lines in the output PDB.

NOTE: modified residues are handled here too: MSE/SEC become
MET/CYS, and other chain-bonded HETATM residues of CCD type `*PEPTIDE LINKING`
become ATOM records under their own resname, so they are neither mined as ligands
nor accepted as slots. Needs resources/ccd_polymer_types.tsv
(scripts/fetch_ccd_polymer_types.py). See _prep_filters.modified_residues_to_protein.
'''

import gzip
import io
import os
import traceback
import json
import numpy as np
import prody as pr

from ligand_vdgs.functions.interactions import add_pdb_to_nr_db_dict
from ligand_vdgs.functions.utils import set_up_outdir
from ligand_vdgs.functions import parent_db
from ligand_vdgs.preprocessing._prep_filters import (
    drop_prepwizard_hazard_residues, load_ccd_polymer_types, modified_residues_to_protein)
from ligand_vdgs.preprocessing._chain_ids import (ChainIdOverflow, remap_chain_ids,
                                                  remap_remarks)

origin_dir = '/home/sophia/DockDesign/databases/consolidated_BioLiP2_split'
target_dir = '/home/sophia/DockDesign/databases/consolidated_BioLiP2_trimmed'
lig_avg_bfactor_cutoff = 40 # default is 40
output_database_dict_name = '20240809_database_dict.json'
overwrite_pdbs = False # If set to true NEED TO ADD code for checking whether a pdbfile already exists (TODO)
skip_to_output_pdbs = True
# Sphere radius around each ligand, in A. Measured over 124 deposited entries
# (contact shell = 4.5 A): 10 A retains 93.2% of a contact vdM's +/-2 flanking
# residues, 83.5% of the ligand's second shell and 8.2% of its third; 16 A
# completes the second shell, and 20 A is the first radius that also reaches
# the third (99.5%). 20 A costs ~3.0x the atoms of 10 A. Third-shell context is
# what a protonation program needs to place a buried His/Asp correctly, so the
# radius is set for that rather than for the flank walk.
radius = 20

pdbs_already_done = 'log_processed_pdbs.txt'
previous_checkpoint_dict = 'checkpoint.pkl'


def parse_pdb_with_single_char_chains(pdbpath, ccd_types):
    '''ProDy AtomGroup for *pdbpath* after the textual fixes that must precede parsing:
    multi-character chain IDs remapped to single characters, MSE/SEC renamed to MET/CYS,
    HETATM amino acids and chain-bonded modified residues made ATOM records
    (_prep_filters.modified_residues_to_protein; *ccd_types* from
    load_ccd_polymer_types). Returns (atoms, {old: new} chain mapping, stats dict from
    modified_residues_to_protein). Raises ChainIdOverflow if the structure has more
    chains than single characters.'''
    opener = gzip.open if pdbpath.endswith('.gz') else open
    with opener(pdbpath, 'rt') as f:
        lines = f.readlines()
    lines, mapping = remap_chain_ids(lines)
    lines, stats = modified_residues_to_protein(lines, ccd_types)
    atoms = pr.parsePDBStream(io.StringIO(''.join(lines)))
    return atoms, mapping, stats


def write_pdb_with_remap_remarks(output_path, atoms, mapping):
    '''pr.writePDB, with the chain remap recorded as REMARK lines up front.'''
    buf = io.StringIO()
    pr.writePDBStream(buf, atoms)
    with open(output_path, 'w') as f:
        f.writelines(remap_remarks(mapping))
        f.write(buf.getvalue())


def main():

    checkpoint_name = output_database_dict_name + '.checkpoint'

    # Iterate through pdb files
    pdbfile_ix = 0

    # Initialize nested dict to keep track of pdbs, segs, chains, resnums, resnames, etc.
    # Nested database_ligands_dict format: lig resname:
    #                                           pdb: 
    #                                                lig res: 
    #                                                      list of interacting residues 
    #                                                      (segment, chain, rensum, resname)
    # See interactions.add_pdb_to_nr_db_dict() for more details.
    database_dict = {}

    # CCD types decide which chain-bonded HETATM residues are modified amino acids
    # (-> ATOM records) rather than covalent cofactors (stay ligands).
    ccd_types = load_ccd_polymer_types()
    unknown_resnames = set()

    #if restart: # the name restart could be confusng bc now it should be reload
    #    with open(pdbs_already_done) as inF:
    #        pdbs_done = [i.strip() for i in inF.readlines()]
    #    print("Starting json load from previous run", flush=True) # TODO: problem; there might be >1 prev runs
    #    database_dict = json.load(open(previous_checkpoint_dict, 'r'))
    #    print('Num of PDBs processed in a previous run: ', 
    #        len([i for i in pdbs_done if i.endswith('.pdb')]))

    if not skip_to_output_pdbs:
        for subdir in os.listdir(origin_dir):
            subdir_path = os.path.join(origin_dir, subdir)

            pdbfiles = os.listdir(subdir_path)

            for pdbfile in pdbfiles:
                pdbfile_ix += 1
                #if restart:
                #    if pdbfile in pdbs_done:
                #        continue
                print(pdbfile, flush=True)
                try:

                    pdbpath = os.path.join(subdir_path, pdbfile)
                    atoms, _, stats = parse_pdb_with_single_char_chains(pdbpath, ccd_types)
                    unknown_resnames.update(stats['unknown_resnames'])

                    # Add the ligand(s) and interacting residues to the ligand dict
                    # if nonredundant. 
                    # Redundancy here is a quick and dirty check (computationally
                    # cheap); it is refined later in the vdG creation process.
                    database_dict = add_pdb_to_nr_db_dict(
                        database_dict, pdbpath,
                        lig_bfactor_cutoff=lig_avg_bfactor_cutoff, atoms=atoms)
                except ChainIdOverflow as e:
                    print(f'[ERROR] Skipping {pdbfile}: {e}')
                except Exception as e:
                    print(f'[ERROR] PDB failed: {pdbfile}')
                    print(e)
                    traceback.print_exc()

                if pdbfile_ix % 1000 == 0:
                    print(pdbfile_ix)
                    print(f"Dumping on pdbfile_ix: {pdbfile_ix}", flush=True)
                    json.dump(database_dict, open(checkpoint_name, "w"))

    
        print("Starting final JSON dump", flush=True)
        json.dump(database_dict, open(output_database_dict_name, "w"))
        print("Done with final JSON dump", flush=True)

    # Now that the ligand instances have been deduplicated, print out only the binding sites
    # (*radius* A around the lig)
    
    AAs = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 
            'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
    
    database_dict = json.load(open(output_database_dict_name, 'r'))
    
    set_up_outdir(target_dir, overwrite_pdbs)
    # First, get all the lig and vdm residues in each pdb, so you only have to process each 
    # pdb once.
    bindingsite_dict = {} # key=pdb, values=list of ligs, where lig=[segment, chain, resnum]
    for ligresn, lig_interactions in database_dict.items():
        # SKIP OVER LIGAND IF IT'S ACTUALLY AN AA
        if len(ligresn) == 4 and ligresn[:3] in AAs:
            continue
        for lig_res in lig_interactions.keys():
            # Keys are ' '-joined 'pdbfile seg chain resnum resname' (see interactions.py);
            # seg and chain may be empty strings.
            pdbname, _seg, _ch, _resnum, _resname = lig_res.split(' ')
            bindingsite_dict.setdefault(pdbname, []).append((_seg, _ch, int(_resnum)))

    # Then, select *radius* around the lig residues, and write out the pdb.
    for pdb_name, list_ligs in bindingsite_dict.items():
        # First, determine if this pdb was already written out in a previous run (if this script is
        # resuming from a previous incomplete run).
        original_pdb_subdir = parent_db.shard(parent_db.stem_of(pdb_name))
        output_subdir = os.path.join(target_dir, original_pdb_subdir)
        output_path = os.path.join(output_subdir, pdb_name)
        if not overwrite_pdbs and os.path.exists(output_path):
            continue

        # Load pdb
        original_pdbpath = os.path.join(origin_dir, original_pdb_subdir, pdb_name)
        try:
            # Same remap as the dict-building pass, so the seg/chain/resnum keys in
            # database_dict refer to the remapped chains.
            parsed, chain_mapping, stats = parse_pdb_with_single_char_chains(
                original_pdbpath, ccd_types)
            unknown_resnames.update(stats['unknown_resnames'])
            if parsed is None:
                print(f'[ERROR] Could not parse {original_pdbpath}')
                continue

            # Select all ligand residues. Match on the seg/chain/resnum arrays rather than
            # building a selection string: a blank segname or chid, or a negative resnum,
            # can't be round-tripped through ProDy's selection syntax.
            segnames, chids, resnums = parsed.getSegnames(), parsed.getChids(), parsed.getResnums()
            lig_mask = np.zeros(len(parsed), dtype=bool)
            for _seg, _ch, _resnum in list_ligs:
                lig_mask |= (segnames == _seg) & (chids == _ch) & (resnums == _resnum)
            if not lig_mask.any():
                print(f'[ERROR] No ligand atoms found in {original_pdbpath} for {list_ligs}')
                continue
            ligs = parsed[np.flatnonzero(lig_mask)]

            # Select sphere around all ligand residues (entire residues)
            around = parsed.select(f'within {radius} of ligobj', ligobj=ligs)
            # Drop residues that would provoke prepwizard in s02 (see _prep_filters).
            keep_resindices = drop_prepwizard_hazard_residues(
                parsed, set(around.getResindices()))
            resindices_around_lig = ' '.join(str(i) for i in keep_resindices)
            around_lig = parsed.select(f'resindex {resindices_around_lig}')

            # Write out PDB
            if not os.path.isdir(output_subdir):
                os.makedirs(output_subdir)
            write_pdb_with_remap_remarks(output_path, around_lig, chain_mapping)
            print(output_path)
        except ChainIdOverflow as e:
            print(f'[ERROR] Skipping {output_path}: {e}')
        except Exception as e:
            print('--------------------------------------')
            print(f'[ERROR] Failed to output {output_path}')
            print(e)
            traceback.print_exc()
            print('--------------------------------------')
    if unknown_resnames:
        print(f'[WARNING] {len(unknown_resnames)} HETATM resname(s) not in the CCD type '
              f'table were left as written (a chain-bonded modified amino acid among '
              f'them would be mined as a ligand); refresh with '
              f'scripts/fetch_ccd_polymer_types.py: {sorted(unknown_resnames)}')
    # Check number of PDBs that were actually output
    num_output_pdbs = 0
    for pdb_output_subdir in os.listdir(target_dir):
        for pdb_file in os.listdir(os.path.join(target_dir, pdb_output_subdir)):
            if pdb_file.endswith('.pdb') or pdb_file.endswith('.pdb.gz'):
                num_output_pdbs += 1
    print('Number of deduplicated PDBs processed:', len(bindingsite_dict.keys()))
    print('Number of PDBs successfully output:', num_output_pdbs)





    print('Script successfully completed.')


if __name__ == "__main__":
    main()
