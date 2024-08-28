from rdkit import Chem


def find_matches(prody_pkl, patt, mols, molnames, validation_df, pdb_dataset_path, 
                 pdb_ligs_dict, outputfile, debug, halt_early=False):
    """Find ligands that match a SMARTS pattern in a pickled ProDy object.
    patt may also be re-defined depending on if it needed to be converted to a different form.

    Parameters
    ----------
    prody_pkl : str
        Path to pickled ProDy object.
    patt : rdkit.Chem.rdmol.Mol
        Mol object generated from a SMARTS pattern by RDKit.
    mols : list
        List of RDKit Mol objects that may match the SMARTS pattern.
    molnames : list
        List of three-letter accession codes for the matching ligands.
    validation_df : pandas.DataFrame
        Dataframe containing validation information on all protein 
        structures in the PDB.
    halt_early : bool, optional
        If true, halt after the first match is found.

    Returns
    -------
    match_mols : list
        List of RDKit Mol objects extracted from the pickled ProDy object 
        that match the SMARTS pattern.
    match_molnames : list
        List of the three-letter ligand names of each matching Mol object.
    segs_chains_resnums : list
        List containing tuples of segments, chains, and resnums for each 
        ligand in the pickled Prody object that matches the SMARTS pattern.
    new_patts : list
        Re-do the patt Mol object for each match_mol because in some cases, it's
        converted to smarts/smiles and then converted back to a Mol object in order
        for HasSubstruct() to find the match.
    """
    
    pdbpath = os.path.join(pdb_dataset_path, f'{pdb_acc.lower()}.pdb')
    struct = pr.parsePDB(pdbpath)

    match_mols = []
    match_molnames = []
    segs_chains_resnums = []
    new_patts = []
    ligs = [lig for lig in ligs if lig != '']
    selection = '(resname ' + ' or resname '.join(ligs) + \
                ') within 5 of protein' 
    all_segs_chains_resnums = get_segs_chains_resnums(struct, selection)
    for segname, chain, resnum in all_segs_chains_resnums:
        try:
            if segname:
                sel = struct.select(('segment {} and chain {} '
                                     'and resnum {}').format(
                                    segname, chain, resnum))
            else:
                sel = struct.select('chain {} and resnum {}'.format(
                                    chain, resnum))
            resname = sel.getResnames()[0]
            resname_idx = molnames.index(resname)
            template = mols[resname_idx]
            match_molname = molnames[resname_idx]
            pdbio = StringIO()
            pr.writePDBStream(pdbio, sel, secondary=False)
            
            block = pdbio.getvalue()
            mol = Chem.rdmolfiles.MolFromPDBBlock(block, removeHs=False, 
                                                  sanitize=False)
            

            has_substruct_bool, new_patt = has_substruct(
                mol, patt, match_molname, pdb_acc)
            if has_substruct_bool:
                match_mols, match_molnames, segs_chains_resnums, new_patts = update_list_of_matches(
                    mol, match_molname, segname, chain, resnum, match_mols, match_molnames,
                    segs_chains_resnums, new_patts, new_patt)
            else: # Try to sanitize mol
                Chem.SanitizeMol(mol)
                has_substruct_bool, new_patt = has_substruct(
                    mol, patt, match_molname, pdb_acc)
                if has_substruct_bool:
                    match_mols, match_molnames, segs_chains_resnums, new_patts = update_list_of_matches(
                        mol, match_molname, segname, chain, resnum, match_mols, match_molnames,
                        segs_chains_resnums, new_patts, new_patt)
                else: # Assign bond orders (sometimes necessary for detecting aromaticity)
                    mol = assign_bond_orders_from_template(template, mol)

                    has_substruct_bool, new_patt = has_substruct(
                        mol, patt, match_molname, pdb_acc)
                    if has_substruct_bool:
                        match_mols, match_molnames, segs_chains_resnums, new_patts = update_list_of_matches(
                            mol, match_molname, segname, chain, resnum, match_mols, match_molnames,
                            segs_chains_resnums, new_patts, new_patt)
                    else: # Likely does not have a match because of unresolved density
                        #print(f'{match_molname} in {pdbpath} does not have a substruct match.')
                        with open(outputfile, 'a') as outF: # write in output file
                            outF.write('======================================= \n')
                            outF.write(f'{match_molname} in {pdbpath} does not have a substruct match. \n')
                        continue # to next lig
                
        except Exception as e:
            print('Exception raised for ', pdbpath, match_molname)
            with open(outputfile, 'a') as outF: # write in output file
                outF.write('======================================= \n')
                outF.write('Exception raised for {} {} \n'.format(pdbpath, match_molname))
                outF.write('smiles: \n')
                outF.write(Chem.MolToSmiles(mol)+'\n')
                outF.write('smarts: \n')
                outF.write(Chem.MolToSmarts(mol)+'\n')
                outF.write(f'{e} \n')
            pr.writePDB(f'errored_{pdb_acc}.pdb', struct)

            # Add atom indices to mol for outputting png of failed ligand: 
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx())
            img = Draw.MolToImage(mol)
            img.save(f'{pdb_acc}.png')
            # If debug mode, terminate program to investigate. Otherwise, pass
            if debug:
                raise Exception(e)
            else: 
                pass

        if halt_early and len(match_mols):
            break
    

    if len(match_mols) == 0:
        pr.writePDB(f'errored_{pdb_acc}.pdb', struct)
        if debug:
            raise ValueError(f'PDB file {pdbpath} was not found to contain a smarts-containing ligand although it does; investigate why.') # We know it must contain our smarts because it's in pdb_ligs_dict.

    #if new_patt is None: # append a unique new_patt for each mol
    #    patt_none_error = f'patt is None in find_matches() for PDB {pdbpath}'
    #    print(patt_none_error)
    #    with open(outputfile, 'a') as outF: # write in output file
    #        outF.write('======================================= \n')
    #        outF.write(patt_none_error + '\n')
    return match_mols, match_molnames, segs_chains_resnums, new_patts


