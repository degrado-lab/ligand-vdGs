import numpy as np
import prody as pr
from ligand_vdgs.functions.redundancy import check_networks, check_pdbnames


def get_nr_res_interactions_with_ligs_in_pdb(ligand, atoms, dist_cutoff=4.8,
                                             lig_bfactor_cutoff=40):
    '''
    Identifies the ligands and the *protein* residues interacting with those ligands in a
    single PDB. Unrefined redundancy refers to a crude/quick and dirty way of measuring
    redundancy. Later on in the vdg creation process, there will be an opportunity to refine
    redundancy checks. This quick/dirty approximation is to reduce the size of the database in
    the first place (for prepwizard, etc.)

    This is the database-*selection* step and deliberately not the contact gate. It
    asks only "does this ligand have any protein near it at all", on the whole
    ligand, before anything is prepared; the vdG membership criterion is buried CG
    surface (`ligand_vdgs.functions.sasa`), measured per chemical group after
    preparation. Do not collapse the two: loosening this cutoff changes which
    structures enter the database, not which residues enter a vdG.

    dist_cutoff refers to distance between ligand heavy atom and protein heavy atom.

    Ligands with no protein residue within dist_cutoff are dropped: they have no vdG to
    contribute.
    '''

    ligand_networks = {} # key = 'lig segment, lig chain, lig resnum, lig resname'
                         # value = ligand "networks", where each "network" is a list of
                         # interacting residues (comprising vdg). Each interacting res is
                         # represented by its (segment, chain, resnum, resname)

    ligand_resinds = set(ligand.getResindices())
    lig_residues = [] # each element is (lig segment, lig chain, lig resnum, lig resindex)
    # Collect and sort all lig residues before processing any of them, so that which copy of a
    # ligand is kept doesn't depend on the order resindices happen to appear in. Otherwise, in a
    # single PDB, the first resind of ligand CAH might be chain B, so we put chain B in the nr
    # dict and deem chain A's CAH redundant, while the first resind of HEM is chain A, so we put
    # HEM's chain A in the nr dict and deem HEM's chain B irrelevant; we'd write out both chains
    # A + B even though they're equivalent (see 9c/4c9n.pdb). Happens fairly often.

    for lig_resind in sorted(ligand_resinds):
        lig_obj = atoms.select(f'resindex {lig_resind}')
        if lig_obj is None:
            continue
        # `atoms` is already hydrogen-stripped by the caller, so lig_obj is heavy atoms only.
        # Remove instances where the lig is a single heavy atom (like O)
        if len(lig_obj) == 1:
            continue
        # Remove instances where the average b-factor of the ligand is > the lig_bfactor_cutoff.
        lig_avg_b = np.mean(lig_obj.getBetas())
        if lig_avg_b > lig_bfactor_cutoff:
            continue

        firstligatom = lig_obj[0]
        lig_residues.append((firstligatom.getSegname(), firstligatom.getChid(),
                             int(firstligatom.getResnum()), lig_resind))

    for lig_seg, lig_chain, lig_resnum, lig_resind in sorted(lig_residues):
        # Select by resindex rather than by seg/chain/resnum: a blank segname or chid, or a
        # negative resnum, can't be round-tripped through a selection string.
        lig_pr_obj = atoms.select(f'resindex {lig_resind}')
        lig_resname = lig_pr_obj.getResnames()[0]
        candidate_lig_id = ' '.join([lig_seg, lig_chain, str(lig_resnum), lig_resname])

        # Determine which segments, chains, and residues each lig interacts with.
        # 'protein' is the complement of the ligand selection in add_pdb_to_nr_db_dict
        # ('not (water or ion or protein)'), so ligand and neighbours stay disjoint by
        # construction (exwithin is belt-and-braces). Note this also excludes nucleic acids
        # and any residue ProDy's 'protein' flag doesn't recognize.
        candidate_interacting_residues = []
        neighbs = atoms.select(f'protein and exwithin {dist_cutoff} of ligobj',
                               ligobj=lig_pr_obj)
        if neighbs is None: # no protein contacts -> nothing to mine; drop the ligand
            continue
        for neighb in neighbs:
            seg, chain, resnum, resname = neighb.getSegname(), neighb.getChid(), \
                neighb.getResnum(), neighb.getResname()
            vdm_id = [str(seg), str(chain), int(resnum), str(resname)] # numpy int and str not recognized by json
            if vdm_id not in candidate_interacting_residues:
                candidate_interacting_residues.append(vdm_id)

        # Does the vdg/network (a ligand and its vdms) have the same resnums and resnames as
        # others within this same PDB? If yes, then redundant.
        is_redundant = False
        for already_deduplicated_ligs, already_deduplicated_network in ligand_networks.items():
            # Keys are ' '-joined 'seg chain resnum resname'; split to compare
            # fields (indexing the string would compare single characters). Split
            # from the right, and only twice: a blank segname or chid contributes
            # an empty field, and a segname may even contain a space, but resnum
            # and resname never do.
            _, dedup_resnum, dedup_resname = already_deduplicated_ligs.rsplit(' ', 2)
            # Check if lig resnums and resname match
            if lig_resnum == int(dedup_resnum) and lig_resname == dedup_resname:
                # If yes, check if vdm resnums and resnames match. Tolerate 2 "missing" residues
                # because this is checking within the same PDB (intra-pdb)
                is_redundant = check_networks(candidate_interacting_residues,
                                              already_deduplicated_network, tol=2)
                if is_redundant:
                    break
        if not is_redundant:
            ligand_networks[candidate_lig_id] = candidate_interacting_residues

    return ligand_networks


def add_pdb_to_nr_db_dict(database_dict, pdbpath, lig_bfactor_cutoff, atoms=None):
    '''
    Select each ligand and determine what chains it interacts with. 

    *atoms* is the parsed structure; when None it is parsed from pdbpath. s01 passes
    it pre-parsed so that chain IDs are remapped first (see preprocessing/_chain_ids).
    
    Store the ligs (seg/chain/resnum) and interacting residues in a dict to further
    reduce the database size by making a guess at whether monomers within a pdb are redundant
    by looking up the lig and vdm resnums and seeing if they're the same across the different 
    monomers (intra-pdb redundancy). You can additionally make a guess about whether 2 pdbs are
    redundant by seeing if their acc. codes are similar (i.e. in the same series), and their ligs 
    and vdms are on the same chains/resnums (inter-pdb redundancy).
    
    Nested database_dict format: lig resname:
                                   (pdbfile, (lig seg, lig chain, lig resnum, ligresname)): 
                                      nested list of interacting res (seg, chain, rensum, resname)
    
    '''

    pdbfile = pdbpath.split('/')[-1]
    if not (pdbfile.endswith('.pdb') or pdbfile.endswith('.pdb.gz')):
        print('NOT A PDB:', pdbpath)
        return database_dict
    
    # Identify ligand(s)
    if atoms is None:
        atoms = pr.parsePDB(pdbpath)
    if atoms is None:
        print('[ERROR] ProDy could not parse:', pdbpath)
        return database_dict
    atoms = atoms.select('not element H D')
    ligands = atoms.select('not (water or ion or protein)') # prody "hetero" selection not suitable
    if ligands is None: # could be None if the "ligand" is actually a noncanonical AA
        return database_dict
    
    pdb_dict = get_nr_res_interactions_with_ligs_in_pdb(
        ligands, atoms, lig_bfactor_cutoff=lig_bfactor_cutoff)
    # Add every vdg instance (lig + its interacting residues) to database_dict if not redundant
    for candidate_lig_res_id, candidate_network in pdb_dict.items():
        lig_resname = candidate_lig_res_id.split(' ')[-1]
        # Is this vdg (ligand + its interacting residues) redundant to something already in the
        # database dict?
        vdg_is_redundant = False
        existing_ligs = database_dict.setdefault(lig_resname, {})
        for deduplicated_lig_res, deduplicated_lig_network in existing_ligs.items():
            if check_networks(candidate_network, deduplicated_lig_network):
                # Are the PDBs part of the same series?
                if check_pdbnames(pdbfile, deduplicated_lig_res.split(' ')[0]):
                    vdg_is_redundant = True
                    break

        if not vdg_is_redundant:
            existing_ligs[' '.join([pdbfile, candidate_lig_res_id])] = candidate_network

    return database_dict
