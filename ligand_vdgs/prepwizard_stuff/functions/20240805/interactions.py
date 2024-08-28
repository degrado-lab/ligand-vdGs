import numpy as np
import prody as pr
from smart_vdms.functions.redundancy import check_networks, check_pdbnames


def get_nr_res_interactions_with_ligs_in_pdb(ligand, complex, unrefined, dist_cutoff=4.8, 
                                             lig_bfactor_cutoff=40):
    '''
    Identifies the ligands and the protein residues interacting with those ligands in a single PDB.
    Unrefined redundancy refers to a crude/quick and dirty way of measuring redundancy. 
    Later on in the vdg creation process, there will be an opportunity to refine redundancy checks.
    This quick/dirty approximation is to reduce the size of the database in the first place (for prepwizard, probe, etc.) 
    
    dist_cutoff refers to distance between ligand heavy atom and protein heavy atom.

    '''
    
    ligand_networks = {} # key = (ligand segment, lig chain, lig resnum, lig resname)
                         # value = ligand "networks", where each "network" is a list of 
                         # interacting residues (comprising vdg). Each interacting res is represented
                         # by its (segment, chain, resnum, resname)

    ligand_resinds = set(ligand.getResindices())
    all_lig_residues = [] # each element is (lig segment, lig chain, lig resnum)
    # First, collect all lig res seg, chain, and resnums to avoid cases where, for example, in a single 
    # PDB, the first resind of ligand CAH is chain B, so we put chain B in the nr dict and deem chain A's 
    # CAH redundant, but the first resind of HEM is chain A, so we put HEM's chain A in the nr dict and 
    # deem HEM's chain B irrelevant; in this case, we would be writing out both chains A + B, but in 
    # actuality, they're equivalent (see 9c/4c9n.pdb). Happens fairly often.

    for lig_resind in sorted(ligand_resinds): # but don't store the lig resind b/c it will get re-indexes
                                      # in the trimmed pdb
        lig_obj = complex.select(f'resindex {lig_resind}')
        # Remove instances where the lig is a single heavy atom (like O)
        if len(lig_obj.select('not element H D')) == 1:
            continue
        # Remove instances where the average b-factor of the ligand is > the lig_bfactor_cutoff.
        lig_avg_b = np.mean(lig_obj.select('not element H D').getBetas())
        if lig_avg_b > lig_bfactor_cutoff:
            continue

        firstligatom = lig_obj[0]
        _ligseg_, _ligch_, _ligresnum_ = firstligatom.getSegname(), firstligatom.getChid(), firstligatom.getResnum()
        all_lig_residues.append(tuple([_ligseg_, _ligch_, _ligresnum_]))

    sorted_lig_residues = sorted(set(all_lig_residues))
    for ligseg_, ligch_, ligresnum_ in sorted_lig_residues:
        if ligresnum_ < 0: # need grave accent to escape negative sign
            ligresnum_ = f"`{ligresnum_}`"
        if ligseg_ == '': 
            lig_pr_obj = complex.select(f'chain {ligch_} and resnum {ligresnum_}')
        else:
            lig_pr_obj = complex.select(f'segment {ligseg_} and chain {ligch_} and resnum {ligresnum_}')

        candidate_lig_seg, candidate_lig_chain, candidate_lig_resnum, candidate_lig_resname = \
            lig_pr_obj.getSegnames()[0], lig_pr_obj.getChids()[0], lig_pr_obj.getResnums()[0], \
            lig_pr_obj.getResnames()[0]
        candidate_lig_id = [candidate_lig_seg, candidate_lig_chain, candidate_lig_resnum, 
                            candidate_lig_resname]
        
        # Determine which segments, chains, and residues each lig interacts with
        candidate_interacting_residues = []
        neighbs = complex.select(f'within {dist_cutoff} of ligobj', ligobj=lig_pr_obj)
        for neighb in neighbs:
            seg, chain, resnum, resname = neighb.getSegname(), neighb.getChid(), neighb.getResnum(), \
                neighb.getResname()
            vdm_id = tuple([seg, chain, resnum, resname])
            if vdm_id not in candidate_interacting_residues:
                candidate_interacting_residues.append(vdm_id)

        # Does the vdg/network (a ligand and its vdms) have the same resnums and resnames as others
        # within this same PDB? If yes, then redundant.  

        if len(ligand_networks.keys()) == 0:
            ligand_networks[tuple(candidate_lig_id)] = tuple(candidate_interacting_residues)
        else: 
            candidate_network = {}
            is_redundant = False
            for already_deduplicated_ligs, already_deduplicated_network in ligand_networks.items():
                # Check if lig resnums and resname match
                if (candidate_lig_resnum == already_deduplicated_ligs[2] and 
                    candidate_lig_resname == already_deduplicated_ligs[3]):
                    # If yes, check if vdm resnums and resnames match. Tolerate 2 "missing" residues
                    # because this is checking within the same PDB (intra-pdb)
                    is_redundant = check_networks(candidate_interacting_residues, 
                                                  already_deduplicated_network, tol=2)
                    if is_redundant:
                        break
            if not is_redundant:
                candidate_network[tuple(candidate_lig_id)] = tuple(candidate_interacting_residues)
            # Update the dict after iteration
            ligand_networks.update(candidate_network)

    return ligand_networks


def add_pdb_to_nr_db_dict(database_dict, pdbpath, unrefined, lig_bfactor_cutoff):
    '''
    Select each ligand and determine what chains it interacts with. 
    
    Store the ligs (seg/chain/resnum) and interacting residues in a dict to further
    reduce the database size by making a guess at whether monomers within a pdb are redundant
    by looking up the lig and vdm resnums and seeing if they're the same across the different 
    monomers (intra-pdb redundancy). You can additionally make a guess about whether 2 pdbs are
    redundant by seeing if their acc. codes are similar (i.e. in the same series), and their ligs 
    and vdms are on the same chains/resnums (inter-pdb redundancy).
    
    Nested database_dict format: lig resname:
                                   (pdbfile, (lig seg, lig chain, lig resnum, ligresname)): 
                                      nested tuple of interacting res (seg, chain, rensum, resname)
    
    '''

    pdbfile = pdbpath.split('/')[-1]
    if not pdbfile.endswith('pdb'):
        print('NOT A PDB:', pdbpath)
        return database_dict
    
    # Identify ligand(s)
    complex = pr.parsePDB(pdbpath).select('not element H D')
    ligands = complex.select('not (water or ion or protein)') # prody "hetero" selection not suitable
    if ligands is None: # could be None if the "ligand" is actually a noncanonical AA
        return database_dict
    
    pdb_dict = get_nr_res_interactions_with_ligs_in_pdb(ligands, complex, unrefined, 
                                                        lig_bfactor_cutoff=lig_bfactor_cutoff)
    # Add every vdg instance (lig + its interacting residues) to database_dict if not redundant
    for candidate_lig_res_id, candidate_network in pdb_dict.items():
        lig_resname = candidate_lig_res_id[-1]
        # Is this vdg (ligand + its interacting residues) redundant to something already in the database dict?
        vdg_is_redundant = False
        
        if len(database_dict) == 0:
            database_dict[lig_resname] = {}
            database_dict[lig_resname][(pdbfile, candidate_lig_res_id)] = candidate_network
        else:
            if lig_resname not in database_dict.keys():
                database_dict[lig_resname] = {}
                database_dict[lig_resname][(pdbfile, candidate_lig_res_id)] = candidate_network
            else: # Need to check against existing networks, which is quite deep in the nested structure
                for deduplicated_lig_res, deduplicated_lig_network in database_dict[lig_resname].items():
                    network_is_redundant = check_networks(candidate_network, deduplicated_lig_network)
                    if network_is_redundant:
                        # Are the PDBs part of the same series?
                        same_pdb_series = check_pdbnames(pdbfile, deduplicated_lig_res[0])
                        if same_pdb_series:
                            vdg_is_redundant = True
                            break

        if not vdg_is_redundant:
            database_dict[lig_resname][(pdbfile, candidate_lig_res_id)] = candidate_network

    return database_dict




