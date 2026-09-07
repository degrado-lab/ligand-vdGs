'''
Quick/dirty redundancy checks used while trimming the parent PDB database.

These are deliberately crude: they compare residue numbers/names rather than
geometry, to shrink the database cheaply before prepwizard/probe. Redundancy is
refined properly later, during vdG generation.
'''


def check_pdbnames(pdbfile1, pdbfile2):
    '''True if two PDB filenames share >=3 of their 4 accession characters,
    i.e. they were likely deposited as part of the same project.'''
    pdbacc1, pdbacc2 = pdbfile1[:4], pdbfile2[:4]
    num_matches = sum(1 for c1, c2 in zip(pdbacc1, pdbacc2) if c1 == c2)
    return num_matches >= 3


def check_networks(candidate_interacting_residues, already_deduplicated_network, tol=0):
    '''
    True if two interaction networks are redundant, judged by whether their
    (resnum, resname) sets coincide.

    tol is the number of "missing" residues to tolerate: use 2 for intra-PDB
    (a monomer may drop a contact), 0 for inter-PDB (a real mutant or an extra
    ligand contact should not be collapsed away).
    '''
    a = [(i[2], i[3]) for i in candidate_interacting_residues]
    b = [(i[2], i[3]) for i in already_deduplicated_network]

    num_matches = len(set(a) & set(b))
    num_vdms_between_both_sets = len(set(a) | set(b))

    return num_matches >= num_vdms_between_both_sets - tol
