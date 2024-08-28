'''
Given a query structure (a crystal structure of a protein-ligand complex) and its binding site 
residues, determine whether any vdGs can correctly place the FGs within the ligand of the query
structure.
'''

from rdkit import Chem

query_path = '/wynton/home/degradolab/skt/docking/query_structs/1acj__A_999_1.pdb'
smarts = 'n1ccccc1'

bindingsite_residues = []
bindingsite_residues.append(['', 'A', 84])
bindingsite_residues.append(['', 'A', 440])
bindingsite_residues.append(['', 'A', 330])
bindingsite_residues.append(['', 'A', 432])


# Try all singles and pairs of residues in the binding site
