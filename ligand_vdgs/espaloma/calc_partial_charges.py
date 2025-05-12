import prody as pr
from io import StringIO
import espaloma as esp
from rdkit import Chem
from openff.toolkit import topology

pdbpath = 'espaloma/tsa.pdb'
''
def calc_partial_charges(molecule, frag_names):
    # Use espaloma to calculate partial charges, given an openff molecule obj
    molecule_graph = esp.Graph(molecule)
    espaloma_model = esp.get_model("latest") # load pre-trained model
    espaloma_model(molecule_graph.heterograph) # apply the model
    # create an OpenMM System and calculate the partial charges
    try:
        openmm_system = esp.graphs.deploy.openmm_system_from_graph(molecule_graph)
        # retrieve partial charges from molecule obj
        dict_from_esp = {} # key = atom name, val = partial_charge 
        for a in molecule.atoms:
            name, partial_charge = (a.name, a.partial_charge.magnitude)
            if name not in frag_names:
                continue
            assert name not in dict_from_esp.keys()
            dict_from_esp[name] = partial_charge
        return dict_from_esp
    except Exception as e:
        print(f"Error: {e}")
        return None

def prody_to_rdmol_and_openff_mol(prody_obj):
    # Convert prody obj to an rdkit mol and openff Topology obj to run espaloma. 
    # To avoid writing and reading from disk, write to pdb stream.

    with StringIO() as pdb_stream: # write to memory buffer
        pr.writePDBStream(pdb_stream, prody_obj)
        pdb_data = pdb_stream.getvalue()

    # Convert to RDKit Mol
    rdmol = Chem.MolFromPDBBlock(pdb_data, sanitize=False, removeHs=False)
    num_atoms_in_rdmol = rdmol.GetNumAtoms()
    
    # Convert to openFF Molecule
    ff_mol = topology.Molecule.from_rdkit(rdmol, hydrogens_are_explicit=True,
                                          allow_undefined_stereo=True)
    num_atoms_in_ff_mol = len(ff_mol.atoms) 
    if num_atoms_in_rdmol != num_atoms_in_ff_mol:
        print('WARNING: this lig was processed incorrectly. The # of atoms in the rdkit '
              'mol must be the same as the # of atoms in the openff mol.')
    return rdmol, ff_mol

par = pr.parsePDB(pdbpath)
rdmol, ff_mol = prody_to_rdmol_and_openff_mol(par)
print(par.getNames())
calc_partial_charges(ff_mol, par.getNames())
