import os
import espaloma as esp
from openff.toolkit import Topology
#from openff.toolkit import topology
import openmm as mm

def calc_partial_charges(pdbpath):
    # Use espaloma to calculate partial charges, given a smiles 
    print(os.path.exists(pdbpath))
    mm_molecule = Topology.from_pdb(pdbpath)
    #mm_molecule = topology.Molecule.from_smiles(smiles)
    molecule_graph = esp.Graph(mm_molecule)
    espaloma_model = esp.get_model("latest") # load pre-trained model
    espaloma_model(molecule_graph.heterograph) # apply the model
    # create an OpenMM System 
    openmm_system = esp.graphs.deploy.openmm_system_from_graph(molecule_graph)
    list_atoms = [] # list of openff Atom objs. Append in order of atom_inds, in case 
                    # the order of molecule.atoms is inconsistent.
    for atom in mm_molecule.atoms:
        atomic_num, charge, is_arom, in_ring = (atom.atomic_number, atom.partial_charge, 
                                                atom.is_aromatic, atom.is_in_ring)

    #for force_type in openmm_system.getForces():
    #    if isinstance(force_type, mm.openmm.NonbondedForce):
    #        nonbonded = force_type
    #        break

pdbpath = '/wynton/group/degradolab/skt/docking/databases/vdg_lib/aliphatic_bromide/vdg_pdbs/6qgn__A_301_1.pdb.gz'
calc_partial_charges(pdbpath)