import os
import io
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import Counter
import prody as pr
from ligand_vdgs.functions import frag_enumeration

KEY_SCHEMA = 'annot-1'
_LIGAND_EXCLUDE_CLAUSE = 'ion or water or resname SEP or resname TPO or resname MSE'

def _heavy_element_counts(sel):
    heavy = sel.select("not element H D") if sel is not None else None
    return Counter(el.upper() for el in heavy.getElements()) if heavy is not None else Counter()

def _template_element_counts(template, n=1):
    counts = Counter(atom.GetSymbol().upper() for atom in template.GetAtoms())
    return counts if n == 1 else Counter({el: c * n for el, c in counts.items()})

def _disambiguate_ligand_pool(pool, query_lig_smiles):
    resnames = sorted(set(pool.getResnames()))
    if len(resnames) == 1:
        return pool, resnames
    template = Chem.MolFromSmiles(query_lig_smiles)
    if template is None:
        return pool, resnames
    matches = []
    for name in resnames:
        sel = pool.select(f"resname {name}")
        if sel is None:
            continue
        n_residues = len(set(sel.getResindices().tolist()))
        if n_residues > 0 and _heavy_element_counts(sel) == _template_element_counts(template, n_residues):
            matches.append(name)
    if len(matches) == 1:
        return pool.select(f"resname {matches[0]}"), matches
    return pool, resnames

def identify_ligand_selection(query_struct, query_lig_smiles):
    hetatms = query_struct.select(f'hetatm and not ({_LIGAND_EXCLUDE_CLAUSE})')
    if hetatms is not None:
        sel, resnames = _disambiguate_ligand_pool(hetatms, query_lig_smiles)
        if len(resnames) == 1:
            return sel, resnames[0]

    all_pool = query_struct.select(f'not (protein or nucleic or {_LIGAND_EXCLUDE_CLAUSE})')
    if all_pool is None:
        raise ValueError("get_query_ligand_mol: could not find ligand atoms.")
    sel, resnames = _disambiguate_ligand_pool(all_pool, query_lig_smiles)
    if len(resnames) != 1:
        raise ValueError(
            f"get_query_ligand_mol: expected exactly one ligand residue, got: {resnames}")
    return sel, resnames[0]

def classify_ligand_instances(query_struct, query_lig_smiles):
    sel, resname = identify_ligand_selection(query_struct, query_lig_smiles)
    template = Chem.MolFromSmiles(query_lig_smiles)
    if template is None:
        return resname, []
    resindices = sorted(set(sel.getResindices().tolist()))
    instances = []
    for resindex in resindices:
        res = sel.select(f"resindex {resindex}")
        if res is not None and _heavy_element_counts(res) == _template_element_counts(template):
            instances.append((resindex,
                f"{res.getChids()[0]}{int(res.getResnums()[0])}{res.getIcodes()[0]}".strip()))
    if len(instances) != len(resindices):
        return resname, []
    return resname, instances

def select_struct_for_ligand_instance(query_struct, resname, other_resindices):
    if not other_resindices:
        return query_struct
    return query_struct.select(
        f"not (resname {resname} and resindex {' '.join(str(i) for i in other_resindices)})")

def get_query_ligand_mol(query_struct_or_path, query_lig_smiles):
    if isinstance(query_struct_or_path, str):
        if query_struct_or_path.endswith(('.pdb', '.pdb.gz')):
            query_struct = pr.parsePDB(query_struct_or_path)
        elif query_struct_or_path.endswith(('.cif', '.cif.gz')):
            query_struct = pr.parseCIF(query_struct_or_path)
        else:
            raise ValueError(f"Unsupported file format: {query_struct_or_path}")
    else:
        query_struct = query_struct_or_path

    query_hetatms, ligname = identify_ligand_selection(query_struct, query_lig_smiles)

    n_query_ligand_residues = len(set(query_hetatms.getResindices().tolist()))
    if n_query_ligand_residues != 1:
        raise ValueError(
            f"get_query_ligand_mol: expected exactly one ligand "
            f"residue instance of {ligname!r}, got {n_query_ligand_residues} "
            "residues.")

    buf = io.StringIO()
    pr.writePDBStream(buf, query_hetatms)
    query_pdb_mol = Chem.MolFromPDBBlock(buf.getvalue(), removeHs=True)
    query_lig_template = Chem.MolFromSmiles(query_lig_smiles)
    if query_pdb_mol is None or query_lig_template is None:
        print(f"[ERROR] Processing {query_struct_or_path}: RDKit failed to parse molecule",
              flush=True)
        return None

    if query_pdb_mol.GetNumAtoms() != query_lig_template.GetNumAtoms():
        print(f"[WARNING] Processing {query_struct_or_path}: ligand selection has "
              f"{query_pdb_mol.GetNumAtoms()} heavy atoms but the query SMILES has "
              f"{query_lig_template.GetNumAtoms()}; bond orders are assigned only to the "
              "atoms the template matches, and the rest stay single-bonded.",
              flush=True)

    try:
        stripped = frag_enumeration.manually_remove_Hs(
            AllChem.AssignBondOrdersFromTemplate(query_lig_template, query_pdb_mol),
            return_single_mol_or_perms='single')
        if stripped is None:
            print(f"[ERROR] Processing {query_struct_or_path}: could not strip hydrogens "
                  "from the ligand (see the warning above)", flush=True)
            return None
        (query_pdb_mol_assigned_bonds_no_H, _pdb_mol_perm_inds), _smiles_no_H = stripped
        Chem.GetSymmSSSR(query_pdb_mol_assigned_bonds_no_H)
        return query_pdb_mol_assigned_bonds_no_H
    except Exception as e:
        print(f"[ERROR] Processing {query_struct_or_path}: "
              f"ligand hydrogen removal failed: {type(e).__name__}: {e}", flush=True)
        return None

def submol_from_match(mol, match):
    conf_src = mol.GetConformer()
    em = Chem.RWMol()
    conf = Chem.Conformer(len(match))
    for new_idx, orig_idx in enumerate(match):
        em.AddAtom(Chem.Atom(mol.GetAtomWithIdx(orig_idx)))
        conf.SetAtomPosition(new_idx, conf_src.GetAtomPosition(orig_idx))
    positions = {orig_idx: new_idx for new_idx, orig_idx in enumerate(match)}
    for i, orig_i in enumerate(match):
        for bond in mol.GetAtomWithIdx(orig_i).GetBonds():
            j = positions.get(bond.GetOtherAtomIdx(orig_i))
            if j is not None and j > i:
                em.AddBond(i, j, bond.GetBondType())
    sub = em.GetMol()
    sub.AddConformer(conf, assignId=True)
    return sub

def check_vdg_job_status(sub_smiles, vdg_lib_dir):
    vdg_log_file = os.path.join(vdg_lib_dir, sub_smiles, f'{sub_smiles}_log')
    if not os.path.exists(vdg_log_file):
        return False
    with open(vdg_log_file, 'r') as f:
        return 'Job completed.' in f.read()

def summarize_frags(frags_in_lib, frags_to_exclude, frags_to_include, logfile_fh):
    groups = {
        "Excluded": [],
        "Not in include list": [],
        "In vdg lib but incomplete": [],
        "In vdg lib and in include list": []}

    for sub_smiles in frags_in_lib:
        if sub_smiles in frags_to_exclude:
            groups["Excluded"].append(sub_smiles)
        elif (frags_to_include != 'all' and isinstance(frags_to_include, list) and
              sub_smiles not in frags_to_include):
            groups["Not in include list"].append(sub_smiles)
        elif not frags_in_lib[sub_smiles]:
            groups["In vdg lib but incomplete"].append(sub_smiles)
        else:
            groups["In vdg lib and in include list"].append(sub_smiles)

    print("\n--- Fragment vdG Library Summary ---", file=logfile_fh)
    for category, frags in groups.items():
        if frags:
            print(f"{category}:", file=logfile_fh)
            print("   ", ", ".join(frags), file=logfile_fh)
    print("\n", file=logfile_fh)
