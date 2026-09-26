"""Check current-schema vdG buckets and fragment symmetry sidecars for 
   checking the vdg library after it's built."""
import json, multiprocessing, os, sys
import numpy as np
from ligand_vdgs.functions import parent_db, utils, vdg_npz_utils as npz

H_CLASS_FIELDS = ('nr_cg_heavy_degree', 'nr_cg_num_h', 'mem_cg_heavy_degree', 'mem_cg_num_h')

def contract_problems(data, npz_path):
    """Shape/dtype/sentinel problems in a bucket's per-observation H-class fields."""
    n_cg, n_nr, n_mem = len(data['cg_elements']), len(data['cluster_id']), len(data['mem_cluster_id'])
    problems = []
    for field in H_CLASS_FIELDS:
        shape = (n_nr if field.startswith('nr_') else n_mem, n_cg)
        if field not in data:
            problems.append(f'{npz_path}: missing field {field}; expected int8 {shape}')
            continue
        arr = data[field]
        if arr.shape != shape:
            problems.append(f'{npz_path}: {field} has shape {arr.shape}, expected {shape}')
        if not np.issubdtype(arr.dtype, np.integer):
            problems.append(f'{npz_path}: {field} has dtype {arr.dtype}, expected an integer dtype')
        elif arr.dtype != np.int8:
            print(f'[WARNING] {npz_path}: {field} is {arr.dtype}, not int8; read anyway.', file=sys.stderr)
    for prefix in ('nr', 'mem'):
        deg, nh = f'{prefix}_cg_heavy_degree', f'{prefix}_cg_num_h'
        if deg not in data or nh not in data or data[deg].shape != data[nh].shape:
            continue
        if n_zero := int((data[deg] == 0).sum()):
            problems.append(f'{npz_path}: {deg} has {n_zero} entries == 0; '
                            'the contract requires a NEGATIVE sentinel for unreadable atoms.')
        if n_rev := int(((data[deg] < 0) & (data[nh] >= 0)).sum()):
            print(f'[WARNING] {npz_path}: {n_rev} entries have {deg} < 0 but {nh} >= 0.', file=sys.stderr)
    return problems

def inspect_library(lib, fragments, zero_match=(), workers=1):
    """Check every fragment in parallel; return (bucket count, sorted errors)."""
    with multiprocessing.Pool(workers) as pool:
        results = pool.starmap(inspect_fragment,
                               [(lib, frag, frag in zero_match) for frag in fragments], chunksize=1)
    errors = sorted(e for _, _, frag_errors, _ in results for e in frag_errors)
    provenance = {}
    for frag, _, _, source in results:
        if source is not None:
            provenance.setdefault(json.dumps(source, sort_keys=True, default=str), []).append(frag)
    if len(provenance) > 1:
        groups = sorted(sorted(frags) for frags in provenance.values())
        errors.append(f'provenance differs across fragments ({len(groups)} variants, '
                      f'e.g. {groups[0][0]} vs {groups[1][0]})')
    buckets = sum(n for _, n, _, _ in results)
    if not buckets:
        errors.append('no buckets checked')
    return buckets, errors

def inspect_fragment(lib, frag, zero_match):
    """Check one fragment; return (frag, bucket count, errors, provenance)."""
    errors, provenance, buckets = [], None, 0
    root = os.path.join(lib, frag, 'nr_vdgs')
    if zero_match:
        if os.path.exists(root):
            errors.append(f'{frag}: zero-match fragment has nr_vdgs output')
        return frag, 0, errors, None
    try:
        smarts, perms = npz.load_cg_symmetry(lib, frag)
        n_cg = utils.fragment_key_query_mol(smarts).GetNumAtoms()
        if utils.smiles_to_filename(smarts) != frag:
            raise ValueError('stored SMARTS does not match directory label')
        utils.validate_atom_permutations(perms, n_cg)
    except (OSError, ValueError, TypeError, AttributeError, KeyError) as exc:
        errors.append(f'{frag}: symmetry: {exc}')
        n_cg = None
    if not os.path.isdir(root):
        errors.append(f'{frag}: missing nr_vdgs')
        return frag, 0, errors, None
    for size in os.listdir(root):
        size_path = os.path.join(root, size)
        if size == 'cg_symmetry.npz':
            continue
        if size not in ('1', '2') or not os.path.isdir(size_path):
            errors.append(f'{frag}: unexpected nr_vdgs entry {size}')
            continue
        for sign in os.listdir(size_path):
            sign_path = os.path.join(size_path, sign)
            if sign not in npz.CHARGE_SIGNS or not os.path.isdir(sign_path):
                errors.append(f'{frag}: unexpected charge entry {size}/{sign}')
                continue
            for name in os.listdir(sign_path):
                path = os.path.join(sign_path, name)
                if not name.endswith('.npz') or not os.path.isfile(path):
                    errors.append(f'{frag}: unexpected bucket entry {size}/{sign}/{name}')
                    continue
                buckets += 1
                try:
                    # NpzFile decompresses on every access, so read each array once.
                    with np.load(path, allow_pickle=False) as handle:
                        data = {key: handle[key] for key in handle.files}
                    source = check_bucket(data, path, int(size), sign, name[:-4], n_cg)
                    source.pop('build_date', None)
                    if provenance is None:
                        provenance = source
                    elif source != provenance:
                        raise ValueError('provenance differs from another bucket')
                except npz.CORRUPT_NPZ_ERRORS + (IndexError,) as exc:
                    errors.append(f'{path}: {exc}')
    return frag, buckets, errors, provenance

def check_bucket(data, path, size, sign, name, n_cg):
    required = ('schema', 'aa_bucket_parts', 'charge_sign', 'parent_pdb_dir',
                'cluster_id', 'cluster_size', 'mem_cluster_id', 'cluster_num_parents',
                'first_stage_cluster_id', 'second_stage_cluster_id', 'cluster_pose_radius',
                'nr_parent_biounit', 'mem_parent_biounit', 'cg_elements', 'nr_cg_coords',
                'nr_vdm_bb_coords', 'nr_vdm_o_coords', 'nr_cg_names', 'mem_cg_names',
                'nr_scrr_resname', 'mem_scrr_resname', 'nr_slot_flag')
    missing = set(required) - set(data)
    if missing:
        raise ValueError(f'missing schema fields: {sorted(missing)}')
    schema = json.loads(str(data['schema']))
    if schema['schema_version'] != npz.BUCKET_SCHEMA_VERSION:
        raise ValueError('wrong bucket schema version')
    parts = [str(x) for x in data['aa_bucket_parts']]
    if len(parts) != size or npz.make_aa_bucket(parts) != name:
        raise ValueError('AA bucket parts disagree with path')
    if str(data['charge_sign']) != sign:
        raise ValueError('charge sign disagrees with path')
    if str(data['parent_pdb_dir']) != schema['parent_pdb_dir']:
        raise ValueError('parent PDB directory disagrees with provenance')
    ids, sizes = data['cluster_id'], data['cluster_size']
    members = data['mem_cluster_id']
    count = len(ids)
    if not count or not np.array_equal(ids, np.arange(1, count + 1)):
        raise ValueError('cluster IDs are empty or nonconsecutive')
    if np.any(sizes < 1) or len(members) != int(np.sum(sizes - 1)):
        raise ValueError('cluster sizes disagree with member count')
    if len(members) and (np.min(members) < 1 or np.max(members) > count):
        raise ValueError('member cluster ID outside cluster range')
    if not np.array_equal(np.bincount(members, minlength=count + 1)[1:], sizes - 1):
        raise ValueError('member cluster IDs disagree with cluster sizes')
    for prefix, rows in (('nr_', count), ('mem_', len(members))):
        for key, value in data.items():
            if key.startswith(prefix) and (not value.shape or value.shape[0] != rows):
                raise ValueError(f'{key} row count differs from {prefix} rows')
    if n_cg is not None and data['nr_cg_coords'].shape != (count, n_cg, 3):
        raise ValueError('CG coordinate shape disagrees with symmetry SMARTS')
    if data['nr_vdm_bb_coords'].shape != (count, size, 3, 3):
        raise ValueError('vdM backbone coordinate shape disagrees with path')
    if data['nr_vdm_o_coords'].shape != (count, size, 3):
        raise ValueError('vdM O coordinate shape disagrees with path')
    problems = contract_problems(data, path)
    if problems:
        raise ValueError(problems[0])
    if np.any(data['cluster_num_parents'] < 1) or np.any(data['cluster_num_parents'] > sizes):
        raise ValueError('parent support outside cluster size')
    parents = [{parent_db.entry_of(stem)} for stem in data['nr_parent_biounit']]
    for cid, stem in zip(members, data['mem_parent_biounit']):
        parents[int(cid) - 1].add(parent_db.entry_of(stem))
    if not np.array_equal([len(group) for group in parents], data['cluster_num_parents']):
        raise ValueError('parent support disagrees with biounits')
    return schema
