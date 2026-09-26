"""Check a completed vdG library."""
# Reads every bucket once; wall time scales with 1/--workers.

import argparse, os
from ligand_vdgs.functions import utils
from ligand_vdgs.functions.ligand_structure import check_vdg_job_status
from ligand_vdgs.generate_vdgs.extract_fragment_smiles import load_frags_dict
from ligand_vdgs.generate_vdgs.post_build_buckets import inspect_library

def library_dirs(lib):
    return sorted(d for d in os.listdir(lib) if not d.startswith('.') and os.path.isdir(os.path.join(lib, d)))

def check_aliases(lib, dict_keys, fragments):
    path = os.path.join(lib, 'fragment_aliases.tsv')
    if not os.path.isfile(path):
        return 0, [f'missing {path}']
    aliases, errors = {}, []
    with open(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 3 or not all(fields):
                errors.append(f'alias line {line_no}: expected three nonempty fields')
                continue
            alias, rep, kind = fields
            if alias in aliases:
                errors.append(f'alias line {line_no}: duplicate {alias}')
            aliases[alias] = rep
            if alias not in dict_keys:
                errors.append(f'alias line {line_no}: {alias} absent from fragment dict')
            if kind not in ('neutral_twin', 'promoted') or (rep in dict_keys) != (kind == 'neutral_twin'):
                errors.append(f'alias line {line_no}: invalid kind {kind} for {rep}')
            if utils.smiles_to_filename(rep) not in fragments:
                errors.append(f'alias line {line_no}: missing representative {rep}')
    if not aliases:
        errors.append('alias table has no rows')
    if not any('+' in key or '-' in key for key in aliases):
        errors.append('alias table has no charged alias')
    errors += [f'unknown fragment directory {d}' for d in sorted(
        set(fragments) - {utils.smiles_to_filename(k) for k in dict_keys | set(aliases.values())})]
    return len(aliases), errors

def check_completion(lib, fragments):
    errors, zero_match = [], set()
    for frag in fragments:
        path = os.path.join(lib, frag, f'{frag}_log')
        if not os.path.isfile(path):
            errors.append(f'{frag}: missing completion log')
            continue
        if not check_vdg_job_status(frag, lib):
            errors.append(f'{frag}: no completion marker')
            continue
        with open(path, errors='replace') as handle:
            lines = {line.strip() for line in handle}
            if 'Job completed.' not in lines:
                errors.append(f'{frag}: no completion marker')
            if 'No ligand matches; no vdG buckets.' in lines:
                zero_match.add(frag)
    return errors, zero_match

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vdg-lib-dir', required=True)
    parser.add_argument('--frags-dict', default='resources/database_frags_dict.pkl')
    parser.add_argument('--workers', type=int, default=int(os.environ.get('NSLOTS', 4)),
                        help='parallel processes (default: $NSLOTS, else 4)')
    args = parser.parse_args()
    lib = os.path.abspath(os.path.expanduser(args.vdg_lib_dir))
    if not os.path.isdir(lib):
        parser.error(f'no such library: {lib}')
    frags, _, _ = load_frags_dict(args.frags_dict)
    fragments = library_dirs(lib)
    checks = []
    checks.append(('fragments', [] if fragments else ['library has no fragments']))
    alias_count, alias_errors = check_aliases(lib, {key for stem_map in frags.values()
                                                   for key in stem_map}, fragments)
    checks.append(('aliases', alias_errors))
    completion_errors, zero_match = check_completion(lib, fragments)
    checks.append(('completion', completion_errors))
    buckets, bucket_errors = inspect_library(lib, fragments, zero_match, args.workers)
    checks.append(('buckets', bucket_errors))
    for name, errors in checks:
        print(f'{name}: {"FAIL" if errors else "OK"}')
        for error in errors[:10]:
            print(f'  [ERROR] {error}')
        if len(errors) > 10:
            print(f'  [ERROR] {len(errors) - 10} more')
    print(f'{len(fragments)} fragments, {alias_count} aliases, {buckets} buckets checked')
    if any(errors for _, errors in checks):
        return 1
    print('All post-build checks passed.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
