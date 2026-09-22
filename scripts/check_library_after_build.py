"""Post-build invariants for a vdG library. Read-only; exits non-zero on failure.

Run after every full library build or rebuild. Each check is independent and
reports its own verdict, so one failure does not hide the others.

Checks:
  aliases-resolve   every library directory is a fragment-dict key OR a
                    kind=promoted representative in fragment_aliases.tsv.
  charged-forward   every charged dict key maps forward to a directory that
                    exists. This is the direction that prevents "I searched for
                    the drawn form and found nothing".
  completeness      every directory carries 'Job completed.' in its <cg>_log.
  h-class-fields    the per-observation H-class columns the H-class diagnostic
                    needs are present (delegates to h_class_diagnostic --dry-run;
                    docs/database_generation_guide.md says to run it after every
                    build).

Fragment keys contain'[' and ']' (glob metacharacters), and 185
shipped keys strictly prefix another. Everything here uses os.listdir and maps
fragment -> directory FORWARD through smiles_to_filename. Never prefix-match.
"""
import argparse, os, subprocess, sys
from ligand_vdgs.functions import utils
from ligand_vdgs.generate_vdgs.extract_fragment_smiles import load_frags_dict

def _library_dirs(lib):
    return sorted(d for d in os.listdir(lib)
                  if os.path.isdir(os.path.join(lib, d)) and not d.startswith('.'))

def _read_aliases(lib):
    """{alias: (representative, kind)} from fragment_aliases.tsv, or {} if absent."""
    path = os.path.join(lib, 'fragment_aliases.tsv')
    if not os.path.isfile(path):
        return {}, path
    out = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                out[parts[0]] = (parts[1], parts[2])
    return out, path

def check_aliases_resolve(lib, dict_keys, aliases):
    """Every directory is attributable to a dict key or a promoted representative."""
    reps = {rep for rep, kind in aliases.values()}
    by_dir = {utils.smiles_to_filename(k): k for k in dict_keys}
    by_dir.update({utils.smiles_to_filename(r): r for r in reps})
    orphans = [d for d in _library_dirs(lib) if d not in by_dir]
    if orphans:
        return False, (f'{len(orphans)} library directory(ies) match neither a fragment-dict '
                       f'key nor an alias-table representative, e.g. {orphans[:3]}')
    return True, f'all {len(_library_dirs(lib))} directories attributable'

def check_charged_forward(lib, dict_keys, aliases):
    """Every representative named in the alias table exists as a directory.

    SCOPE, and the bug this replaced: the first version demanded a directory for
    every CHARGED KEY IN THE DICT. The dict holds 106,926 keys; a library built at
    --min-support 250 holds 2,639, so ~8,800 charged keys legitimately have no
    directory and the check reported them all as failures. "In the dict" is not
    "in the library". The alias table is written per-library, so it is the right
    scope: if a row names a representative, that representative was built.

    This is the direction that prevents "I searched for the drawn charged form and
    found nothing" -- holding [N+](=O)[O-], the table sends you to the neutral
    directory that actually exists.
    """
    reps = {rep for rep, _kind in aliases.values()}
    charged_aliases = [a for a in aliases if '+' in a or '-' in a]
    if not reps:
        return False, 'alias table empty -- this check would be vacuous'
    if not charged_aliases:
        return False, ('no CHARGED aliases in the table, so nothing exercises the '
                       'charge-stripping path -- vacuous, investigate')
    missing = [r for r in reps
               if not os.path.isdir(os.path.join(lib, utils.smiles_to_filename(r)))]
    if missing:
        return False, (f'{len(missing)} of {len(reps)} alias representative(s) have no '
                       f'directory, e.g. {missing[:3]}')
    promoted = {rep for rep, kind in aliases.values() if kind == 'promoted'}
    return True, (f'all {len(reps)} representatives resolve '
                  f'({len(charged_aliases)} charged aliases, {len(promoted)} promoted)')

def check_completeness(lib):
    """Every fragment directory carries the completion marker."""
    dirs = _library_dirs(lib)
    partial = []
    for name in dirs:
        log = os.path.join(lib, name, f'{name}_log')
        if not os.path.isfile(log):
            partial.append(name)
        elif 'Job completed.' not in open(log, errors='replace').read():
            partial.append(name)
    if partial:
        return False, (f'{len(partial)} of {len(dirs)} fragment(s) incomplete (no '
                       f"'Job completed.'), e.g. {partial[:3]}. If jobs are still in "
                       f'qstat this is expected -- rerun once the fleet drains.')
    return True, f'all {len(dirs)} fragments complete'

def check_h_class_fields(lib):
    """Delegate to the diagnostic the generation guide already mandates."""
    tool = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        'ligand_vdgs', 'tools', 'h_class_diagnostic.py')
    if not os.path.isfile(tool):
        return False, f'h_class_diagnostic.py not found at {tool}'
    # The flag is --lib, not --vdg-lib-dir; guessing it produced an argparse
    # error that this check reported as a missing-column failure.
    proc = subprocess.run([sys.executable, tool, '--lib', lib, '--dry-run'],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout).strip().splitlines()
        return False, f'h_class_diagnostic --dry-run failed: {tail[-1] if tail else "?"}'
    return True, 'H-class columns present'

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vdg-lib-dir', required=True)
    parser.add_argument('--frags-dict', default='resources/database_frags_dict.pkl')
    parser.add_argument('--skip', default='', help='comma-separated check names to skip')
    args = parser.parse_args()
    lib = os.path.expanduser(args.vdg_lib_dir)
    if not os.path.isdir(lib):
        raise SystemExit(f'[ERROR]: no such library: {lib}')
    frags, _support, _meta = load_frags_dict(args.frags_dict)
    dict_keys = set()
    for stem_map in frags.values():
        dict_keys.update(stem_map.keys())
    aliases, alias_path = _read_aliases(lib)
    if not aliases:
        print(f'[WARNING] no alias table at {alias_path}; alias checks will be weak.')
    skip = {s.strip() for s in args.skip.split(',') if s.strip()}
    checks = [('aliases-resolve', lambda: check_aliases_resolve(lib, dict_keys, aliases)),
              ('charged-forward', lambda: check_charged_forward(lib, dict_keys, aliases)),
              ('completeness', lambda: check_completeness(lib)),
              ('h-class-fields', lambda: check_h_class_fields(lib))]
    print(f'Post-build checks on {lib}')
    print(f'  fragment-dict keys: {len(dict_keys)}   alias rows: {len(aliases)}')
    failed = []
    for name, fn in checks:
        if name in skip:
            print(f'  {name:<18} SKIPPED')
            continue
        ok, detail = fn()
        print(f'  {name:<18} {"OK  " if ok else "FAIL"}  {detail}')
        if not ok:
            failed.append(name)
    if failed:
        print(f'\nFAILED: {", ".join(failed)}')
        return 1
    print('\nAll post-build checks passed.')
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
