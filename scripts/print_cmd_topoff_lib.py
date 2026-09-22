#!/usr/bin/env python
"""top-off diff: build_set = select_fragments(N_new) \\
already_built, already_built via Frags.check_vdg_job_status, never directory
existence. Prints the diff, confirms, then prints (never runs) the --mode
include-only command to submit it. Refuses when --n-new >= --n-old. --help."""
import argparse
import os

from ligand_vdgs.generate_vdgs.extract_fragment_smiles import (
    load_frags_dict, prepare_fragments, select_fragments)
from ligand_vdgs.generate_vdgs.make_sge_scripts_for_frags import partition_by_completion

def compute_build_set(frags_dict, support_pooled, max_size, n_old, n_new, vdg_lib_dir):
    if n_new >= n_old:
        raise SystemExit(
            f"[ERROR]: --n-new ({n_new}) >= --n-old ({n_old}).")
    build_set, already_built, partial = partition_by_completion(
        select_fragments(frags_dict, n_new, max_size, support=support_pooled,
                        prepared=prepare_fragments(frags_dict, max_size)),
        vdg_lib_dir)
    if set(build_set) & set(already_built):
        raise SystemExit(
            "[ERROR]: fragment already marked 'Job completed.'.")
    return build_set, already_built, partial

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--frags-dict', required=True)
    ap.add_argument('--vdg-lib-dir', required=True)
    ap.add_argument('--max-size', type=int, required=True)
    ap.add_argument('--n-old', type=int, required=True,
                     help='The --min-support the live library was built at.')
    ap.add_argument('--n-new', type=int, required=True,
                     help='The lower --min-support to top off to.')
    ap.add_argument('--yes', action='store_true', help='Skip the interactive confirmation.')
    args = ap.parse_args()

    frags_dict, support_pooled, _frags_meta = load_frags_dict(args.frags_dict)
    build_set, already_built, partial = compute_build_set(
        frags_dict, support_pooled, args.max_size, args.n_old, args.n_new, args.vdg_lib_dir)

    print(f'top-off: --min-support {args.n_old} -> {args.n_new}')
    print(f'  {len(already_built)} fragment(s) at the new threshold are already built '
          f"('Job completed.' in their log) -- excluded.")
    print(f'  {len(partial)} fragment(s) have a directory but did NOT complete -- '
          f'these ARE in build_set (an h_rt kill, not a finished job).')
    print(f'  build_set: {len(build_set)} fragment(s) to submit.')
    if not build_set:
        print('Nothing to submit.')
        return

    if not args.yes and input(f'Submit {len(build_set)} job(s)? [y/N] ').strip().lower() != 'y':
        print('Aborted -- nothing submitted, nothing printed to run.')
        return

    print('\nRun this to actually build build_set (goes through --include-only\'s '
          'own provenance check, which this script deliberately does not '
          'duplicate):\n')
    quoted_frags = ' '.join("'" + s + "'" for s in build_set)
    print(f"INCLUDE={quoted_frags} \\\n"
          f"  {os.path.dirname(os.path.dirname(os.path.abspath(__file__)))}"
          f"/run_production_frags.sh --mode include-only")

if __name__ == '__main__':
    main()
