#!/usr/bin/env python
"""Build the CCD template sqlite that ligand_perception reads.

Run once on a login node (compute nodes have no network, but this only reads a
local file). ~30 s for the 51k-component CCD.

    python scripts/build_ccd_templates.py

Writes `templates.sqlite` beside `components.cif.gz` in the CCD directory, so
the two share one FETCHED_ON provenance. Neither belongs in git.
"""
import argparse
import os
import sys

from ligand_vdgs.functions import ccd_templates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-c', '--components',
                        default=ccd_templates.components_path(),
                        help='path to components.cif.gz')
    parser.add_argument('-o', '--out', default=ccd_templates.template_db_path(),
                        help='output sqlite path')
    args = parser.parse_args()

    if not os.path.isfile(args.components):
        sys.exit(f'{args.components} not found; fetch it on a login node first.')
    n = ccd_templates.build_template_db(args.components, args.out)
    print(f'Wrote {n:,} templates to {args.out}')


if __name__ == '__main__':
    main()
