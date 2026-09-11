'''Build resources/ccd_polymer_types.tsv (CCD id, type, parent) from the wwPDB CCD.

s01 uses the type column to tell a chain-bonded modified amino acid (`L-PEPTIDE
LINKING`, e.g. KCX/SEP/LLP: becomes an ATOM record, neither slot nor ligand) from a
covalently attached cofactor (`NON-POLYMER`, e.g. PLP/HEM: stays HETATM and is mined
as a ligand). See _prep_filters.modified_residues_to_protein. Compute nodes have no
network, so the table is fetched here once and committed; it covers every CCD entry
(~45k rows) so a new database needs no refetch unless it uses entries newer than the
table. Streams components.cif.gz (~120 MB); a few minutes on a login node.

    python scripts/fetch_ccd_polymer_types.py [--out resources/ccd_polymer_types.tsv]
'''
import argparse
import gzip
import os
import sys
import urllib.request

CCD_URL = 'https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz'


def _unquote(value):
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in '"\'':
        return value[1:-1]
    return value


def iter_chem_comp(stream):
    '''Yield (id, type, parent) per data block of a components.cif text stream.'''
    cur = None
    for raw in stream:
        line = raw.rstrip('\n')
        if line.startswith('data_'):
            if cur is not None:
                yield cur['id'], cur.get('type', ''), cur.get('parent', '')
            cur = {'id': line[5:].strip()}
        elif cur is None:
            continue
        elif line.startswith('_chem_comp.type'):
            cur['type'] = _unquote(line[len('_chem_comp.type'):])
        elif line.startswith('_chem_comp.mon_nstd_parent_comp_id'):
            cur['parent'] = _unquote(line[len('_chem_comp.mon_nstd_parent_comp_id'):])
    if cur is not None:
        yield cur['id'], cur.get('type', ''), cur.get('parent', '')


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--out', default=os.path.join(os.path.dirname(__file__), '..',
                                                 'resources', 'ccd_polymer_types.tsv'))
    p.add_argument('--url', default=CCD_URL)
    args = p.parse_args()

    tmp = args.out + '.tmp'
    n = 0
    with urllib.request.urlopen(args.url, timeout=120) as resp, \
            gzip.open(resp, 'rt', encoding='utf-8', errors='replace') as text, \
            open(tmp, 'w') as out:
        out.write(f'# CCD id\ttype\tmon_nstd_parent_comp_id  (source: {args.url})\n')
        for cid, ctype, parent in iter_chem_comp(text):
            out.write(f'{cid}\t{ctype}\t{parent}\n')
            n += 1
    if n < 30000:
        os.remove(tmp)
        sys.exit(f'only {n} entries parsed; the CCD has ~45k, so the download or the '
                 'parser is broken. Nothing written.')
    os.replace(tmp, args.out)
    print(f'wrote {n} CCD entries to {args.out}')


if __name__ == '__main__':
    main()
