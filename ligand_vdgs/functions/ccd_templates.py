"""The CCD as the source of ligand chemistry, keyed by atom name.

`components.cif` gives bond orders, aromaticity, formal charge and hydrogen
count per atom name. That is what the pipeline needs and what neither RDKit nor
OpenBabel can supply reliably from coordinates alone (docs/pitfalls.md), so the
template is authoritative and perception is only the fallback.
"""
import functools
import gzip
import os
import pickle
import sqlite3
from collections import namedtuple

# Bumped when the parsed fields or the blob layout change; the reader refuses a
# database built by a different parser rather than reading fields that moved.
TEMPLATE_PARSER_VERSION = 2

DEFAULT_CCD_DIR = '/wynton/group/degradolab/skt/docking/databases/ccd'
CCD_DIR_ENV = 'LIGAND_VDGS_CCD_DIR'
TEMPLATE_DB_NAME = 'templates.sqlite'
COMPONENTS_NAME = 'components.cif.gz'

# `leaving` matters: OXT/HXT and friends are absent from the PDB by design when
# the ligand is covalently attached, so their absence is not a name-mapping
# failure.
TemplateAtom = namedtuple('TemplateAtom',
                          'name alt_name element charge aromatic leaving xyz')
TemplateBond = namedtuple('TemplateBond', 'a b order aromatic')
Template = namedtuple('Template', 'comp_id atoms bonds')

_BOND_ORDER = {'SING': 1, 'DOUB': 2, 'TRIP': 3, 'QUAD': 4}


class TemplateDbMismatch(Exception):
    """The template database was built by a different parser version."""


def ccd_dir():
    return os.environ.get(CCD_DIR_ENV) or DEFAULT_CCD_DIR


def template_db_path():
    return os.path.join(ccd_dir(), TEMPLATE_DB_NAME)


def components_path():
    return os.path.join(ccd_dir(), COMPONENTS_NAME)


def split_cif_row(line):
    """Split an mmCIF loop row into values, honouring quoted fields.

    `line.split()` is wrong and fails silently: CCD `alt_atom_id` values are
    written `"N A"` for chlorophyll and heme nitrogens, so a naive split
    produces one field too many, the row is rejected as malformed, and those
    atoms vanish from the template -- which sent every CLA and HEC ligand to
    OpenBabel perception (4.5% of ligands in a 200-structure census).

    A quote only opens a value at the start of a token and only closes it when
    followed by whitespace or end of line. That second half matters: nucleotide
    atom names are written `O5'` unquoted, and a naive quote scanner swallows
    the rest of the row from the apostrophe onwards.
    """
    out, i, n = [], 0, len(line)
    while i < n:
        while i < n and line[i].isspace():
            i += 1
        if i >= n:
            break
        quote = line[i]
        if quote in "'\"":
            j, closed = i + 1, False
            while True:
                k = line.find(quote, j)
                if k == -1:
                    break
                if k + 1 >= n or line[k + 1].isspace():
                    out.append(line[i + 1:k])
                    i, closed = k + 1, True
                    break
                j = k + 1
            if closed:
                continue
        j = i
        while j < n and not line[j].isspace():
            j += 1
        out.append(line[i:j])
        i = j
    return out


def _unquote(value):
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in '"\'':
        value = value[1:-1]
    return value


def iter_templates(stream):
    """Parse `components.cif` text into Template objects, one per component.

    Two syntaxes, both of which occur: a `loop_` for multi-row categories, and
    plain `_category.key value` items when a category has exactly one row --
    which is how every monatomic ion (ZN, NA, ...) writes its single atom.
    Handling only the loop form drops all of them.

    Columns are located by loop header, never by fixed index: the header order
    is the file's contract and the atom loop alone has 21 columns, most of them
    coordinates this does not read.
    """
    comp_id = None
    atoms, bonds = [], []
    header, loop = None, None
    single = {}

    def take_single():
        """Flush a non-loop `_chem_comp_atom.*` / `_chem_comp_bond.*` item."""
        nonlocal single
        if not single:
            return
        if 'atom_id' in single:
            atoms.append(_atom_from_row(single))
        elif 'atom_id_1' in single:
            bonds.append(_bond_from_row(single))
        single = {}

    for line in stream:
        if line.startswith('data_'):
            take_single()
            if comp_id is not None:
                yield Template(comp_id, atoms[:], bonds[:])
            comp_id = line[len('data_'):].strip()
            atoms, bonds, header, loop, single = [], [], None, None, {}
            continue
        if line.startswith('loop_'):
            take_single()
            header, loop = [], None
            continue
        if line.startswith('_chem_comp_atom.') or line.startswith('_chem_comp_bond.'):
            body = line.strip().split('.', 1)[1]
            kind = 'atom' if line.startswith('_chem_comp_atom.') else 'bond'
            parts = body.split(None, 1)
            if len(parts) == 1:
                # Header line of a loop: key with no value on the same line.
                if header is None:
                    header = []
                header.append(parts[0])
                loop = kind
            else:
                # Single-row item; the value is on the line.
                if loop is not None and header:
                    take_single()
                    header, loop = None, None
                single[parts[0]] = _unquote(parts[1])
            continue
        if line.startswith('_') or line.startswith('#'):
            take_single()
            if not line.startswith('_'):
                loop = None
            continue
        if loop is None or not line.strip():
            continue
        parts = split_cif_row(line)
        if len(parts) != len(header):
            # A quoted value containing a space, or a multi-line entry. These
            # occur in descriptor loops, not the two read here; skip rather
            # than mis-assign a column.
            continue
        row = dict(zip(header, parts))
        atoms.append(_atom_from_row(row)) if loop == 'atom' else \
            bonds.append(_bond_from_row(row))
    take_single()
    if comp_id is not None:
        yield Template(comp_id, atoms[:], bonds[:])


def _atom_from_row(row):
    charge = _unquote(row.get('charge', '0'))
    return TemplateAtom(
        name=_unquote(row['atom_id']),
        alt_name=_unquote(row.get('alt_atom_id', row['atom_id'])),
        element=_unquote(row['type_symbol']).capitalize(),
        charge=int(charge) if charge.lstrip('-+').isdigit() else 0,
        aromatic=_unquote(row.get('pdbx_aromatic_flag', 'N')) == 'Y',
        leaving=_unquote(row.get('pdbx_leaving_atom_flag', 'N')) == 'Y',
        # Idealised coordinates, used to build test fixtures and never for
        # science: observed geometry always comes from the structure.
        xyz=_ideal_xyz(row))


def _ideal_xyz(row):
    try:
        return tuple(float(row[f'pdbx_model_Cartn_{axis}_ideal'])
                     for axis in 'xyz')
    except (KeyError, ValueError):
        return (0.0, 0.0, 0.0)


def _bond_from_row(row):
    return TemplateBond(
        a=_unquote(row['atom_id_1']), b=_unquote(row['atom_id_2']),
        order=_BOND_ORDER.get(_unquote(row.get('value_order', 'SING')), 1),
        aromatic=_unquote(row.get('pdbx_aromatic_flag', 'N')) == 'Y')


def build_template_db(components_cif_gz, db_path):
    """Parse components.cif.gz into a comp_id -> template sqlite. Returns count."""
    tmp = f'{db_path}.{os.getpid()}.tmp'
    if os.path.exists(tmp):
        os.remove(tmp)
    conn = sqlite3.connect(tmp)
    ok = False
    try:
        conn.execute('CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)')
        conn.execute('CREATE TABLE templates (comp_id TEXT PRIMARY KEY, blob BLOB)')
        n = 0
        with gzip.open(components_cif_gz, 'rt') as handle:
            rows = ((t.comp_id, pickle.dumps((t.atoms, t.bonds), protocol=5))
                    for t in iter_templates(handle))
            for batch in _batched(rows, 5000):
                conn.executemany('INSERT OR REPLACE INTO templates VALUES (?, ?)',
                                 batch)
                n += len(batch)
        fetched_on = ''
        marker = os.path.join(os.path.dirname(components_cif_gz), 'FETCHED_ON')
        if os.path.isfile(marker):
            with open(marker) as handle:
                fetched_on = handle.read().strip()
        conn.executemany('INSERT OR REPLACE INTO meta VALUES (?, ?)', [
            ('parser_version', str(TEMPLATE_PARSER_VERSION)),
            ('source', os.path.abspath(components_cif_gz)),
            ('fetched_on', fetched_on)])
        conn.commit()
        ok = True
    finally:
        conn.close()
        if not ok and os.path.exists(tmp):
            # A crash part-way would otherwise leave a partial sqlite in the
            # CCD directory, where nothing cleans it up.
            os.remove(tmp)
    os.replace(tmp, db_path)
    return n


def _batched(iterable, size):
    batch = []
    for item in iterable:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


_CONNECTION = None


def connection():
    """One read-only connection per process, opened lazily."""
    global _CONNECTION
    if _CONNECTION is None:
        path = template_db_path()
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f'{path} not found. Build it once on a login node with '
                f'scripts/build_ccd_templates.py, or point {CCD_DIR_ENV} at a '
                'directory that has it.')
        conn = sqlite3.connect(f'file:{path}?mode=ro', uri=True)
        version = conn.execute(
            "SELECT value FROM meta WHERE key = 'parser_version'").fetchone()
        if version is None or int(version[0]) != TEMPLATE_PARSER_VERSION:
            conn.close()
            raise TemplateDbMismatch(
                f'{path} was built by parser version '
                f'{version[0] if version else "unknown"}, this code reads '
                f'{TEMPLATE_PARSER_VERSION}; rebuild it.')
        _CONNECTION = conn
    return _CONNECTION


@functools.lru_cache(maxsize=4096)
def get_template(comp_id):
    """Template for one CCD id, or None if the CCD does not have it."""
    row = connection().execute(
        'SELECT blob FROM templates WHERE comp_id = ?', (str(comp_id).strip(),)
    ).fetchone()
    if row is None:
        return None
    atoms, bonds = pickle.loads(row[0])
    return Template(str(comp_id).strip(), atoms, bonds)


def index_by_name(template):
    """{atom name: TemplateAtom}, including alt names that do not collide.

    Older depositions write the alt form (`1HB` for `HB1`, primes vs stars in
    nucleotides), so matching on `atom_id` alone turns a perfectly good template
    into a fallback.
    """
    by_name = {atom.name: atom for atom in template.atoms}
    for atom in template.atoms:
        if atom.alt_name and atom.alt_name not in by_name:
            by_name[atom.alt_name] = atom
    return by_name


def template_h_counts(template):
    """{heavy atom name: number of H bonded to it in the template}.

    The template's nominal protonation, which is the point: it does not depend
    on whether hydrogens were placed, modelled or stripped (rebuild-notes 6).

    Leaving hydrogens (HXT, and the second amide H of a free amino acid) are
    counted, so a peptide-linked N or a C-terminus whose OXT is absent reads
    one H high. Left as is deliberately: leaving hydrogens sit on heteroatoms,
    and only carbon's binary H flag is key material, so no key is affected --
    only `cg_num_h`, which is metadata and already provenance-bound.
    """
    elements = {atom.name: atom.element for atom in template.atoms}
    counts = {name: 0 for name, el in elements.items() if el not in ('H', 'D')}
    for bond in template.bonds:
        for near, far in ((bond.a, bond.b), (bond.b, bond.a)):
            if elements.get(far) in ('H', 'D') and near in counts:
                counts[near] += 1
    return counts
