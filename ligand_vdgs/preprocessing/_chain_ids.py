'''Force single-character chain IDs on parent PDBs before anything parses them.

The PDB fixed-column format gives the chain ID one column (22). Every reader in
this pipeline honours that: ProDy (``line[21]``), probe, vdG-miner's line hashing
(``line[12:26]``) and its ligand scan (``line[21]``). Some sources (BioLiP
receptor files of large assemblies) write a two-character chain ID across
columns 21-22 instead, e.g. ``HETATM 1440  N   ASPA4  66`` for chain ``A4``.
Column 21 is then dropped by every reader, so chain ``A4`` is read as chain ``4``
and its residue 66 merges with residue 66 of the real chain ``4``. ProDy >= 2.x
gives the merged atoms a single resindex, so the collision survives
parsePDB/writePDB with no trace in the output text. In the 2026-09 database this
affected 99 structures, all of which vdG-miner then dropped whole (probe and PDB
lines no longer hash to the same key).

The remap is deterministic per file: single-character chains keep their ID,
each multi-character chain gets the first unused character from ``CHAIN_POOL``
in order of first appearance. The mapping is written back as REMARK lines so it
can be recovered from the file itself.
'''
import string

CHAIN_POOL = string.ascii_uppercase + string.digits + string.ascii_lowercase
CHAIN_REMAP_REMARK = 'REMARK 900 CHAIN ID REMAPPED'
_CHAIN_RECORDS = ('ATOM  ', 'HETATM', 'ANISOU', 'TER   ')


class ChainIdOverflow(ValueError):
    '''More distinct chains than single characters available to name them.'''


def chain_id(line):
    '''The chain ID a writer put in columns 21-22 (what readers see is column 22 only).'''
    return line[20:22].strip()


def _has_chain_columns(line):
    return line[:6] in _CHAIN_RECORDS and len(line) > 21


def multichar_chain_ids(lines):
    '''Distinct chain IDs longer than one character, in order of first appearance.'''
    seen = []
    for line in lines:
        if _has_chain_columns(line):
            cid = chain_id(line)
            if len(cid) > 1 and cid not in seen:
                seen.append(cid)
    return seen


def chain_id_mapping(lines):
    '''{multichar_id: single_char} for *lines*; empty if nothing needs remapping.'''
    used = set()
    for line in lines:
        if _has_chain_columns(line):
            used.add(chain_id(line))
    mapping = {}
    pool = iter(c for c in CHAIN_POOL if c not in used)
    for cid in multichar_chain_ids(lines):
        try:
            mapping[cid] = next(pool)
        except StopIteration:
            raise ChainIdOverflow(
                f'{len(used)} distinct chain IDs; only {len(CHAIN_POOL)} single '
                f'characters are available') from None
    return mapping


def remap_chain_ids(lines):
    '''Rewrite multi-character chain IDs so column 21 is blank on every record.

    Returns (new_lines, mapping). Lines are returned unchanged (same object) when
    nothing needs remapping. Raises ChainIdOverflow if the file has more distinct
    chains than CHAIN_POOL can name; the caller should skip the structure rather
    than write it with a collision.
    '''
    mapping = chain_id_mapping(lines)
    if not mapping:
        return lines, mapping
    out = []
    for line in lines:
        if _has_chain_columns(line):
            new = mapping.get(chain_id(line))
            if new is not None:
                line = line[:20] + ' ' + new + line[22:]
        out.append(line)
    return out, mapping


def remap_remarks(mapping):
    '''REMARK lines recording *mapping*, to prepend to the written PDB.'''
    return [f'{CHAIN_REMAP_REMARK} {old} -> {new}\n' for old, new in mapping.items()]


def parse_remap_remarks(lines):
    '''Inverse of remap_remarks: {old: new} from a file's REMARK lines.'''
    mapping = {}
    for line in lines:
        if line.startswith(CHAIN_REMAP_REMARK):
            old, _, new = line[len(CHAIN_REMAP_REMARK):].split()
            mapping[old] = new
    return mapping


def assert_single_char_chains(lines, name=''):
    '''Raise ValueError if any coordinate record still has a non-blank column 21.'''
    bad = multichar_chain_ids(lines)
    if bad:
        raise ValueError(f'{name or "PDB"} has multi-character chain IDs {bad}; '
                         f'run remap_chain_ids first')
