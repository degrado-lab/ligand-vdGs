'''Identify a parent database by its contents rather than by its path.

A built library, and the cost estimate that decided what to build, are both tied to one
parent database. The path is not that database: the library is meant to be copyable
between machines, where the same parent database legitimately sits somewhere else, while a
different database can just as legitimately be written to the same path. So consumers
compare this identity and treat the path as informational.

The identity is a sha256 over the sorted ``(relative path, byte size)`` of every
structure in the parent db, plus the structure count. That is one directory walk with a
stat per file -- no file contents are read -- so it is cheap enough to compute on demand
and is cached in ``<pdb_dir>/DB_IDENTITY``.

Its limit, stated rather than implied: it detects added, removed, renamed and
resized files, not an edit that leaves every size unchanged. It is a
same-database check, not a tamper check.
'''
import hashlib
import json
import os
import time

from ligand_vdgs.functions import parent_db

IDENTITY_FILENAME = 'DB_IDENTITY'
IDENTITY_VERSION = 1


def compute_identity(pdb_dir):
    '''{'sha256', 'n_structures', 'version'} from a walk of *pdb_dir*.

    Only the parent DB's structures count, so the DB_IDENTITY file itself, logs and
    manifests sitting alongside them do not change the identity.
    '''
    entries = []
    for stem, path in parent_db.iter_structures(pdb_dir):
        entries.append((os.path.relpath(path, pdb_dir), os.path.getsize(path)))
    entries.sort()
    digest = hashlib.sha256()
    digest.update(f'{IDENTITY_VERSION}\n{len(entries)}\n'.encode())
    for rel, size in entries:
        digest.update(f'{rel}\t{size}\n'.encode())
    return {'version': IDENTITY_VERSION, 'sha256': digest.hexdigest(),
            'n_structures': len(entries)}


def identity_path(pdb_dir):
    return os.path.join(pdb_dir, IDENTITY_FILENAME)


def read_identity(pdb_dir):
    '''The cached record, or None if absent or unreadable.

    Unreadable is treated as absent on purpose: a truncated or hand-edited file must
    fall back to recomputing, not abort a build that was otherwise fine.
    '''
    path = identity_path(pdb_dir)
    if not os.path.isfile(path):
        return None
    try:
        with open(path) as handle:
            record = json.load(handle)
    except (ValueError, OSError):
        return None
    if record.get('version') != IDENTITY_VERSION or not record.get('sha256'):
        return None
    return record


def write_identity(pdb_dir, source_dir=None, record=None):
    '''Cache *record* (default: computed now) to <pdb_dir>/DB_IDENTITY. Returns it.

    *source_dir* is recorded for provenance only and is deliberately not part of the
    hash -- a copy of a database is the same database, so its identity must not depend
    on where it was produced.
    '''
    record = dict(record or compute_identity(pdb_dir))
    record['written'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    if source_dir:
        record['repaired_from'] = os.path.abspath(source_dir)
    final = identity_path(pdb_dir)
    tmp = final + '.tmp'
    # tmp + rename, like the library provenance: a crash mid-dump would otherwise leave
    # JSON that read_identity has to discard, silently costing a full walk every run.
    with open(tmp, 'w') as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write('\n')
    os.replace(tmp, final)
    return record


def identity_of(pdb_dir, cache=True):
    '''The cached identity if there is one, else compute it and try to cache it.

    Caching is best effort: a read-only db (a shared database another group owns)
    must still be usable, so a failed write is not an error.
    '''
    record = read_identity(pdb_dir)
    if record is not None:
        return record
    record = compute_identity(pdb_dir)
    if cache:
        try:
            record = write_identity(pdb_dir, record=record)
        except OSError:
            pass
    return record
