"""Parent-structure identity: stems, entry IDs, and the mirror layout.

The parent database is an RCSB-style mirror, ``<dir>/<shard>/<stem>.pdb`` with
``shard = stem[1:3].lower()``. A stem is ``1abc`` or ``1abc_2`` (assembly 2);
its *entry* is the text before the first underscore, which also holds for
PLINDER system IDs (``1abc__1__1.A__1.B``). Nothing else in the pipeline slices
a stem or builds one of these paths.
"""
import glob
import os

STRUCTURE_EXT = ".pdb"
KNOWN_STRUCT_EXTS = (".pdb.gz", ".pdb", ".cif.gz", ".cif", ".gz")


def shard(stem):
    return stem[1:3].lower()


def structure_path(pdb_dir, stem):
    return os.path.join(pdb_dir, shard(stem), stem + STRUCTURE_EXT)


def stem_of(path):
    """Basename with a structure-file extension removed (``x/1abc.pdb.gz`` -> ``1abc``)."""
    base = os.path.basename(str(path))
    for ext in KNOWN_STRUCT_EXTS:
        if base.endswith(ext):
            return base[:-len(ext)]
    return os.path.splitext(base)[0]


def entry_of(stem):
    """The deposition accession shared by every assembly/system of one entry."""
    return str(stem).split("_", 1)[0]


def iter_structures(pdb_dir):
    """Sorted ``(stem, path)`` for every structure in the mirror."""
    pattern = os.path.join(pdb_dir, "*", "*" + STRUCTURE_EXT)
    for path in sorted(glob.glob(pattern)):
        yield stem_of(path), path


def is_mirror(pdb_dir):
    """Whether at least one structure sits in the shard its stem prescribes."""
    if not os.path.isdir(pdb_dir):
        return False
    for stem, path in iter_structures(pdb_dir):
        if os.path.basename(os.path.dirname(path)) == shard(stem):
            return True
    return False
