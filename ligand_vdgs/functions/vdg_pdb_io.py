# vdg_pdb_io.py

"""Naming and file-writing helpers shared by the vdG PDB writers.

Used by ``generate_vdgs/materialize_vdg_pdbs.py`` and ``tools/write_vdg_hit_pdbs.py``.
Both produce directories loaded into PyMOL wholesale, so both need the same answer to
"what is this file called and is that name already taken".

Name grammar: fields are joined with ``_``, and a residue tag is exactly four of them,
``seg_chain_resnum_resname`` (empty segment kept, so the usual tag reads ``_A_370_ASP``).
The fixed field count lets a consumer split a basename on ``_`` and take residue tags as
fixed-size groups **from the right** -- the head isn't positionally parseable, since a
sanitized SMILES and an AA bucket label each contribute a variable number of fields.
``NameRegistry``'s duplicate suffix uses ``~`` for the same reason.

Output-collision policy: never overwrite, never skip. Each leaf output directory must be
absent or empty when the run first reaches it (``fresh_dir``), and names within a run are
made unique by a ``NameRegistry`` -- scoped per run if basenames must stay unique across
the directories one PyMOL session loads together.
"""

import gzip
import os
import shutil
import tempfile

import prody as pr

from ligand_vdgs.functions import parent_db
from ligand_vdgs.functions.utils import set_up_outdir


DUP_SUFFIX = "~"


def sanitize(x):
    return str(x).replace("/", "_").replace("\\", "_").replace(" ", "_")


def strip_known_exts(path):
    return parent_db.stem_of(path)


def normalize_seg(seg):
    """Segment as a string, mapping every spelling of "absent" (``None``, ``""``,
    ``"None"``, depending on whether it came from ProDy, an npz ``U*`` column or a
    parsed tag) to ``""``."""
    return "" if seg is None or str(seg) in ("", "None") else str(seg)


def _tag_field(value):
    """One residue-tag field, guaranteed ``_``-free.

    A literal ``_`` (segment names may carry one) and the ``/``, ``\\``, space that
    ``sanitize`` would turn into ``_`` all become ``-``; either would add a field and
    break the fixed count. Whitespace is stripped, not substituted, so a blank-padded
    segment reads back as the empty field rather than a run of dashes.
    """
    out = str(value).strip()
    for ch in ("_", "/", "\\", " "):
        out = out.replace(ch, "-")
    return out


def format_residue_tag(resname, seg, chain, resnum):
    """``seg_chain_resnum_resname``, always exactly four ``_``-separated fields.

    The empty segment field is kept deliberately: dropping it would make the tag 3
    fields sometimes and 4 others, defeating right-anchored parsing.
    """
    return "_".join((_tag_field(normalize_seg(seg)), _tag_field(chain),
                     str(int(resnum)), _tag_field(resname)))


class NameRegistry:
    """Hands out output paths, numbering every file of a repeated stem from ~1.

    A collision is only discovered when the second file arrives, so occurrence 1 is
    renamed from ``stem`` to ``stem~1`` at that moment: a bare stem always means the
    stem occurred once. Write each file immediately after claiming it -- the rename
    follows the recorded path, and a caller that claims twice before writing gets a
    FileNotFoundError instead of a silently un-renumbered occurrence 1.

    The suffix avoids ``_`` (it would read back as an extra field). ``~`` is a
    SMARTS bond primitive, but fragment keys come from MolToSmiles output, which
    never emits one, and nothing else in a name does either.
    """

    def __init__(self, ext=".pdb.gz"):
        self._ext = ext
        self._counts = {}       # stem -> times claimed
        self._first_path = {}   # stem -> path of occurrence 1, until it is renamed

    def _path(self, out_dir, stem, n=None):
        suffix = "" if n is None else f"{DUP_SUFFIX}{n}"
        return os.path.join(out_dir, f"{stem}{suffix}{self._ext}")

    def claim(self, stem, out_dir):
        """Path to write ``stem``'s next file to, renaming occurrence 1 if needed."""
        n = self._counts.get(stem, 0) + 1
        self._counts[stem] = n
        if n == 1:
            path = self._path(out_dir, stem)
            self._first_path[stem] = path
            return path
        if n == 2:
            self._renumber_first(stem)
        return self._path(out_dir, stem, n)

    def _renumber_first(self, stem):
        """Move occurrence 1 from ``stem`` to ``stem~1`` now that a second exists."""
        first = self._first_path.pop(stem)
        target = self._path(os.path.dirname(first), stem, 1)
        if os.path.exists(target):
            # Leaves are empty per fresh_dir, so this is not a name this run
            # assigned -- only a stem that itself contains "~" (a hand-written
            # CG SMARTS) reaches here. Refuse rather than clobber.
            raise FileExistsError(
                f"Cannot renumber {first} to {target}: target already exists.")
        if not os.path.exists(first):
            raise FileNotFoundError(
                f"Cannot renumber {first} to {target}: it was claimed but never "
                "written. Write each file immediately after claiming its name.")
        os.replace(first, target)


def fresh_dir(path, made_dirs):
    """Claim ``path`` as an output directory this run may write into.

    Creates it, requiring it to be empty on first claim so leftovers from an earlier run
    can't be mistaken for this run's output (no --overwrite, no resume; clear by hand).
    ``made_dirs`` records claims so the emptiness check runs once per directory.
    """
    path = os.path.abspath(path)
    if path in made_dirs:
        return path
    set_up_outdir(path)
    made_dirs.add(path)
    return path


def write_pdb_gz(ag, out_path):
    """Write an AtomGroup as gzipped PDB. Assumes ``fresh_dir`` prepared the parent.

    Written to a temp name in the destination directory and renamed into place; a kill
    mid-copy would otherwise leave an undetected truncated .pdb.gz on the final name.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_pdb = os.path.join(tmpdir, "tmp.pdb")
        pr.writePDB(tmp_pdb, ag)
        # Same filesystem as out_path, so os.replace is atomic.
        tmp_gz = f"{out_path}.{os.getpid()}.tmp"
        try:
            with open(tmp_pdb, "rb") as fin, gzip.open(tmp_gz, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            os.replace(tmp_gz, out_path)
        except BaseException:
            try:
                os.remove(tmp_gz)
            except OSError:
                pass
            raise
