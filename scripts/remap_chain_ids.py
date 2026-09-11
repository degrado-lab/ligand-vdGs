'''Apply s01's text-level repairs to an already-trimmed parent database.

s01 applies them at trim time (`_prep_filters.text_repairs`: chain-ID remap, then
MSE/SEC renames, HETATM amino acids to ATOM, chain-bonded peptide-linking residues to
ATOM). Databases trimmed before those fixes -- `prepwizard_BioLiP2` -- need them applied
after the fact, text-level only: no re-trim, nothing parsed before the fixes.

Two modes:

    # survey: what would change, write nothing (one line per affected structure)
    python scripts/remap_chain_ids.py --pdb-dir <db> --list > affected.txt

    # repair the whole database into a new mirror-layout dir, the pdb_dir every
    # later step points at. Unchanged structures are copied verbatim, so the output
    # is self-contained.
    python scripts/remap_chain_ids.py --pdb-dir <db> --out-dir <repaired> --all \
        --manifest <repaired>.manifest.tsv

    # once every shard has finished: stamp the finished parent database with its identity
    python scripts/remap_chain_ids.py --pdb-dir <db> --out-dir <repaired> --write-identity

`--shard i/--num-shards n` splits the deterministic stem order for an SGE array job;
each task writes its own manifest. Without `--all` only structures with
multi-character chain IDs are written (the original behaviour, for a spot repair).

Structures with more distinct chains than single characters (ChainIdOverflow) are
skipped and not written, as in s01 -- a chain collision makes them unmineable anyway.

Safe to re-run: every output is written to a pid-tagged .tmp and renamed into place, so
an interrupted task leaves no half-written PDB and a repeated shard simply overwrites.
Each structure's outcome is independent -- a failed read or a failed write costs that one
structure a manifest row of `skipped_error`, never the rest of the shard.
'''
import argparse
import os
import shutil
import sys

from ligand_vdgs.functions import parent_db
from ligand_vdgs.functions.db_identity import IDENTITY_FILENAME, write_identity
from ligand_vdgs.preprocessing._chain_ids import ChainIdOverflow, remap_remarks
from ligand_vdgs.preprocessing._prep_filters import load_ccd_polymer_types, text_repairs

MANIFEST_HEADER = ('stem\tchanged\tchain_remap\trenamed\tamino_acids_to_atom\t'
                   'modres_to_atom\tunknown_resnames\n')
# `changed` values. A skipped structure gets a row too, so the manifest accounts for
# every structure in the source database and the output count can be reconciled against
# it without re-reading the job logs.
SKIPPED_OVERFLOW = 'skipped_overflow'
SKIPPED_ERROR = 'skipped_error'


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--pdb-dir', required=True, help='Parent database (mirror layout).')
    p.add_argument('--out-dir', help='Where repaired copies go (mirror layout).')
    p.add_argument('--all', action='store_true',
                   help='Repair every structure and copy the unchanged ones, so --out-dir '
                        'is a complete database. Without it only structures with '
                        'multi-character chain IDs are written.')
    p.add_argument('--manifest', help='TSV of per-structure outcomes (one row per '
                                      'structure processed).')
    p.add_argument('--list', action='store_true',
                   help='Print affected stems and what would change; write nothing.')
    p.add_argument('--write-identity', action='store_true',
                   help=f'Repair nothing: stamp --out-dir with its {IDENTITY_FILENAME} '
                        f'and exit. Run this once after every shard has finished, since '
                        f'the identity covers the whole parent database.')
    p.add_argument('--shard', type=int, default=0,
                   help='0-based index of this shard of the stem order.')
    p.add_argument('--num-shards', type=int, default=1)
    p.add_argument('stems', nargs='*', help='Restrict to these stems (default: all).')
    args = p.parse_args(argv)
    if not args.list and not args.out_dir:
        p.error('--out-dir is required unless --list')
    if args.write_identity and args.list:
        p.error('--write-identity and --list do nothing together')
    if not 0 <= args.shard < args.num_shards:
        p.error(f'--shard must be in [0, {args.num_shards})')
    return args


def repair_lines(lines, ccd_types):
    '''(new_lines, mapping, stats, changed). *changed* is byte-level, not derived from
    the stats: a rewrite the stats do not account for must still be written out.'''
    new_lines, mapping, stats = text_repairs(lines, ccd_types)
    return new_lines, mapping, stats, new_lines != lines


def _manifest_row(stem, changed, mapping, stats):
    remap = ','.join(f'{k}->{v}' for k, v in mapping.items()) or '-'
    unknown = ','.join(stats['unknown_resnames']) or '-'
    changed = changed if isinstance(changed, str) else int(changed)
    return (f"{stem}\t{changed}\t{remap}\t{stats['renamed']}\t"
            f"{stats['amino_acids_to_atom']}\t{stats['modres_to_atom']}\t{unknown}\n")


def _skipped_row(stem, kind, reason):
    return _manifest_row(stem, kind, {}, {'renamed': 0, 'amino_acids_to_atom': 0,
                                          'modres_to_atom': 0,
                                          'unknown_resnames': [reason]})


def _ensure_dir(path):
    '''os.makedirs(exist_ok=True) is not atomic on BeeGFS.

    CPython implements exist_ok by catching the mkdir error and re-checking
    os.path.isdir; with 100 array tasks writing into the same two-letter shard
    directories, a stale negative metadata cache right after another node created the
    directory makes that re-check fail and FileExistsError escapes. It killed tasks 7,
    27 and 92 of the first full run. The directory does exist in that case, so swallowing
    the error is correct; if it genuinely does not, the write below fails loudly.
    '''
    try:
        os.makedirs(path, exist_ok=True)
    except FileExistsError:
        pass


def _tmp_path(out):
    '''Scratch name for an atomic write. Pid-tagged so a rerun overlapping a live shard
    cannot collide, and suffixed .tmp so a killed task's leftovers are invisible to the
    `*.pdb` glob every reader (and the database identity) uses.'''
    return f'{out}.{os.getpid()}.tmp'


def _write_repaired(out, remarks, lines):
    '''Write then rename, so *out* is either absent or complete. A task killed mid-write
    would otherwise leave a truncated PDB that looks like a valid structure, and reruns
    overwrite cleanly instead of tripping over a partial file.'''
    tmp = _tmp_path(out)
    with open(tmp, 'w') as handle:
        handle.writelines(remarks)
        handle.writelines(lines)
    os.replace(tmp, out)


def _copy_structure(src, out):
    '''Same guarantee for the unchanged structures, which are the bulk of the database.'''
    tmp = _tmp_path(out)
    shutil.copyfile(src, tmp)
    os.replace(tmp, out)


def main(argv=None):
    args = parse_args(argv)
    if args.out_dir and os.path.abspath(args.out_dir) == os.path.abspath(args.pdb_dir):
        sys.exit('--out-dir must differ from --pdb-dir (never repair in place)')
    if args.write_identity:
        # After the array job, not during it: a shard sees only its own slice, so an
        # identity stamped mid-run would describe a partial database.
        record = write_identity(args.out_dir, source_dir=args.pdb_dir)
        print(f"[INFO] {os.path.join(args.out_dir, IDENTITY_FILENAME)}: "
              f"{record['n_structures']} structures, sha256 {record['sha256']}")
        return

    ccd_types = load_ccd_polymer_types()
    wanted = set(args.stems)

    # Deterministic order so a shard is reproducible and shards partition the database.
    entries = sorted(parent_db.iter_structures(args.pdb_dir))
    entries = [e for i, e in enumerate(entries) if i % args.num_shards == args.shard]
    if wanted:
        entries = [e for e in entries if e[0] in wanted]

    manifest = open(args.manifest, 'w') if args.manifest else None
    if manifest:
        manifest.write(MANIFEST_HEADER)

    n_seen = n_changed = n_copied = n_written = n_overflow = n_failed = 0
    unknown_all = set()
    for stem, path in entries:
        n_seen += 1
        try:
            with open(path) as f:
                lines = f.readlines()
            new_lines, mapping, stats, changed = repair_lines(lines, ccd_types)
        except ChainIdOverflow as e:
            n_overflow += 1
            print(f'[ERROR] {stem}: {e}', file=sys.stderr)
            if manifest:
                manifest.write(_skipped_row(stem, SKIPPED_OVERFLOW, str(e)))
            continue
        except Exception as e:  # one bad file must not abort a 65k-structure pass
            n_failed += 1
            print(f'[ERROR] {stem}: {e!r}', file=sys.stderr)
            if manifest:
                manifest.write(_skipped_row(stem, SKIPPED_ERROR, repr(e)))
            continue
        unknown_all.update(stats['unknown_resnames'])
        n_changed += changed

        # The row is written only once the structure's fate is settled, so a manifest
        # row always means "this structure was accounted for", never "this structure
        # would have been written if the write had not then failed".
        row = _manifest_row(stem, changed, mapping, stats)
        if args.list:
            if changed:
                print(row, end='')
        elif changed or args.all:
            out = parent_db.structure_path(args.out_dir, stem)
            try:
                _ensure_dir(os.path.dirname(out))
                if changed:
                    # The chain remap is recoverable only from the file itself.
                    _write_repaired(out, remap_remarks(mapping), new_lines)
                    n_written += 1
                else:
                    _copy_structure(path, out)
                    n_copied += 1
            except Exception as e:
                # A failed write is one structure's problem, not the shard's: the read
                # path above already had this guarantee, the write path did not, and a
                # single FileExistsError took out three whole tasks because of it.
                n_failed += 1
                print(f'[ERROR] {stem}: writing {out}: {e!r}', file=sys.stderr)
                row = _skipped_row(stem, SKIPPED_ERROR, repr(e))
        if manifest:
            manifest.write(row)

    if manifest:
        manifest.close()
    if unknown_all:
        print(f'[WARNING] {len(unknown_all)} HETATM resname(s) not in the CCD type table '
              f'were left as written; refresh with scripts/fetch_ccd_polymer_types.py: '
              f'{sorted(unknown_all)}', file=sys.stderr)
    print(f'[INFO] shard {args.shard}/{args.num_shards}: {n_seen} structures processed, '
          f'{n_changed} changed, {n_written} repaired, {n_copied} copied unchanged, '
          f'{n_overflow} skipped (too many chains), {n_failed} failed',
          file=sys.stderr)


if __name__ == '__main__':
    main()
