# write_vdg_hit_pdbs.py

import os
import sys
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

from ligand_vdgs.functions import vdg_npz_utils as vdg_npz
from ligand_vdgs.functions import utils
from ligand_vdgs.functions.vdg_pdb_io import (
    NameRegistry, format_residue_tag, fresh_dir, sanitize, strip_known_exts,
    write_pdb_gz)


MAX_PER_FRAG_BSR = 10


def parse_bsr_combo(bsr_combo):
    out = []
    for token in str(bsr_combo).split(";"):
        seg, chain, resnum = token.split(":")
        out.append((str(seg), str(chain), int(resnum)))
    return out


def make_R_t(row):
    R = np.array([[row[f"R{i}{j}"] for j in range(3)] for i in range(3)], dtype=float)
    t = np.array([row[f"t{i}"] for i in range(3)], dtype=float)
    return R, t


def slot_order_from_perm_idx(bucket_parts, aa_perm_idx):
    """
    Library vdM slot indices in query BSR position order.

    `aa_perm_idx` indexes `vdg_npz_utils.aa_perm_indices(bucket_parts)`; the hit
    finder applies the chosen permutation as `vdg_bb[perm]`, so perm[i] is the
    library slot superposed onto BSR position i. The permutation group is
    inverse-closed, so the index is meaningless without that direction: reading
    it the other way silently transposes the residue correspondence for every
    bucket with repeated labels (bb_bb, ASP_ASP, ...).

    The hit finder always writes a valid index, so anything else means the TSV
    and the library disagree; raise rather than fall back to stored order, which
    would mislabel the correspondence silently.
    """
    idx = int(aa_perm_idx)
    perms = vdg_npz.aa_perm_indices(list(bucket_parts))
    if not 0 <= idx < len(perms):
        raise ValueError(
            f"aa_perm_idx {idx} out of range for bucket {list(bucket_parts)} "
            f"({len(perms)} permutations)")
    return list(perms[idx])


def load_centroid(npz_path, idx, slot_order):
    """build_kwargs, parent name, and vdM residue tags in query BSR order."""
    with np.load(npz_path) as data:
        build_kwargs = vdg_npz.nr_build_kwargs(data, idx)
    source_name = strip_known_exts(build_kwargs["parent_pdb_path"])

    # scrr_* are stored in library slot order; slot_order (from aa_perm_idx)
    # reorders them into query BSR position order, so tag i names the library
    # residue matched to residue i of the enclosing bsr folder. Without it, any
    # bucket with repeated labels mislabels that correspondence.
    #
    # Not deduplicated: one tag per vdM slot keeps the count equal to the subset size.
    # Two slots of one vdG are always different residues anyway. Unlike a materialized
    # library vdG, a hit name carries no ligand tag -- the CG here is the query's own,
    # and the query is already named at the front of the file name. Note the tags are
    # NOT the last fields of a hit name: make_hit_stem appends the RMSD after them, so
    # a consumer strips one trailing field before taking 4 * subset_size from the right.
    tags = [
        format_residue_tag(*fields)
        for fields in zip(build_kwargs["scrr_resname"], build_kwargs["scrr_seg"],
                          build_kwargs["scrr_chain"], build_kwargs["scrr_resnum"])
    ]
    return build_kwargs, source_name, [tags[j] for j in slot_order]


def make_bsr_folder(bsr_combo, aa_bucket):
    """Directory naming the query residues a hit was matched against.

    ``aa_bucket`` supplies the residue *labels* (library side) and ``bsr_combo`` the
    query positions; they are zipped, so they must be the same length. zip would
    otherwise truncate to the shorter and emit a folder naming only some of the
    residues -- two different BSR combos could then collide into one directory, and
    nothing downstream would notice, since fresh_dir short-circuits a directory it
    has already claimed.
    """
    aas = str(aa_bucket).split("_")
    bsr = parse_bsr_combo(bsr_combo)
    if len(aas) != len(bsr):
        raise ValueError(
            f"aa_bucket {aa_bucket!r} has {len(aas)} label(s) but bsr_combo "
            f"{bsr_combo!r} has {len(bsr)} residue(s); they must describe the same "
            "vdG. Check that the hits TSV and the vdG library are from the same run.")
    return "_".join(
        format_residue_tag(aa, seg, chain, resnum)
        for aa, (seg, chain, resnum) in zip(aas, bsr)
    )


def make_hit_stem(row, source_name, scrr_tags):
    """Output basename stem (no extension, not yet uniquified) for one hit.

    Ends ``..._<vdM tags>_<rmsd>``. The trailing RMSD is one field (it contains a
    ``.``, not a ``_``), so parsing the tags means dropping it first.

    The stem names the parent PDB and the library residues but not which CG site
    of that ligand the vdG came from, so two hits from the same deposition collide
    whenever their RMSDs also tie at 2 dp. Nothing checks the path on write, so the
    caller's NameRegistry is the only thing keeping the second hit from clobbering
    the first.
    """
    query = sanitize(strip_known_exts(row["pdbfile"]))
    frag = utils.smiles_to_filename(row["frag"])
    return (f"{query}_{frag}_{sanitize(source_name)}_{'_'.join(scrr_tags)}"
            f"_{float(row['vdg_rmsd']):.2f}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits-tsv", required=True)
    parser.add_argument("--vdg-lib-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--max-per-frag-bsr", type=int, default=MAX_PER_FRAG_BSR)
    args = parser.parse_args()

    df = pd.read_csv(args.hits_tsv, sep="\t")
    counts = defaultdict(int)
    made_dirs, names_by_dir = set(), defaultdict(NameRegistry)
    skipped = []

    def skip(row_idx, row, detail):
        msg = (f"Row {row_idx}: {row['pdbfile']} | {row['frag']} | "
               f"{row['bsr_combo']}\n    {detail}")
        skipped.append(msg)
        print(f"[SKIP] {msg}")

    for row_idx, (_, row) in enumerate(df.iterrows()):
        key = (row["pdbfile"], row["frag"], row["bsr_combo"])
        if counts[key] >= args.max_per_frag_bsr:
            continue
        counts[key] += 1

        npz_path = vdg_npz.vdg_npz_path(
            args.vdg_lib_dir, row["frag"], int(row["subset_size"]), row["aa_bucket"])
        vdg_idx = int(row["vdg_index"])

        if not os.path.exists(npz_path):
            skip(row_idx, row, f"NPZ file not found: {npz_path}")
            continue

        with np.load(npz_path) as data:
            num_vdgs = len(data["nr_cg_coords"])
            bucket_parts = [str(x) for x in data["aa_bucket_parts"]]
        if not 0 <= vdg_idx < num_vdgs:
            skip(row_idx, row, f"vdg_index {vdg_idx} out of range "
                               f"[0, {num_vdgs-1}] in {npz_path}")
            continue

        slot_order = slot_order_from_perm_idx(bucket_parts, row["aa_perm_idx"])
        build_kwargs, source_name, scrr_tags = load_centroid(
            npz_path, vdg_idx, slot_order)
        ag, _ = vdg_npz.build_vdg_atomgroup_from_npz(**build_kwargs)
        ag = vdg_npz.apply_rigid_transform(ag.copy(), *make_R_t(row))

        out_dir = os.path.join(
            args.outdir,
            sanitize(strip_known_exts(row["pdbfile"])),
            'vdg_matches',
            utils.smiles_to_filename(row["frag"]),
            make_bsr_folder(row["bsr_combo"], row["aa_bucket"]),
        )
        # No --overwrite here: one run writes many <query>/<frag>/<bsr> leaves and
        # they are only discovered as rows stream past, so a blanket 'replace what's
        # there' would delete an unpredictable subset of an earlier run. Clear the
        # --outdir by hand instead.
        fresh_dir(out_dir, made_dirs)
        out_path = names_by_dir[out_dir].claim(
            make_hit_stem(row, source_name, scrr_tags), out_dir)
        write_pdb_gz(ag, out_path)

    if skipped:
        print(f"\n[ERROR] {len(skipped)} row(s) were skipped:")
        for msg in skipped:
            print(f"  {msg}")
        print("\nPossible causes:")
        print("  1. The TSV and NPZ files are from different runs (regenerate one or both)")
        print("  2. The NPZ file is corrupted or truncated")
        print("  3. There's a bug in vdg_hit_finder.py that generated invalid indices")
        sys.exit(1)


if __name__ == "__main__":
    main()
