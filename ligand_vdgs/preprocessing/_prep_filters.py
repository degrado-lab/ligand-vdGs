'''
Filters applied to parent PDBs before Schrodinger prepwizard (s02) sees them.

Prepwizard silently mangles residues it cannot interpret as amino acids, and the damage
is invisible downstream because the output is still a well-formed PDB. Everything here
exists to remove the *input* that provokes that, so it must run before s02. It is called
from both s01 (which is already rewriting PDBs) and s02 itself, because s01 is optional.
'''

import os

import prody as pr

AAs = ('ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
       'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR')

# A vdM's local frame is built from these three, so a residue missing any of them is
# unusable downstream regardless of what prepwizard does to it.
BACKBONE = ('N', 'CA', 'C')


def incomplete_backbone_resindices(parsed):
    '''Resindices of amino acid residues missing any of N, CA or C.

    These are the residues that provoke prepwizard. A residue modelled with its backbone
    N alone (density ran out mid-chain) cannot be built as an amino acid, so prepwizard
    reclassifies the orphan N as a free molecule, protonates it to ammonium and re-emits
    it carrying the resname/chain/resnum of a *ligand* -- giving that ligand a stray atom
    tens of A away (2y1x SAH A:1001, whose second atom named N is really THR D:478's).
    A residue reduced to a lone CA becomes methane the same way (1c53 GLY A:3).

    Dropping them costs nothing: a vdM needs all three backbone atoms for its frame.
    Measured on 4,000 prepared PDBs, this is 159 of 472,225 protein residues (0.034%).

    Single-atom HET groups are handled separately by single_carbon_het_resindices.
    Neither covers HET residues that prepwizard *renames* to a standard amino acid
    (CYT -> CYS in 5buv); that damage happens to residues this filter has no reason to
    remove, so it cannot be prevented here.
    '''
    names, resindices, resnames = (parsed.getNames(), parsed.getResindices(),
                                   parsed.getResnames())
    is_aa = {i for i, rn in zip(resindices, resnames) if rn in AAs}
    complete = set(resindices[names == BACKBONE[0]])
    for bb_name in BACKBONE[1:]:
        complete &= set(resindices[names == bb_name])
    return is_aa - complete


def single_carbon_het_resindices(parsed):
    '''Resindices of non-amino-acid residues consisting of one carbon heavy atom.

    Prepwizard mislabels these the same way, absorbing the atom into a nearby ligand
    residue tens of A away (CF0 in 5d61/2h6m/2h9h/2hal, 0QE in 6bfl/6bgr/4qu0/4qu9/
    4qub/4qul -- 11 structures in a 65,600-PDB database). They are real deposited
    density: the terminal methyl of a covalent inhibitor, LINKed to a catalytic
    cysteine. Dropping them loses nothing mineable, because the smallest CG is four
    connected heavy atoms, so a lone carbon can never be part of one.

    Restricted to carbon on purpose. Monatomic ions -- metals, halides, LU in 1plu --
    are single-atom residues too and must survive.
    '''
    names, resindices, resnames = (parsed.getNames(), parsed.getResindices(),
                                   parsed.getResnames())
    elements = parsed.getElements()
    per_res = {}
    for i, rn, nm, el in zip(resindices, resnames, names, elements):
        el = (el or '').strip().upper() or nm.strip('0123456789')[:1].upper()
        if el in ('H', 'D'):
            continue
        per_res.setdefault(i, (rn, []))[1].append(el)
    return {i for i, (rn, els) in per_res.items()
            if rn not in AAs and rn != 'HOH' and els == ['C']}


def resindices_to_drop(parsed):
    '''Every resindex this module would remove before prepwizard.'''
    return incomplete_backbone_resindices(parsed) | single_carbon_het_resindices(parsed)


def drop_prepwizard_hazard_residues(parsed, resindices=None):
    '''Filter *resindices* (default: all of them) down to the residues worth keeping.'''
    if resindices is None:
        resindices = set(parsed.getResindices())
    return sorted(set(resindices) - resindices_to_drop(parsed))


def clean_pdb_file(in_path, out_path):
    '''Write *in_path* to *out_path* without the residues prepwizard would mangle.

    Returns the number of residues dropped, or None if the file could not be parsed (in
    which case nothing is written and the caller should fall back to the original).
    '''
    parsed = pr.parsePDB(in_path)
    if parsed is None:
        return None
    keep = drop_prepwizard_hazard_residues(parsed)
    n_dropped = len(set(parsed.getResindices())) - len(keep)
    if not n_dropped:
        return 0
    pr.writePDB(out_path, parsed.select(f"resindex {' '.join(str(i) for i in keep)}"))
    return n_dropped


# --- Protein chains written as HETATM ---------------------------------------------


def embedded_amino_acid_hetatm_to_atom(lines):
    '''Rewrite HETATM records of standard amino acids that are peptide-bonded to another
    amino acid as ATOM records. Returns (lines, number of residues rewritten).

    The 2026-09 database carried 670 structures whose two-character-chain segments
    (see _chain_ids) were written as HETATM throughout: 113,720 protein residues, and
    every HETATM standard amino acid in the database sat on such a chain. vdG-miner's
    ligand scan reads HETATM only, so once the chain IDs are repaired those residues
    would be mined as ligands (TRP indole as an aromatic CG, ASP as a carboxylate...).
    A free amino acid -- HETATM, bonded to nothing, carboxylate complete (OXT) -- is a
    real ligand and is left alone. Rewritten: residues with a peptide bond (N or C within
    COVALENT_BOND_CUTOFF of another amino acid's C or N), and unbonded residues with no
    OXT, which are chain fragments whose neighbours fell outside the s01 sphere (8any:
    ARG 120 alone on chain AF). A Cys adduct or a single-residue ligand keeps its records.
    '''
    import numpy as np
    key_of = lambda line: (line[21], line[22:27])
    bb = {}      # reskey -> {'N': xyz, 'C': xyz} over every standard amino acid
    het = set()  # reskeys of HETATM standard amino acids
    has_oxt = set()
    for line in lines:
        if line[:6] not in ('ATOM  ', 'HETATM') or line[17:20] not in AAs:
            continue
        name = line[12:16].strip()
        if name in ('N', 'C'):
            bb.setdefault(key_of(line), {})[name] = (
                float(line[30:38]), float(line[38:46]), float(line[46:54]))
        elif name == 'OXT':
            has_oxt.add(key_of(line))
        if line.startswith('HETATM'):
            het.add(key_of(line))
    if not het:
        return lines, 0

    keys = list(bb)
    n_xyz = np.array([bb[k].get('N', (9e9,) * 3) for k in keys])
    c_xyz = np.array([bb[k].get('C', (9e9,) * 3) for k in keys])
    # residue i's N bonded to residue j's C (i != j), or vice versa
    d = np.sqrt(((n_xyz[:, None, :] - c_xyz[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(d, np.inf)
    bonded = (d <= COVALENT_BOND_CUTOFF).any(axis=1) | (d <= COVALENT_BOND_CUTOFF).any(axis=0)
    convert = {k for k, b in zip(keys, bonded) if k in het and (b or k not in has_oxt)}
    if not convert:
        return lines, 0
    out = [('ATOM  ' + line[6:]) if line.startswith('HETATM') and key_of(line) in convert
           else line for line in lines]
    return out, len(convert)


# --- Restoring ligands prepwizard renames into protein -----------------------------



def _coord_key(line):
    '''Rounded xyz from a PDB line. Prepwizard runs with -noimpref and does not move
    heavy atoms, so this matches a residue across s02 without trusting chain/resnum.'''
    return line[30:54]


COVALENT_BOND_CUTOFF = 1.8  # A; longer than any C-C/C-N/C-S bond, shorter than contact


def _atoms(path):
    '''Per-atom (reskey, resname, is_het, name, is_heavy, xyz) from a PDB file.'''
    out = []
    with open(path) as fh:
        for line in fh:
            if line[:6] not in ('ATOM  ', 'HETATM'):
                continue
            out.append(((line[21], line[22:27]), line[17:20], line.startswith('HETATM'),
                        line[12:16].strip(), line[76:78].strip().upper() not in ('H', 'D'),
                        (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                        line[30:54]))
    return out


def snapshot_restorable_ligands(path):
    '''Heavy-atom coordinate key -> original resname, for HET residues that are real
    ligands rather than part of a polymer chain.

    Prepwizard rewrites some HET residues into ATOM records under a standard amino acid
    resname. find_cg_matches reads only HETATM, so an affected ligand vanishes from
    mining with no warning (CYT -> CYS in 5buv/5epu, HSE -> HIS in 6a0s, SRO -> SER in
    7bs2). This snapshot is what lets s02 put them back.

    "Real ligand" = written as HETATM, not named for an amino acid, and **not covalently
    bonded to a standard amino acid**. That last test is the discriminator, and it is
    geometric on purpose: SEQRES would say the same thing but does not survive s01
    (ProDy's writePDB drops it), and an atom-name backbone test misses chain-embedded
    chromophores whose backbone is named N1/CA1/C1 (CR8 in 3tmr, MDO in 2o7d).

    Chain-embedded residues are therefore left as prepwizard wrote them: they are
    protein, the X slot label already covers them, and turning them back into HETATM
    would make them mineable as ligands -- a change in library composition, not a bug
    fix. Covalent-inhibitor warheads (CF0, 0QE) test as bonded for the same reason and
    are dropped upstream anyway.
    '''
    import numpy as np
    atoms = _atoms(path)
    if not atoms:
        return {}
    aa_xyz = np.array([a[5] for a in atoms if a[1].strip() in AAs and a[4]] or [[9e9] * 3])

    by_res = {}
    for a in atoms:
        by_res.setdefault(a[0], []).append(a)

    out = {}
    for reskey, res_atoms in by_res.items():
        resname, is_het = res_atoms[0][1], res_atoms[0][2]
        if not is_het or resname.strip() in AAs:
            continue
        heavy = [a for a in res_atoms if a[4]]
        if not heavy:
            continue
        xyz = np.array([a[5] for a in heavy])
        if np.sqrt(((xyz[:, None, :] - aa_xyz[None, :, :]) ** 2).sum(-1)).min() \
                <= COVALENT_BOND_CUTOFF:
            continue  # bonded into the chain -- protein, not a ligand
        for a in heavy:
            out[a[6]] = resname
    return out


def restore_renamed_ligands(prepped_path, snapshot):
    '''Undo prepwizard's rename in place. Returns the number of residues restored.

    Only the resname and the ATOM/HETATM record type are put back. Renamed *atoms*
    (6a0s HSE: N -> NA) are left as prepwizard wrote them -- CG matching and the later
    per-name selection both read this same file, so internal consistency is all that is
    needed.
    '''
    if not snapshot:
        return 0
    with open(prepped_path) as fh:
        lines = fh.readlines()

    votes, heavy = {}, {}
    for line in lines:
        if line[:6] not in ('ATOM  ', 'HETATM') or line[76:78].strip().upper() in ('H', 'D'):
            continue
        key = (line[21], line[22:27])
        heavy[key] = heavy.get(key, 0) + 1
        orig = snapshot.get(_coord_key(line))
        if orig is not None:
            votes.setdefault(key, {})
            votes[key][orig] = votes[key].get(orig, 0) + 1

    restore = {}
    for key, counts in votes.items():
        resname, n = max(counts.items(), key=lambda kv: kv[1])
        # Majority of the residue's heavy atoms must come from that one original ligand,
        # which is what stops an absorbed orphan atom from renaming its host.
        if n * 2 >= heavy[key]:
            restore[key] = resname

    changed = set()
    for i, line in enumerate(lines):
        if line[:6] not in ('ATOM  ', 'HETATM'):
            continue
        key = (line[21], line[22:27])
        resname = restore.get(key)
        if resname is None or resname == line[17:20]:
            continue
        lines[i] = 'HETATM' + line[6:17] + resname + line[20:]
        changed.add(key)
    if changed:
        with open(prepped_path, 'w') as fh:
            fh.writelines(lines)
    return len(changed)


# --- Modified residues without prepwizard --------------------------------------------


# Isosteric selenium variants: same atom names apart from SE for SD/SG, so the residue
# is a full slot once renamed (vdg_struct_utils._ALTERNATE_HEAVY_ATOMS admits SE).
MODIFIED_RESIDUE_RENAMES = {'MSE': 'MET', 'SEC': 'CYS'}

CCD_POLYMER_TYPES_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..',
                                      'resources', 'ccd_polymer_types.tsv')


def load_ccd_polymer_types(path=CCD_POLYMER_TYPES_PATH):
    '''{CCD id: _chem_comp.type} from scripts/fetch_ccd_polymer_types.py's table.
    Missing table -> {} with a warning: every non-standard HETATM residue then counts
    as unknown and is left alone (prepwizard-era behaviour minus the renames).'''
    if not os.path.isfile(path):
        print(f'[WARNING] {path} not found; run scripts/fetch_ccd_polymer_types.py. '
              'Chain-bonded modified residues will be left as HETATM (mined as ligands).')
        return {}
    types = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            types[fields[0]] = fields[1] if len(fields) > 1 else ''
    return types


def is_peptide_linking(ccd_type):
    '''`L-PEPTIDE LINKING`, `D-PEPTIDE LINKING`, `PEPTIDE LINKING` (any case); not
    `peptide-like` (a ligand class) and not the gamma/beta-peptide variants.'''
    return ccd_type.strip().upper().endswith('PEPTIDE LINKING')


def modified_residues_to_protein(lines, ccd_types):
    '''Make chain-embedded modified residues protein records, so that without
    prepwizard's renaming they are neither dropped as slots nor scanned as ligands.

    Prepwizard used to rename every modified residue to its parent (LLP -> LYS,
    KCX -> LYS, SEP -> SER...), which vdG-miner then accepted as a slot (labelled X
    when the modification is within 4.5 A). Without it they stay HETATM under their
    own resname: the 20-resname slot whitelist rejects them and the HETATM-only ligand
    scan would mine KCX/SEP/PTR/CSO as carbamate/phosphate/sulfenic-acid ligand CGs.
    Audit of the 2026-09 database: 3322 such residues, 1686 within 4.5 A of a ligand;
    MSE alone 4881 residues, 1148 near a ligand.

    Three rewrites, in order:
    1. MODIFIED_RESIDUE_RENAMES (MSE -> MET, SEC -> CYS): isosteres whose slots are
       worth keeping; record type unchanged, atom names unchanged (SE stays SE).
    2. embedded_amino_acid_hetatm_to_atom: standard amino acids written as HETATM
       (including the ones just renamed) become ATOM when peptide-bonded or when an
       unbonded chain fragment; a free amino-acid ligand with OXT stays HETATM.
    3. Every other HETATM residue whose CCD type is `*PEPTIDE LINKING` and which has
       a heavy atom within COVALENT_BOND_CUTOFF of a standard amino acid's heavy atom
       becomes an ATOM record under its own resname: not a slot, not a ligand. The
       CCD type is the discriminator, not backbone atom names -- in the audited
       database 353 chain-bonded cofactor residues (SAM, SAH, NXL, 0G6) carry N/CA/C
       names and must stay ligands, while 23 peptide-linking residues do not carry
       them. A peptide-linking residue bonded only to other non-standard residues (a
       peptide ligand made of ncAAs) stays HETATM; one bonded to a standard residue
       of a peptide ligand is lost as a ligand (tiny numbers: 85 chain-bonded
       peptide-linking residues in 65,600 structures, 15 near a ligand, 10 of them
       SEC). Unknown resnames (not in *ccd_types*) are left alone and reported.

    Returns (lines, stats) with stats = {'renamed': n residues, 'amino_acids_to_atom':
    n, 'modres_to_atom': n, 'unknown_resnames': sorted list}.
    '''
    import numpy as np
    coord_records = ('ATOM  ', 'HETATM')
    key_of = lambda line: (line[21], line[22:27])

    # 1. isosteric renames (also on ANISOU so the file stays self-consistent)
    renamed = set()
    out = []
    for line in lines:
        if line[:6] in coord_records + ('ANISOU',):
            new = MODIFIED_RESIDUE_RENAMES.get(line[17:20])
            if new is not None:
                if line[:6] in coord_records:
                    renamed.add(key_of(line))
                line = line[:17] + new + line[20:]
        out.append(line)
    lines = out

    # 2. HETATM standard amino acids (incl. the renamed ones) -> ATOM
    lines, n_aa = embedded_amino_acid_hetatm_to_atom(lines)

    # 3. chain-bonded peptide-linking non-standard residues -> ATOM
    aa_xyz, cand, unknown = [], {}, set()
    for line in lines:
        if line[:6] not in coord_records or line[76:78].strip().upper() in ('H', 'D'):
            continue
        resname = line[17:20]
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if resname in AAs:
            aa_xyz.append(xyz)
        elif line.startswith('HETATM') and resname != 'HOH':
            ccd_type = ccd_types.get(resname.strip())
            if ccd_type is None:
                unknown.add(resname.strip())
            elif is_peptide_linking(ccd_type):
                cand.setdefault(key_of(line), []).append(xyz)
    convert = set()
    if cand and aa_xyz:
        aa_xyz = np.asarray(aa_xyz)
        for key, xyz in cand.items():
            d2 = ((np.asarray(xyz)[:, None, :] - aa_xyz[None, :, :]) ** 2).sum(-1)
            if d2.min() <= COVALENT_BOND_CUTOFF ** 2:
                convert.add(key)
    if convert:
        lines = [('ATOM  ' + line[6:]) if line.startswith('HETATM') and key_of(line) in convert
                 else line for line in lines]
    return lines, {'renamed': len(renamed), 'amino_acids_to_atom': n_aa,
                   'modres_to_atom': len(convert), 'unknown_resnames': sorted(unknown)}
