'''
Shared helpers for identify_bioisosteres: AA/AA-pair enrichment and
similarity-matrix clustering. Used by compute_aa_profiles.py,
compare_aa_profiles.py, compare_cg_geometries.py.

Contact categories partition by (residue, moiety), not just residue: <AA>,
<AA>-bb, and pooled bb. prior = bg_freq * size / sum(bg_freq * size), which
sums to exactly 1 in every backbone mode. `X` (non-canonical contacts) is
never a category and is reported separately; see category_priors/
load_bucket_counts docstrings for the per-category bg_freq/size table.
'''
import math
import os
import numpy as np

from ligand_vdgs.functions.Frags import check_vdg_job_status
from ligand_vdgs.functions.vdg_npz_utils import CHARGE_SIGNS
from ligand_vdgs.functions.vdg_struct_utils import BB_LABEL, BB_LABELS, NONCANONICAL_AA_LABEL
from ligand_vdgs.tools.reference_values import AA_background_freq, AA_sidechain_size, BB_HEAVY_ATOMS

# Suffix marking a per-residue backbone category, e.g. 'GLY-bb'. Deliberately not
# '_': bucket file names join labels with '_' and callers split pair labels on it.
BB_CATEGORY_SUFFIX = '-' + BB_LABEL

# How backbone-mediated contacts enter the category set.
BB_MODES = ('per-residue', 'pooled', 'off')


def is_bb_category(label):
    '''True for the pooled backbone category and for any <AA>-bb category.'''
    return label in BB_LABELS or label.endswith(BB_CATEGORY_SUFFIX)


def category_table(bb_mode='per-residue'):
    '''Return category -> (background frequency, heavy-atom size).'''
    if bb_mode not in BB_MODES:
        raise ValueError(f'Unknown bb_mode: {bb_mode!r}; expected one of {BB_MODES}')

    table = {aa: (AA_background_freq[aa], AA_sidechain_size[aa]) for aa in AA_background_freq
             if AA_sidechain_size.get(aa, 0) > 0}
    if bb_mode == 'pooled':
        table[BB_LABEL] = (1.0, BB_HEAVY_ATOMS)
    elif bb_mode == 'per-residue':
        table.update({aa + BB_CATEGORY_SUFFIX: (freq, BB_HEAVY_ATOMS)
                      for aa, freq in AA_background_freq.items()})
    return table


def category_priors(bb_mode='per-residue', norm_method='bg_weighted'):
    '''Return category priors for one normalization method.'''
    table = category_table(bb_mode)
    if norm_method == 'none':
        if bb_mode != 'off':
            raise ValueError("norm_method='none' has no size term, so a backbone category "
                             "has no defensible prior against a sidechain one. Use "
                             "bb_mode='off' with it, or use 'bg_weighted'.")
        denom = sum(bg for bg, _ in table.values())
        return {label: bg / denom for label, (bg, _) in table.items()}

    if norm_method != 'bg_weighted':
        raise ValueError(f'Unknown norm_method: {norm_method!r}')
    denom = sum(bg * size for bg, size in table.values())
    return {label: bg * size / denom for label, (bg, size) in table.items()}


def load_bucket_counts(npz_dir, bb_mode='off', allow_incomplete=False):
    '''Return (support by category, excluded-contact counts).'''
    if bb_mode not in BB_MODES:
        raise ValueError(f'Unknown bb_mode: {bb_mode!r}; expected one of {BB_MODES}')

    npz_dir_abs = os.path.abspath(npz_dir)
    frag_dir = os.path.dirname(os.path.dirname(npz_dir_abs))
    cg_label = os.path.basename(frag_dir)
    vdg_lib_dir = os.path.dirname(frag_dir)
    if os.path.basename(os.path.dirname(npz_dir_abs)) != 'nr_vdgs':
        raise ValueError(
            f'Expected npz_dir of the form <vdg_lib_dir>/<cg_label>/nr_vdgs/'
            f'<subset_size>, got {npz_dir!r}. The completion check derives the '
            'fragment label from this shape and cannot verify an arbitrary '
            'directory; pass a bucket directory from the library.')
    if not check_vdg_job_status(cg_label, vdg_lib_dir):
        if not allow_incomplete:
            raise ValueError(
                f'vdG generation is incomplete for {cg_label!r} (no '
                f"'Job completed.' line in its log under {vdg_lib_dir!r}).")
        print(f'WARNING: counting partial buckets for {cg_label!r}; vdG '
              'generation has not completed.')

    counts, extras = {}, {'X': 0, 'noncanonical_bb': 0}
    sign_dirs = [d for d in (os.path.join(npz_dir, s) for s in CHARGE_SIGNS)
                 if os.path.isdir(d)]
    if not sign_dirs:
        raise ValueError(
            f'{npz_dir!r} contains none of the charge-sign subdirectories '
            f'{CHARGE_SIGNS} a writer produces.')
    for sign_dir in sign_dirs:
        for fname in sorted(os.listdir(sign_dir)):
            if not fname.endswith('.npz'):
                continue
            bucket = fname[:-4]
            parts = bucket.split('_')
            if NONCANONICAL_AA_LABEL in parts:
                with np.load(os.path.join(sign_dir, fname)) as npz:
                    extras['X'] += len(npz['cluster_id'])
                continue
            has_bb = any(p in BB_LABELS for p in parts)
            if has_bb and bb_mode == 'off':
                continue

            with np.load(os.path.join(sign_dir, fname)) as npz:
                num_parents = npz['cluster_num_parents']
                if not (has_bb and bb_mode == 'per-residue'):
                    counts[bucket] = counts.get(bucket, 0) + int(num_parents.sum())
                    continue
                slot_parts = [str(x) for x in npz['aa_bucket_parts']]
                if sorted(slot_parts) != sorted(parts):
                    raise ValueError(f'{fname}: aa_bucket_parts {slot_parts} and the file name '
                                     f'{parts} are not the same labels, so which column holds '
                                     'which residue slot is undetermined.')
                for row, resname_row in enumerate(npz['nr_scrr_resname']):
                    labels = [str(resname_row[i]) + BB_CATEGORY_SUFFIX if p in BB_LABELS else p
                              for i, p in enumerate(slot_parts)]
                    if any(lab[:-len(BB_CATEGORY_SUFFIX)] not in AA_background_freq
                           for lab, p in zip(labels, slot_parts) if p in BB_LABELS):
                        extras['noncanonical_bb'] += 1
                        continue
                    key = '_'.join(sorted(labels))
                    counts[key] = counts.get(key, 0) + int(num_parents[row])
    return counts, extras


def _log_enrichment(count, prior, total):
    '''log(observed/expected) for one category, or None if not computable.'''
    if prior is None or prior <= 0 or count <= 0 or total <= 0:
        return None
    expected = total * prior
    return math.log(count / expected) if expected > 0 else None


def calc_single_aa_propensities(nr_vdgs_dir, norm_method='bg_weighted', counts=None,
                                bb_mode='per-residue'):
    '''Returns (propensities, counts) where propensities maps a contact category to
    its log-enrichment, sorted ascending, and counts is category -> cluster count.

    Pass `counts` from load_bucket_counts with the *same* bb_mode; a mismatch leaves
    categories in the denominator that have no prior in the numerator.'''
    if counts is None:
        single_dir = os.path.join(nr_vdgs_dir, '1')
        if not os.path.isdir(single_dir):
            raise FileNotFoundError(f'No single-AA vdGs found at: {single_dir}')
        counts, _ = load_bucket_counts(single_dir, bb_mode=bb_mode)
    total = sum(counts.values())
    priors = category_priors(bb_mode, norm_method)

    propensities = {label: val for label, count in counts.items()
                    if (val := _log_enrichment(count, priors.get(label), total)) is not None}
    return dict(sorted(propensities.items(), key=lambda item: item[1])), counts


def calc_aa_pair_propensities(nr_vdgs_dir, norm_method='bg_weighted', counts=None,
                              bb_mode='pooled'):
    '''Returns dict of "CAT1_CAT2" -> log-enrichment (unsorted).

    Pair prior is the independent-draw product p(cat1, cat2) = sf * p1 * p2,
    sf = 2 for a heterogeneous bucket (absorbs both orderings), 1 for
    homogeneous; sum(sf * p1 * p2) = (sum p)^2 = 1, so pair priors partition
    exactly like the single ones.'''
    if counts is None:
        pair_dir = os.path.join(nr_vdgs_dir, '2')
        if not os.path.isdir(pair_dir):
            raise FileNotFoundError(f'No AA-pair vdGs found at: {pair_dir}')
        counts, _ = load_bucket_counts(pair_dir, bb_mode=bb_mode)
    total = sum(counts.values())
    priors = category_priors(bb_mode, norm_method)

    propensities = {}
    for bucket, count in counts.items():
        parts = bucket.split('_')
        if len(parts) != 2:
            continue
        p1, p2 = priors.get(parts[0]), priors.get(parts[1])
        if p1 is None or p2 is None:
            continue
        val = _log_enrichment(count, (1 if parts[0] == parts[1] else 2) * p1 * p2, total)
        if val is not None:
            propensities[bucket] = val

    return propensities
