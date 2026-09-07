'''
Shared helpers for the identify_bioisosteres pipeline: AA/AA-pair enrichment
calculation and similarity-matrix clustering. Used by compute_aa_profiles.py,
compare_aa_profiles.py, compare_cg_geometries.py, and visualize_bioisosteres.py.

Contact categories
------------------
A vdG bucket label is a claim about *which moiety* of a residue contacts the CG,
not just which residue. So the null model partitions contacts by (residue,
moiety) rather than by residue:

    category            bg_freq                       size
    ------------------  ----------------------------  --------------------
    <AA>       (e.g. ASP)   PDB frequency of that AA      sidechain heavy atoms
    <AA>-bb    (e.g. GLY-bb) PDB frequency of that AA      4  (N, CA, C, O)
    bb         (pooled)      1.0  -- every residue has one 4

`bg_freq` reads as "fraction of residues presenting this category" and `size` as
"heavy atoms that category exposes". The prior is

    p(category) = bg_freq * size / sum(bg_freq * size)

and the denominator is the background-weighted mean residue size, so the priors
sum to **exactly 1** in every backbone mode -- including the pooled one, where
`sum_r freq(r)*n_sc(r) + 1.0*4 = sum_r freq(r)*AA_size(r)`.

Glycine has no sidechain category by construction; under per-residue backbone
deconvolution its share of the prior is claimed by GLY-bb.

`X` (non-canonical atoms contact the CG) is never a category: it has no
background frequency, and hit finding can never produce it. It is dropped from
numerator and denominator alike, and reported separately.
'''
import os
import math
import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

from ligand_vdgs.functions.vdg_struct_utils import (BB_LABEL, BB_LABELS,
                                                    NONCANONICAL_AA_LABEL)
from ligand_vdgs.tools.reference_values import (AA_sidechain_size,
                                                AA_background_freq, BB_HEAVY_ATOMS)

# Suffix marking a per-residue backbone category, e.g. 'GLY-bb'. Deliberately not
# '_': bucket file names join labels with '_' and callers split pair labels on it.
BB_CATEGORY_SUFFIX = '-' + BB_LABEL

# How backbone-mediated contacts enter the category set.
BB_MODES = ('per-residue', 'pooled', 'off')


def is_bb_category(label):
    '''True for the pooled backbone category and for any <AA>-bb category.'''
    return label in BB_LABELS or label.endswith(BB_CATEGORY_SUFFIX)


def category_table(bb_mode='per-residue'):
    '''dict: category label -> (bg_freq, size). See the module docstring.

    Size is sidechain heavy atoms only, which makes the priors an exact
    partition -- what the label claims. Glycine has none, so it has no
    sidechain category.
    '''
    if bb_mode not in BB_MODES:
        raise ValueError(f'Unknown bb_mode: {bb_mode!r}; expected one of {BB_MODES}')

    table = {aa: (AA_background_freq[aa], AA_sidechain_size[aa])
             for aa in AA_background_freq if AA_sidechain_size.get(aa, 0) > 0}
    if bb_mode == 'pooled':
        table[BB_LABEL] = (1.0, BB_HEAVY_ATOMS)
    elif bb_mode == 'per-residue':
        for aa, freq in AA_background_freq.items():
            table[aa + BB_CATEGORY_SUFFIX] = (freq, BB_HEAVY_ATOMS)
    return table


def category_priors(bb_mode='per-residue', norm_method='bg_weighted'):
    '''
    dict: category label -> prior probability, for one normalization method.

    'bg_weighted'  p = bg*size / sum(bg*size)   -- sums to 1 by construction
    'none'         p = bg / sum(bg)             -- no size term, so there is no
                                                   coherent way to weigh a
                                                   4-atom backbone against a
                                                   whole sidechain; requires
                                                   bb_mode='off'. Renormalized
                                                   over the 19 sidechain
                                                   categories that exist, which
                                                   removes the ~+0.079 nat
                                                   uniform shift glycine's
                                                   unclaimed 7.5% used to leave.
    '''
    table = category_table(bb_mode)
    if norm_method == 'none':
        if bb_mode != 'off':
            raise ValueError(
                "norm_method='none' has no size term, so a backbone category has no "
                "defensible prior against a sidechain one. Use bb_mode='off' with it, "
                "or use 'bg_weighted'.")
        denom = sum(bg for bg, _ in table.values())
        return {label: bg / denom for label, (bg, _) in table.items()}

    if norm_method != 'bg_weighted':
        raise ValueError(f'Unknown norm_method: {norm_method!r}')
    denom = sum(bg * size for bg, size in table.values())
    return {label: bg * size / denom for label, (bg, size) in table.items()}


def load_bucket_counts(npz_dir, bb_mode='off'):
    '''
    Read all .npz files in npz_dir and return (counts, extras) where counts maps a
    contact-category label -> non-redundant vdG cluster count. Cluster count is
    used instead of raw PDB occurrence to avoid bias from overrepresented protein
    families.

    bb_mode:
      'off'          backbone buckets are skipped entirely.
      'pooled'       backbone labels are kept as the single category 'bb'; a
                     mixed bucket keeps its name ('ALA_bb').
      'per-residue'  each backbone *slot* is attributed to a residue via its
                     cluster medoid's nr_scrr_resname and relabelled
                     '<AA>-bb'. A mixed bucket therefore fans out ('ALA_bb' ->
                     'ALA_GLY-bb', ...), with the parts re-sorted so a category
                     pair has one name.

                     The attribution is per *cluster*, not per observation, and
                     it is approximate: backbone clustering uses N/CA/C only, so
                     geometry has no power to separate donor residues and a bb
                     cluster may mix them. Measured on one large bb.npz (1,883
                     clusters / 24,697 members): 87% of clusters are
                     single-resname, and the medoid's resname covers 95% of its
                     cluster's members on average. Do not switch to member-based
                     attribution to tighten this -- it would mix units against
                     the cluster-counted sidechain categories.

    extras carries what was deliberately left out of both numerator and
    denominator: {'X': n, 'noncanonical_bb': n}. 'X' means non-canonical atoms are
    what contact the CG, so counting it as an observation of any amino acid --
    which is what a propensity is -- would be wrong, and it has no background
    frequency to divide by.
    '''
    if bb_mode not in BB_MODES:
        raise ValueError(f'Unknown bb_mode: {bb_mode!r}; expected one of {BB_MODES}')

    counts, extras = {}, {'X': 0, 'noncanonical_bb': 0}
    for fname in sorted(os.listdir(npz_dir)):
        if not fname.endswith('.npz'):
            continue
        bucket = fname[:-4]  # strip .npz
        parts = bucket.split('_')
        if NONCANONICAL_AA_LABEL in parts:
            with np.load(os.path.join(npz_dir, fname)) as npz:
                extras['X'] += len(npz['cluster_id'])
            continue
        has_bb = any(p in BB_LABELS for p in parts)
        if has_bb and bb_mode == 'off':
            continue

        with np.load(os.path.join(npz_dir, fname)) as npz:
            if not (has_bb and bb_mode == 'per-residue'):
                counts[bucket] = counts.get(bucket, 0) + len(npz['cluster_id'])
                continue
            # Deconvolute every bb slot into the residue that donated it.
            # nr_scrr_resname columns align with aa_bucket_parts slots.
            resnames = npz['nr_scrr_resname']
            bb_cols = [i for i, p in enumerate(parts) if p in BB_LABELS]
            for row in range(resnames.shape[0]):
                labels, ok = list(parts), True
                for col in bb_cols:
                    rn = str(resnames[row, col])
                    if rn not in AA_background_freq:
                        ok = False  # no prior for it; not an observation we can use
                        break
                    labels[col] = rn + BB_CATEGORY_SUFFIX
                if not ok:
                    extras['noncanonical_bb'] += 1
                    continue
                key = '_'.join(sorted(labels))
                counts[key] = counts.get(key, 0) + 1
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

    propensities = {}
    for label, count in counts.items():
        val = _log_enrichment(count, priors.get(label), total)
        if val is not None:
            propensities[label] = val

    return dict(sorted(propensities.items(), key=lambda x: x[1])), counts


def calc_aa_pair_propensities(nr_vdgs_dir, norm_method='bg_weighted', counts=None,
                              bb_mode='pooled'):
    '''Returns dict of "CAT1_CAT2" -> log-enrichment (unsorted).

    The pair prior is the independent-draw product of the two single-category
    priors: p(cat1, cat2) = sf * p1 * p2, with sf = 2 for a heterogeneous bucket
    (which absorbs both orderings) and 1 for a homogeneous one. Since each p
    already carries its size term, size normalization is applied once, not twice;
    sum(sf * p1 * p2) = (sum p)^2 = 1, so pair priors partition exactly like the
    single ones.'''
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
        symmetry_factor = 1 if parts[0] == parts[1] else 2
        val = _log_enrichment(count, symmetry_factor * p1 * p2, total)
        if val is not None:
            propensities[bucket] = val

    return propensities


def select_top_cgs(sim, n):
    '''Return indices of the top n CGs by maximum off-diagonal similarity score.'''
    sim_no_diag = sim.copy()
    np.fill_diagonal(sim_no_diag, np.nan)
    max_sim = np.nanmax(sim_no_diag, axis=1)
    top_idx = np.argsort(np.nan_to_num(max_sim, nan=-np.inf))[::-1][:n]
    return sorted(top_idx)


def similarity_linkage(sim, method='average'):
    '''Average-linkage clustering on 1 - sim distance; NaN cells -> maximally far (2.0).
    Returns the scipy linkage matrix Z.'''
    dist = np.clip(1.0 - sim, 0, None)
    np.fill_diagonal(dist, 0.0)
    dist = np.nan_to_num(dist, nan=2.0)
    condensed = squareform(dist, checks=False)
    return linkage(condensed, method=method)
