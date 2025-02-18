'''
Calculate the likelihood of each AA to interact with a given CG.
First, calculate propensities of each individual (single) AA to interact with each 
CG, and plot as a bar graph.
Second, calculate the propensities of each pair of AAs to interact with a given CG,
and plot as a heatmap.

In all cases, the propensity is normalized by the individual probability of the AA
simply appearing in the PDB. By default, the propensity is also normalized by the 
size of the AA (# of atoms), but that can be turned off by setting the 
--norm-AAsize flag to False.

Example usage:
    python calc_aa_propensities.py \
        --CG <CG_NAME> \
        --vdglib-dir <VDGLIB_DIR> \
        --aa-bkgrd-freq <AA_BKGRD_FREQ>
'''
import os
import sys
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

sys.path.append(os.path.join(os.path.dirname(__file__), '../tools'))
from references import AA_size, AA_background_freq

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--CG', type=str, required=True, 
                        help='Name of CG.')
    parser.add_argument('--vdglib-dir', type=str, required=True, 
                        help='Path to directory that contains the vdGs of this CG.')
    parser.add_argument('--norm-AAsize', action='store_true', default=True, 
                        help='Normalize by the size of the AA (# of atoms).')
    parser.add_argument('--plot-single-aa-freqs', action='store_true', default=True,
                        help='Plot bar chart for propensities of individual AAs to '
                        'interact with a given CG.')
    return parser.parse_args()

def main():
    args = parse_args()
    CG = args.CG
    norm_by_AAsize = args.norm_AAsize
    plot_single_aa_freqs = args.plot_single_aa_freqs
    vdglib_dir = args.vdglib_dir
    nr_vdgs_dir = os.path.join(vdglib_dir, 'nr_vdgs')
    
    if plot_single_aa_freqs:
        sorted_single_aa_freqs = calc_stats_single_AA_with_CG(nr_vdgs_dir, CG, 
                                                              norm_by_AAsize)
        plot_bar_chart(sorted_single_aa_freqs, CG)

def calc_stats_single_AA_with_CG(nr_vdgs_dir, CG, norm_by_AAsize):
    # Calculate enrichment factor to determine which AAs are favored vs. disfavored 
    # for interacting with a given CG.
    single_aa_dir = os.path.join(nr_vdgs_dir, '1')
    single_aa_freqs = {}
    # Count all files in nr_vdgs/1/*/*pdb.gz
    num_all_single_aas = sum(len(os.listdir(os.path.join(single_aa_dir, i))) for i in os.listdir(single_aa_dir))
    
    for single_aa in os.listdir(single_aa_dir):
        raw_counts = len(os.listdir(os.path.join(single_aa_dir, single_aa)))
        # Calculate observed
        observed = raw_counts / AA_size[single_aa] if norm_by_AAsize else raw_counts
        # Calculate expected (all interactions of CG [on the single AA level] * P(AA) )
        expected = num_all_single_aas * AA_background_freq[single_aa]
        # Calculate "propensity". Preferred is log(observed/expected) but that puts all 
        # values at negative numbers.
        single_aa_freqs[single_aa] = math.log(observed) / math.log(expected)
    
    sorted_single_aa_freqs = dict(sorted(single_aa_freqs.items(), key=lambda item: item[1]))
    return sorted_single_aa_freqs    

def plot_bar_chart(sorted_single_aa_freqs, CG):
    # Plot
    aas = list(sorted_single_aa_freqs.keys())
    freqs = list(sorted_single_aa_freqs.values())
    norm = plt.Normalize(min(freqs), max(freqs)) # normalize freq for color mapping
    colors = cm.viridis(norm(freqs))
    plt.figure(figsize=(4, 1.5), dpi=400)
    bars = plt.bar(np.arange(len(aas)), freqs, color=colors, edgecolor='black', 
                   width=0.7)
    
    # Format. 
    plt.ylabel('Enrichment', fontsize=5)
    plt.yticks(fontsize=4)
    # Remove the x-axis completely
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.3)
    ax.spines['bottom'].set_linewidth(1.3)
    ax = plt.gca()
    ax.set_xticks([])  # Remove the x-ticks
    ax.set_xticklabels([])  # Remove the x-tick labels

    # Add x-labels as text manually above the bars
    for i, aa in enumerate(aas):
        ax.text(i, 0.013, aa, ha='center', va='bottom', fontsize=5, rotation=90)
    
    # Add a line for the avg freq, because that should be used a reference to see which
    # AAs are favored vs. disfavored, instead of using y=0 as the reference (bc the units,
    # after normalizing by # of atoms in AA, are arbitrary).
    avg_freq = np.mean(freqs)
    plt.axhline(y=avg_freq, color='red', linestyle='--', linewidth=1.5, 
                label=f'Avg enrichment: {avg_freq:.2f}')
    plt.legend(loc='upper left', fontsize=5)
    
    plt.tight_layout()
    plt.savefig(f'{CG}_single_aa_propensities.png', dpi=400)

if __name__ == "__main__":
    main()
