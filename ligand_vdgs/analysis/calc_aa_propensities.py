'''
Calculate the likelihood of each AA to interact with a given CG.
First, calculate propensities of each individual (single) AA to interact with each 
CG, and plot as a bar graph.
Second, calculate the propensities of each pair of AAs to interact with a given CG,
and plot as a heatmap.

In all cases, the propensity is normalized by the individual probability of the AA
simply appearing in the PDB. By default, the propensity is also normalized by the 
the size of the AA (# of atoms), but that can be turned off by using the 
`no-norm-AAsize` flag. If normalizing by AA size: 
- observed value is div. by # atoms in the AA (or AA pair) 
- expected value is div. by the avg # of atoms in an AA (or AA pair)

Example usage:
    python calc_aa_propensities.py \
        --CG <CG_NAME> \
        --vdglib-dir <VDGLIB_DIR> \
        --plot-single-aa-freqs \
        --plot-aa-pair-freqs 
'''
import os
import sys
import argparse
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

sys.path.append(os.path.join(os.path.dirname(__file__), '../../tools'))
from reference_values import AA_size, AA_background_freq, avg_size_of_aa_pair, avg_size_aa

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--CG', type=str, required=True, 
                        help='Name of CG.')
    parser.add_argument('--vdglib-dir', type=str, required=True, 
                        help='Path to directory that contains all vdGs.')
    parser.add_argument('--no-norm-AAsize', action='store_true',
                        help='Normalize by the AA size (# of atoms).')
    parser.add_argument('--plot-single-aa-freqs', action='store_true', 
                        help='Plot bar chart for propensities of individual AAs to '
                        'interact with a given CG.')
    parser.add_argument('--plot-aa-pair-freqs', action='store_true', 
                        help='Plot heatmap for propensities of AA pairs to '
                        'interact with a given CG.')
    return parser.parse_args()

def main():
    args = parse_args()
    CG = args.CG
    no_norm_AAsize = args.no_norm_AAsize
    if no_norm_AAsize:
        norm_by_AAsize= False
    else:
        norm_by_AAsize = True
    plot_single_aa_freqs = args.plot_single_aa_freqs
    plot_aa_pairs = args.plot_aa_pair_freqs
    print(f'norm_by_AAsize: {norm_by_AAsize}')
    print(f'plot_single_aa_freqs: {plot_single_aa_freqs}')
    print(f'plot_aa_pairs: {plot_aa_pairs}')
    vdglib_dir = args.vdglib_dir
    nr_vdgs_dir = os.path.join(vdglib_dir, CG, 'nr_vdgs')
    
    if plot_single_aa_freqs:
        sorted_single_aa_freqs = calc_stats_single_AA_with_CG(nr_vdgs_dir, CG, 
                                                              norm_by_AAsize)
        plot_bar_chart(sorted_single_aa_freqs, CG, norm_by_AAsize)
    if plot_aa_pairs:
        sorted_aa_pair_freqs = calc_stats_AA_pairs_with_CG(nr_vdgs_dir, CG, 
                                                           norm_by_AAsize)
        plot_heatmap(sorted_aa_pair_freqs, CG, norm_by_AAsize)

def calc_stats_AA_pairs_with_CG(nr_vdgs_dir, CG, norm_by_AAsize):
    # Calculate enrichment factors to determine which amino acids pairs are 
    # are favored or disfavored when interacting with a given CG.
    pair_aas_dir = os.path.join(nr_vdgs_dir, '2')
    aa_pair_freqs = {}

    # Count total number of AA pair files
    num_all_aa_pairs = sum(len(os.listdir(os.path.join(pair_aas_dir, i))) 
                           for i in os.listdir(pair_aas_dir))

    for aa_pair in os.listdir(pair_aas_dir):
        raw_counts = len(os.listdir(os.path.join(pair_aas_dir, aa_pair)))
        AA_list = aa_pair.split('_')

        # Calculate observed and expected frequencies
        if norm_by_AAsize:
            avg_aa_pair_size = avg_size_of_aa_pair()
            num_AA_atoms_pair = sum([AA_size.get(aa, 1) for aa in AA_list])  
                                          # .get() to avoids key errors
            if num_AA_atoms_pair == 0:  # Prevent division by zero
                continue
            observed = raw_counts / num_AA_atoms_pair
            expected = (num_all_aa_pairs * AA_background_freq.get(AA_list[0], 0) *
                        AA_background_freq.get(AA_list[1], 0)) / avg_aa_pair_size
        else:
            observed = raw_counts
            expected = num_all_aa_pairs * AA_background_freq.get(
                AA_list[0], 0) * AA_background_freq.get(AA_list[1], 0)

        # Ensure expected is nonzero before taking log 
        if expected > 0:
            aa_pair_freqs[aa_pair] = math.log(observed / expected)

    # Sort results
    sorted_aa_pair_freqs = dict(sorted(aa_pair_freqs.items(), 
                                       key=lambda item: item[1], reverse=True))
    return sorted_aa_pair_freqs

def plot_heatmap(sorted_aa_pair_freqs, CG, norm_by_AAsize):
    # Plot a heatmap of amino acid pair enrichment frequencies.
    # Missing values are set to NaN, appearing as white in the heatmap.
    
    aa_scores = {}  # Store total enrichment per AA
    for pair, value in sorted_aa_pair_freqs.items():
        aa1, aa2 = pair.split('_')
        aa_scores[aa1] = aa_scores.get(aa1, 0) + value
        aa_scores[aa2] = aa_scores.get(aa2, 0) + value

    # Sort amino acids based on their total enrichment score (descending order)
    aas = sorted(aa_scores, key=aa_scores.get, reverse=True)

    # Initialize heatmap matrix with NaN for missing values
    heatmap_data = np.full((len(aas), len(aas)), np.nan)

    # Fill the heatmap matrix
    for pair, freq in sorted_aa_pair_freqs.items():
        aa1, aa2 = pair.split('_')
        i, j = aas.index(aa1), aas.index(aa2)
        heatmap_data[i, j] = freq
        heatmap_data[j, i] = freq  # Ensure symmetry

    # Create colormap that treats NaNs as white
    cmap = plt.cm.get_cmap("YlOrRd").copy()  # Use the "YlOrRd" colormap
    cmap.set_bad(color='white')  # Set NaN values to white

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(4, 3))
    img = ax.imshow(heatmap_data, cmap=cmap, interpolation='nearest')

    # Format
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    cbar = plt.colorbar(img, fraction=0.046, pad=0.02)
    cbar.ax.tick_params(labelsize=5)  # Set colorbar tick label size
    cbar.set_label('Enrichment', fontsize=7)
    ax.set_xticks(np.arange(len(aas)))
    ax.set_yticks(np.arange(len(aas)))
    ax.set_xticklabels(aas, rotation=90, fontsize=6)
    ax.set_yticklabels(aas, fontsize=6)
    ax.set_title(f'AA Pair Propensities for CG: {CG}', pad=10, fontsize=6)
    plt.tight_layout()

    # Save plot
    filename = f'{CG}_aa_pair_freq{"_norm_AA_size" if norm_by_AAsize else ""}.png'
    plt.savefig(filename, dpi=400, transparent=True)


def calc_stats_single_AA_with_CG(nr_vdgs_dir, CG, norm_by_AAsize):
    # Calculate enrichment factor to determine which AAs are favored vs. disfavored 
    # for interacting with a given CG.
    single_aa_dir = os.path.join(nr_vdgs_dir, '1')
    single_aa_freqs = {}
    # Count all files in nr_vdgs/1/*/*pdb.gz
    num_all_single_aas = sum(len(os.listdir(os.path.join(single_aa_dir, i))) for i in 
                             os.listdir(single_aa_dir))
    
    # Calculate propensity scores
    for single_aa in os.listdir(single_aa_dir):
        raw_counts = len(os.listdir(os.path.join(single_aa_dir, single_aa)))
        # Calculate observed and expected (all interactions of CG [on the single 
        # AA level] * P(AA) )
        if norm_by_AAsize:
            # Get average size of AA for normalization 
            observed = raw_counts / AA_size[single_aa]
            avg_aa_size = avg_size_aa()
            expected = num_all_single_aas * AA_background_freq[single_aa] / avg_aa_size
        else:
            observed = raw_counts
            expected = num_all_single_aas * AA_background_freq[single_aa]
        # Calculate "propensity". 
        single_aa_freqs[single_aa] = math.log(observed/expected)
    
    sorted_single_aa_freqs = dict(sorted(single_aa_freqs.items(), key=lambda item: item[1]))
    return sorted_single_aa_freqs    

def plot_bar_chart(sorted_single_aa_freqs, CG, norm_by_AAsize):
    # Plot
    aas = list(sorted_single_aa_freqs.keys())
    freqs = list(sorted_single_aa_freqs.values())
    norm = plt.Normalize(min(freqs)-0.2, max(freqs)) # normalize freq for color mapping.
                                # shift the min -0.2 to avoid the dark purples
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
        # 2nd param of ax.text() sets the y-position of the text
        ax.text(i+.04, 0.03, aa, ha='center', va='bottom', fontsize=7, rotation=90)
    
    # Add a line for the avg freq, because that should be used a reference to see which
    # AAs are favored vs. disfavored, instead of using y=0 as the reference (bc the units,
    # after normalizing by # of atoms in AA, are arbitrary).
    avg_freq = np.mean(freqs)
    plt.axhline(y=avg_freq, color='red', linestyle='--', linewidth=1.5, 
                label=f'Avg enrichment: {avg_freq:.2f}')
    plt.legend(loc='upper left', fontsize=5)
    plt.tight_layout()
    if norm_by_AAsize:
        plt.savefig(f'{CG}_single_aa_freq_norm_AA_size.png', dpi=400, transparent=True)
    else:
        plt.savefig(f'{CG}_single_aa_freq.png', dpi=400, transparent=True)

if __name__ == "__main__":
    main()
