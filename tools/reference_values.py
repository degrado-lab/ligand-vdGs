import numpy as np

AA_size = {'ALA': 5,
           'ARG': 11,
           'ASN': 8,
           'ASP': 8,
           'CYS': 6,
           'GLN': 9,
           'GLU': 9,
           'GLY': 4,
           'HIS': 10,
           'ILE': 8,
           'LEU': 8,
           'LYS': 9,
           'MET': 8,
           'PHE': 11,
           'PRO': 7,
           'SER': 6,
           'THR': 7,
           'TRP': 14,
           'TYR': 12,
           'VAL': 7}

AA_background_freq = \
          {'ALA': 0.08733094390532901,
           'ARG': 0.05087485521087214,
           'ASN': 0.041897950447092915,
           'ASP': 0.05879308405096171,
           'CYS': 0.012710544745734504,
           'GLN': 0.03609622691998023,
           'GLU': 0.06630531499164175,
           'GLY': 0.07549944048308671,
           'HIS': 0.024245591820242822,
           'ILE': 0.05869006081418494,
           'LEU': 0.09261670471415538,
           'LYS': 0.05449929225491284,
           'MET': 0.019573399146303915,
           'PHE': 0.04036970189037058,
           'PRO': 0.04583966858516758,
           'SER': 0.05872967862421906,
           'THR': 0.05488167184449434,
           'TRP': 0.013705138015159568,
           'TYR': 0.03489372480131626,
           'VAL': 0.07244700673477375}

def avg_size_aa():
    return np.mean(list(AA_size.values()))

def avg_size_of_aa_pair():
    size_aa_pairs = {}
    for aa1 in set(AA_size):
        for aa2 in set(AA_size):
            sorted_aa_pair = tuple(sorted([aa1, aa2])) # ensure only upper triangle
            pair_size = AA_size[aa1] + AA_size[aa2]
            if sorted_aa_pair not in size_aa_pairs.keys():
                size_aa_pairs[sorted_aa_pair] = pair_size
    avg_size_of_aa_pair = np.mean(list(size_aa_pairs.values()))
    return avg_size_of_aa_pair
