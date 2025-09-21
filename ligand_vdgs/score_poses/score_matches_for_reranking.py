import os
from pprint import pprint

## didn't work: 
## 8rlp, 7h69

#matches_dir = '/wynton/home/degradolab/skt/docking/vdg_matches/save_singleproc_1rm8_frags'
matches_dir = '/wynton/home/degradolab/skt/docking/vdg_matches/1fax_frags'
#matches_dir = '/wynton/home/degradolab/skt/docking/vdg_matches/1g2m_frags'
#matches_dir = '/wynton/home/degradolab/skt/docking/vdg_matches/4eyr'
#matches_dir = '/wynton/home/degradolab/skt/docking/vdg_matches/7h69'
#matches_dir = '/wynton/home/degradolab/skt/docking/vdg_matches/7hcg'
#matches_dir = '/wynton/home/degradolab/skt/docking/vdg_matches/8rlp'
#matches_dir = '/wynton/home/degradolab/skt/docking/vdg_matches/8sce'
rmsd_cut = 0.75

scores_dict = {} # key = struct name, val = num of matches per frag

if matches_dir.endswith('_frags'):
    structures = [matches_dir]
else:
    structures = sorted(os.listdir(matches_dir))

for structure in structures:
    structure_path = os.path.join(matches_dir, structure)
    structure = os.path.basename(structure_path)
    assert structure not in scores_dict.keys()
    scores_dict[structure] = {}
    for frag in os.listdir(structure_path):
        frag_path = os.path.join(structure_path, frag)
        num_frag_matches = 0
        for bsr_combo in os.listdir(frag_path):
            bsr_combo_path = os.path.join(frag_path, bsr_combo)
            for match in os.listdir(bsr_combo_path):
                # does this vdg match meet the rmsd cutoff?
                match_rmsd = float(match.split('_')[-1].removesuffix('.pdb'))
                if match_rmsd <= rmsd_cut:
                    num_frag_matches += 1
        assert frag not in scores_dict[structure].keys()
        scores_dict[structure][frag] = num_frag_matches

# --- Ensure all frags exist for every sample (fill missing with 0) ---
# Build the complete set of fragment keys across all samples
all_keys = sorted({k for frags in scores_dict.values() for k in frags})

# Fill missing fragments with 0 for every sample
for sample in scores_dict:
    for k in all_keys:
        scores_dict[sample].setdefault(k, 0)

# Step A: Compute average matches for each fragment across all samples
frag_totals = {}
frag_counts = {}

for sample, frags in scores_dict.items():
    for frag, count in frags.items():
        frag_totals[frag] = frag_totals.get(frag, 0) + count
        frag_counts[frag] = frag_counts.get(frag, 0) + 1

# When calculating the averages, make sure to include the samples w/ 0 matches
total_samples = len(scores_dict)
frag_avgs = {
    k: sum(scores_dict[sample][k] for sample in scores_dict) / total_samples
    for k in all_keys
}

# Step 1: Collect total scores
sample_totals = {}
for sample, frags in scores_dict.items():
    # instead of sum(frags.values())
    score = 0
    for frag, count in frags.items():
        avg = frag_avgs.get(frag, 1)  # avoid div by zero
        score += count / avg if avg > 0 else 0
    sample_totals[sample] = score

# Step 2: Sort the scores_dict by total descending
sorted_samples = sorted(scores_dict.keys(), key=lambda s: sample_totals[s], reverse=True)

# Step 3: Get all unique fragment keys
all_keys = set()
for subdict in scores_dict.values():
    all_keys.update(subdict)
all_keys = sorted(all_keys)

# Step 4: Print header with Total column
print(f"{'Pose':<20}" + "".join(f"{k:<12}" for k in all_keys) + f"{'Total':<10}")

# Step 5: Print sorted rows with Total column
for sample in sorted_samples:
    counts = scores_dict[sample]
    total = sample_totals[sample]
    row = f"{sample:<20}"
    for k in all_keys:
        row += f"{counts.get(k, 0):<12}"
    row += f"{total:<10.2f}"
    print(row)

# Step 6: Print Averages row at the bottom
avg_row = f"{'Average':<20}"
avg_total = 0
for k in all_keys:
    avg_val = frag_avgs.get(k, 0)
    avg_row += f"{avg_val:<12.2f}"  # average per fragment (raw counts)
    avg_total += avg_val
avg_row += f"{avg_total:<10.2f}"   # sum of averages as the "Total"
print(avg_row)
