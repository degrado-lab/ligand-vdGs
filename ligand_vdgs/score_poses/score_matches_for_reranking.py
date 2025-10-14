import os
from pprint import pprint

matches_dirs = ['/wynton/home/degradolab/skt/docking/vdg_matches/' + x for x in 
               ['4eyr', '7h69', '7hc7', '7hcg', '8rlp', '8sce', '9foc', '9foe', '9d0s', 
               '9gs6', '9rlq', '9msr', '9gw3', '9rci', '9fqy', '9r8y', '9r30']]

rmsd_cuts = [.7, .6]

vdg_sizes_to_incl_in_count = [{"sizes": [1, 2], "weighted": False}, 
                              {"sizes": [1, 2], "weighted": True}, 
                              {"sizes": [1], "weighted": False}, 
                              {"sizes": [2], "weighted": False}]

ignore_frags = ['CC(C)N', 'CN(C)C', 'Cn(c)c', 'CN(c)C']

print_counts = False

for matches_dir in matches_dirs:
    print('\n===================', matches_dir, '===================')
    for sample in vdg_sizes_to_incl_in_count:
        vdg_size_to_incl_in_count = sample["sizes"]
        weighted = sample["weighted"]
        for rmsd_cut in rmsd_cuts:
            scores_dict = {} # key = struct name, val = num of matches per frag

            if matches_dir.endswith('_frags'):
                structures = [matches_dir]
            else:
                structures = sorted(os.listdir(matches_dir))
            print(f'\nRMSD cutoff: {rmsd_cut}A', 'vdg size to incl in counting:', 
                  vdg_size_to_incl_in_count, 'weighted:', weighted,)
            for structure in structures:
                structure_path = os.path.join(matches_dir, structure)
                structure = os.path.basename(structure_path)
                pdbname = structure[:4]
                assert structure not in scores_dict.keys()
                scores_dict[structure] = {}
                for frag in os.listdir(structure_path):
                    if frag in ignore_frags:
                        continue
                    frag_path = os.path.join(structure_path, frag)
                    num_frag_matches = 0
                    for bsr_combo in os.listdir(frag_path):
                        bsr_combo_path = os.path.join(frag_path, bsr_combo)
                        for match in os.listdir(bsr_combo_path):
                            # does this vdg match meet the rmsd cutoff?
                            match_rmsd = float(match.split('_')[-1].removesuffix('.pdb'))
                            if match_rmsd > rmsd_cut:
                                continue
                            # exclude matches from the same pdb as the target
                            if pdbname in match[4:]: # the first 4 chars is the target pdb
                                #print(f'Excluding {match}')
                                continue
                            # determine which elements in spl can be converted to ints to determine 
                            # how many vdms there are in this vdg.
                            if len(bsr_combo.split('_')) == 9:
                                subset_size = 2
                            elif len(bsr_combo.split('_')) == 5:
                                subset_size = 1
                            else:
                                raise ValueError(f'Cannot determine subset size from {bsr_combo}')
                            if subset_size in vdg_size_to_incl_in_count:
                                if weighted:
                                    if subset_size == 1:
                                        num_frag_matches += 1 * 0.5
                                    elif subset_size == 2:
                                        num_frag_matches += 1 
                                else:
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

            # --- Remove fragments that are all zero across all samples ---
            all_keys = [
                frag for frag in all_keys
                if any(scores_dict[sample][frag] > 0 for sample in scores_dict)
            ]

            # Step A: Compute average matches for each fragment across all samples
            frag_totals = {}
            frag_counts = {}

            for sample, frags in scores_dict.items():
                for frag, count in frags.items():
                    if frag in all_keys:  # ensure only kept fragments are considered
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
                # equally weight each fragment by averaging normalized contributions
                num_frags = len(all_keys) if all_keys else 1
                for frag in all_keys:
                    avg = frag_avgs.get(frag, 1)  # avoid div by zero
                    val = scores_dict[sample][frag]
                    contrib = (val / avg) if avg > 0 else 0
                    score += contrib
                sample_totals[sample] = score / num_frags

            # Step 2: Sort the scores_dict by total descending
            sorted_samples = sorted(scores_dict.keys(), key=lambda s: sample_totals[s], reverse=True)

            # Step 3: Get all unique fragment keys
            all_keys = set()
            for subdict in scores_dict.values():
                all_keys.update(subdict)
            # keep only non-zero fragments
            all_keys = sorted([
                frag for frag in all_keys
                if any(scores_dict[sample][frag] > 0 for sample in scores_dict)
            ])

            # Step 4: Print header with Total column
            if print_counts:
                print(f"{'Pose':<20}" + "".join(f"{k:<12}" for k in all_keys) + f"{'Total':<10}")

            # Step 5: Print sorted rows with Total column
            for sample in sorted_samples:
                counts = scores_dict[sample]
                total = sample_totals[sample]
                row = f"{sample:<20}"
                for k in all_keys:
                    row += f"{counts.get(k, 0):<12}"
                row += f"{total:<10.2f}"
                if print_counts:
                    print(row)

            # Step 6: Print Averages row at the bottom
            avg_row = f"{'Average':<20}"
            avg_total = 0
            for k in all_keys:
                avg_val = frag_avgs.get(k, 0)
                avg_row += f"{avg_val:<12.2f}"  # average per fragment (raw counts)
                avg_total += avg_val
            avg_row += f"{avg_total:<10.2f}"   # sum of averages as the "Total"
            if print_counts:
                print(avg_row)

            # Step 7: print ranks of all "true_pos" samples
            true_pos_samples = [s for s in sorted_samples if "true_pos" in s]

            if true_pos_samples:
                for s in true_pos_samples:
                    rank = sorted_samples.index(s) + 1  # +1 for human-readable rank
                    print(f"{s:<25} ranked #{rank} out of {len(sorted_samples)}")
            else:
                print("\n(No samples with 'true_pos' found.)")
