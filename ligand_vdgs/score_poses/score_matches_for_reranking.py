import os
from pprint import pprint

matches_dirs = ['/wynton/home/degradolab/skt/docking/vdg_matches/' + x for x in 
               ['9r30', '9rci', '9gs6', '9r8y', '9y76', '9d5q', 
                '9gx3', '9fqy', '9d0s', '9foc',
                '9msr', '9rlo', '9jj5', '9rfe', '8yqe', '9ddg', 
                '9g0h', '7hc7', '9ei7', '9fs0', '9psz', 
               '9r9k', '9qff', '9d9i']]

rmsd_cuts = [.7, .6]

vdg_sizes_to_incl_in_count = [{"sizes": [1, 2], "weighted": False}, 
                              {"sizes": [1, 2], "weighted": True}, 
                              {"sizes": [1], "weighted": False}, 
                              {"sizes": [2], "weighted": False}]

ignore_frags = ['CC(C)N', 'CN(C)C', 'Cn(c)c', 'CN(c)C', 'cn(c)C',
                'cn(c)-c', 'cN(c)C', 'cn(c)c', 'cN(C)C']
# ignore CC(C)N?

print_counts = False

# Collect per-condition results across all matches_dirs (PDB structs)
results = []  # list of dicts summarizing each (rmsd, sizes, weighted) over all dirs

for matches_dir in matches_dirs:
    if print_counts:
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
            if print_counts:
                print(f'\nRMSD cutoff: {rmsd_cut}A', 'vdg size to incl in counting:', 
                  vdg_size_to_incl_in_count, 'weighted:', weighted,)
            
            num_rank_1_or_2 = 0
            num_rank_3 = 0
            num_true_positives = 0

            for structure in structures:
                if not os.path.isdir(os.path.join(matches_dir, structure)):
                    continue
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
                            match_rmsd = float(match.split('_')[-1].removesuffix('.pdb.gz'))
                            if match_rmsd > rmsd_cut:
                                continue
                            # exclude matches from the same pdb as the target
                            if pdbname in match[4:]: # the first 4 chars is the target pdb
                                #print(f'Excluding {match}')
                                continue
                            # determine which elements in spl can be converted to ints to 
                            # determine how many vdms there are in this vdg.
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

            # Ensure all frags exist for every sample (fill missing with 0)
            # Build the complete set of fragment keys across all samples
            all_keys = sorted({k for frags in scores_dict.values() for k in frags})

            # Fill missing fragments with 0 for every sample
            for sample_key in scores_dict:
                for k in all_keys:
                    scores_dict[sample_key].setdefault(k, 0)

            # Remove fragments that are all zero across all samples 
            all_keys = [
                frag for frag in all_keys
                if any(scores_dict[sample_key][frag] > 0 for sample_key in scores_dict)
            ]

            # Compute average matches for each fragment across all samples
            frag_totals = {}
            frag_counts = {}

            for sample_key, frags in scores_dict.items():
                for frag, count in frags.items():
                    if frag in all_keys:  # ensure only kept fragments are considered
                        frag_totals[frag] = frag_totals.get(frag, 0) + count
                        frag_counts[frag] = frag_counts.get(frag, 0) + 1

            # When calculating the averages, make sure to include the samples w/ 0 matches
            total_samples = len(scores_dict)
            frag_avgs = {
                k: sum(scores_dict[sample_key][k] for sample_key in scores_dict) / total_samples
                for k in all_keys
            }

            # Step 1: Collect total scores
            sample_totals = {}
            for sample_key, frags in scores_dict.items():
                score = 0
                # equally weight each fragment by averaging normalized contributions
                num_frags = len(all_keys) if all_keys else 1
                for frag in all_keys:
                    avg = frag_avgs.get(frag, 1)  # avoid div by zero
                    val = scores_dict[sample_key][frag]
                    contrib = (val / avg) if avg > 0 else 0
                    score += contrib
                sample_totals[sample_key] = score / num_frags

            # Step 2: Sort the scores_dict by total descending
            sorted_samples = sorted(scores_dict.keys(), key=lambda s: sample_totals[s], 
                                    reverse=True)

            # Step 3: Get all unique fragment keys
            all_keys = set()
            for subdict in scores_dict.values():
                all_keys.update(subdict)
            # keep only non-zero fragments
            all_keys = sorted([
                frag for frag in all_keys
                if any(scores_dict[sample_key][frag] > 0 for sample_key in scores_dict)])

            # Step 4: Print header with Total column
            if print_counts:
                print(f"{'Pose':<20}" + "".join(f"{k:<12}" for k in all_keys) + f"{'Total':<10}")

            # Step 5: Print sorted rows with Total column
            for sample_key in sorted_samples:
                counts = scores_dict[sample_key]
                total = sample_totals[sample_key]
                row = f"{sample_key:<20}"
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
                    num_true_positives += 1
                    rank = sorted_samples.index(s) + 1  # +1 for human-readable rank
                    if print_counts:
                        print(f"{s:<25} ranked #{rank} out of {len(sorted_samples)}")
                    if rank in [1, 2]:
                        num_rank_1_or_2 += 1
                    elif rank == 3:
                        num_rank_3 += 1
            else:
                print("\n(No samples with 'true_pos' found.)")

            # Record per-condition result for this matches_dir 
            results.append({
                "matches_dir": matches_dir,
                "rmsd": rmsd_cut,
                "sizes": tuple(vdg_size_to_incl_in_count),
                "weighted": weighted,
                "rank_1_2": num_rank_1_or_2,
                "rank_3": num_rank_3,
                "true_pos_total": num_true_positives,
                "rank_1_2_3_total": num_rank_1_or_2 + num_rank_3,})

# Aggregate across all matches_dirs and report best condition(s) 
if results:
    # aggregate by (rmsd, sizes, weighted)
    from collections import defaultdict

    agg = defaultdict(lambda: {"rank_1_2": 0, "rank_3": 0, "true_pos_total": 0})
    for r in results:
        key = (r["rmsd"], r["sizes"], r["weighted"])
        agg[key]["rank_1_2"] += r["rank_1_2"]
        agg[key]["rank_3"] += r["rank_3"]
        agg[key]["true_pos_total"] += r["true_pos_total"]

    # Compute combined metric (#1/#2/#3 hits)
    summary = []
    for (rmsd_cut, sizes, weighted), vals in agg.items():
        combined = vals["rank_1_2"] + vals["rank_3"]
        summary.append({
            "rmsd": rmsd_cut,
            "sizes": list(sizes),
            "weighted": weighted,
            "rank_1_2": vals["rank_1_2"],
            "rank_3": vals["rank_3"],
            "rank_1_2_3_total": combined,
            "true_pos_total": vals["true_pos_total"],})

    # pretty print summary
    print("\n=================== SUMMARY ACROSS ALL DIRECTORIES ===================")
    print(f"{'RMSD':<8}{'Sizes':<16}{'Weighted':<10}{'#1-2':<8}{'#3':<8}{'#1-3 total':<12}"
          f"{'TP total':<10}")
    for row in sorted(summary, key=lambda x: (-x["rank_1_2_3_total"], -x["rank_1_2"], 
                                              -x["rank_3"])):
        print(f"{row['rmsd']:<8}{str(row['sizes']):<16}{str(row['weighted']):<10}"
              f"{row['rank_1_2']:<8}{row['rank_3']:<8}{row['rank_1_2_3_total']:<12}"
              f"{row['true_pos_total']:<10}")

else:
    print("\n(No results collected; please check input directories.)")
