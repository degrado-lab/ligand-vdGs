# find ligands with multiple "DONE" frags 
# frags need to have enough nr_vdg pdbs, and the ligands must have 
# >2 of those frags. 

import os
import pickle as pkl

logdir = '/wynton/home/degradolab/skt/docking/frag_sge_logs'
vdgdir = '/wynton/group/degradolab/skt/docking/databases/frag_vdg_lib/'

done_frags = []
for logfile in os.listdir(logdir):
    # Open the file to see if it has the word "DONE"
    with open(os.path.join(logdir, logfile), 'r') as f:
        # Read the file line by line
        lines = f.readlines()
        # Check if the word "DONE" is in the file
        done_count = 0
        for line in lines:
            if 'DONE' in line:
                done_frags.append(logfile.split('.')[0])

x = pkl.load(open('resources/database_frags_dict.pkl', 'rb'))


lignames_dict = {} # key=ligname, value=frags
for elements, v in x.items():
    for smiles, lignames in v.items():
        counts = 0
        if smiles not in done_frags:
            continue
        # count number of nr_vdg pdbs
        for aa_pair in os.listdir(os.path.join(vdgdir, smiles, 'nr_vdgs', str(2))):
            counts += len(os.listdir(os.path.join(vdgdir, smiles, 'nr_vdgs', str(2), aa_pair)))
        if counts < 2500:
            continue
        # if it passes the count, add to lignames_dict
        
        for ligname in lignames:
            if ligname not in lignames_dict.keys():
                lignames_dict[ligname] = []
            lignames_dict[ligname].append(smiles)

for k, v in lignames_dict.items():
    if len(v) > 2:
        print(k, v)
