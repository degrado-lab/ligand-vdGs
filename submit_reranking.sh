#qsub -v YMLFILE=/wynton/home/degradolab/skt/docking/reranking_ymls/4eyr_true_pos.yml run_reranking.qsub
#qsub -v YMLFILE=/wynton/home/degradolab/skt/docking/reranking_ymls/4eyr_pred1.yml run_reranking.qsub
#qsub -v YMLFILE=/wynton/home/degradolab/skt/docking/reranking_ymls/4eyr_pred2.yml run_reranking.qsub
#qsub -v YMLFILE=/wynton/home/degradolab/skt/docking/reranking_ymls/4eyr_pred3.yml run_reranking.qsub
#qsub -v YMLFILE=/wynton/home/degradolab/skt/docking/reranking_ymls/4eyr_pred4.yml run_reranking.qsub
#qsub -v YMLFILE=/wynton/home/degradolab/skt/docking/reranking_ymls/4eyr_pred5.yml run_reranking.qsub
#qsub -v YMLFILE=/wynton/home/degradolab/skt/docking/reranking_ymls/4eyr_pred6.yml run_reranking.qsub
#qsub -v YMLFILE=/wynton/home/degradolab/skt/docking/reranking_ymls/4eyr_pred7.yml run_reranking.qsub
#qsub -v YMLFILE=/wynton/home/degradolab/skt/docking/reranking_ymls/4eyr_pred8.yml run_reranking.qsub
#qsub -v YMLFILE=/wynton/home/degradolab/skt/docking/reranking_ymls/4eyr_pred9.yml run_reranking.qsub

qsub -N _4eyr -v YMLFILE=/wynton/home/degradolab/skt/docking/ligand-vdGs/reranking_ymls/4eyr_all_preds.yml ligand_vdgs/score_poses/run_reranking.qsub

