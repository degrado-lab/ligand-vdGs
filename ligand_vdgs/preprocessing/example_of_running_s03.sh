qsub \
-t 1-65637 \
-o /wynton/home/degradolab/skt/DockDesign/smart-vdms/logfiles/run_probe \
-e /wynton/home/degradolab/skt/DockDesign/smart-vdms/logfiles/run_probe \
pipeline/s01_submit_probe.sh \
/wynton/home/degradolab/skt/DockDesign/databases/prepwizard_BioLiP2 \
/wynton/group/degradolab/skt/DockDesign/databases/probe_output \
/wynton/home/degradolab/skt/DockDesign/smart-vdms /wynton/home/degradolab/rkormos/probe/probe

