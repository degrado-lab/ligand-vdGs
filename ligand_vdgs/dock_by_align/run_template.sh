#!/bin/bash

#$ -S /bin/bash
#$ -o /wynton/home/degradolab/skt/docking/docklogs
#$ -N a
#$ -cwd
#$ -j y               # tells system STDERR and STDOUT should be joined
#$ -R yes             #-- SGE host reservation
            

date
hostname

conda activate vdgs

python ../dock_by_align/dock_fgs.py -y dock_ymls/.yml

date
