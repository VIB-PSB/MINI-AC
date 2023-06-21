#!/bin/bash
#$ -l h_vmem=5G

module load nextflow/x86_64/21.10.6

nextflow -C /group/transreg/niman/MINI-AC_runs/mini_ac.config run /group/transreg/niman/miniac_repo/MINI-AC/mini_ac.nf --species maize_v5 --mode locus_based