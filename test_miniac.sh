#!/bin/bash
#$ -l h_vmem=8G

module load nextflow/x86_64/21.10.6

nextflow -C mini_ac_example.config run mini_ac.nf --mode genome_wide --species maize -resume