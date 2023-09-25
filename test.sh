#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 8G # memory pool for all cores
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

module load nextflow

# fetch test motif mappings
curl -k -o tests/data/zma_v4/zma_v4_genome_wide_motif_mappings_chr1.bed \
           https://floppy.psb.ugent.be/index.php/s/NekMYztyxEnsQiY/download/zma_v4_genome_wide_motif_mappings_chr1.bed

curl -k -o tests/data/zma_v4/zma_v4_locus_based_motif_mappings_5kbup_1kbdown_chr1.bed \
           https://floppy.psb.ugent.be/index.php/s/r2wQmFjPy79qSp7/download/zma_v4_locus_based_motif_mappings_5kbup_1kbdown_chr1.bed

~/nf-test test
