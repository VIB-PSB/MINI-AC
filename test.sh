#!/bin/bash
#
#SBATCH -p all # partition (queue)
#SBATCH -c 1 # number of cores
#SBATCH --mem 8G # memory pool for all cores
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

module load nextflow

# fetch test motif mappings (if not already present)
genome_wide_motif_mapping=tests/data/zma_v4/zma_v4_genome_wide_motif_mappings_chr1.bed
locus_based_motif_mapping=tests/data/zma_v4/zma_v4_locus_based_motif_mappings_5kbup_1kbdown_chr1.bed

if [ ! -f "$genome_wide_motif_mapping" ]; then
  curl -k -o $genome_wide_motif_mapping https://floppy.psb.ugent.be/index.php/s/NekMYztyxEnsQiY/download/zma_v4_genome_wide_motif_mappings_chr1.bed
fi

if [ ! -f "$locus_based_motif_mapping" ]; then
  curl -k -o $locus_based_motif_mapping https://floppy.psb.ugent.be/index.php/s/r2wQmFjPy79qSp7/download/zma_v4_locus_based_motif_mappings_5kbup_1kbdown_chr1.bed
fi

# run MINI-AC test pipeline
~/nf-test test
