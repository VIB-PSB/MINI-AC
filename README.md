# MINI-AC

MINI-AC stands for Motif-Informed Network Inference based on Accessible Chromatin, a method that combines accessible chromatin data from bulk or single-cell experiments with transcription factor binding site enrichment to learn gene regulatory networks (GRNs) in plants.
  
MINI-AC uses an Open Source License.

## **Pipeline summary**
1. Generation of background model for input accessible chromatin regions (ACRs).
2. Use of background ACRs to compare with real ACRs and compute motif enrichment statistics.
3. Inference of a GRN based on motif enrichment results.
4. Generation of a functional GRN by gene ontology (GO) enrichment of the regulons.
5. Integration of data to generate informative, user-friendly output files.

Currently, two species are supported by MINI-AC: *Arabidopsis thaliana* and maize. Additionally, it can be run on two different modes depending on the non-coding genomic space considered for motif mapping:
* **genome-wide**: strategy where the whole non-coding genome is considered for motif mappings. It captures all the ACRs of the input dataset for the GRN prediction, which is adviced when working with species with long intergenic regions and distal regulatory elements, like maize for example.
* **locus-based**: strategy where the neighboring sequences within a pre-defined window of each locus, and introns are considered for motif mapping. It only captures the proximal ACRs of the input dataset within the pre-defined window, which can lead to missing distal ACRs in species with long intergenic regions and distal regulatory elements. However, it has the advantage of having a denser signal of TFBS, which are mostly located close to the genes.


A detailed overview of the necessary input files and expected output files can be found [here](example).


## **Inputs**
* **MINI-AC mode**: genome-wide or locus-based.
* **Species**: Arabidopsis or maize.
* **ACR files**: BED files containing genomic coordinates corresponding to accessible chromatin regions (minimal format of 3 columns: chromosome, start, stop). The ACR files' coordinates **must** correspond to the genome versions of Araport11 for Arabidopsis and AGPv4 for maize.
* **Output folder**: Path where the results will be stored.
* **(Optional) DEGs file**: Tab-separated txt file with differential expression data associated with the input ACRs. The only requirement is that the first column must be gene ID. There is no requirement for the number of columns or content, although it should contain statistics associated to a DE analysis.
* **(Optional) Expressed genes file**: One-column txt file with gene IDs for genes expressed in the biological context of the input ACRs, to filter the inferred GRNs.

The pipeline will run in parallel for multiple ACR BED input files. The two optional input files can be provided individually for all the input ACR datasets, meaning that one single "DEG file" and one "Expressed genes file" can be provided for all the ACRs, or there can be multiple files which need to be paired to each ACR file. This is done through the naming of the files. For further details consult [here](example).  

## **Outputs**
* **enrichment_stats folder** contains tab-separated files with tables summarizing the motif enrichment statistics and numbers.
* **networks folder** contains tab-separated files with the accessibility-based GRNs in an edge list format (minimal format of 2 columns with source node or transcription factor (TF) being column 1 and the target node or target gene (TG) being column 2). If specified in parameters setting, there will be, per ACR file, one extra GRN file filtered for the expressed genes.
* **GO_enrichment folder** contains tab-separated files with GO enrichment results of the previously generated GRN regulons. If filtering for expressed genes was specified, the GO enrichment is done on the filtered network.
* **integrative_outputs folder** contains three xlsx Excel files and two tab-separated txt network files that integrate and summarize all the analysis: (1) integrative analysis per TF, (2) integrative analysis per motif, (3) integrative analysis of GO enrichment, (4) functional network edge list file, (5) nodes attributes files of the functional network. The last two files are meant to be used for network visualization using [Cytoscape](https://cytoscape.org/download.html).

## Requirements:

* [Nextflow](https://www.nextflow.io/)
* [Singularity](https://sylabs.io/guides/3.0/user-guide/index.html)
* [Wget](https://www.gnu.org/software/wget/)
* Motif mapping files. They need to be downloaded by the executing the following commands on the **top-level directory of the repository**:

  ```
  wget https://zenodo.org/record/7836392/files/ath_genome_wide_motif_mappings.bed?download=1 -O data/ath/ath_genome_wide_motif_mappings.bed
  wget https://zenodo.org/record/7836392/files/ath_locus_based_motif_mappings_5kbup_1kbdown.bed?download=1 -O data/ath/ath_locus_based_motif_mappings_5kbup_1kbdown.bed
  wget https://zenodo.org/record/7836392/files/zma_genome_wide_motif_mappings.bed?download=1 -O data/zma/zma_genome_wide_motif_mappings.bed
  wget https://zenodo.org/record/7836392/files/zma_locus_based_motif_mappings_5kbup_1kbdown.bed?download=1 -O data/zma/zma_locus_based_motif_mappings_5kbup_1kbdown.bed
  ```
 
NOTE: MINI-AC was developed using the following versions: Nextflow version 21.10.6, Singularity version 3.8.7-1.el7 and in a Sun Grid Engine (SGE) computer cluster.

## How to run it

* Define the paths with the input files and the desired parameters setting in the [configuration file](docs/configuration_pipeline.md), and run it executing the following Nextflow command:

```shell
nextflow -C mini_ac.config run mini_ac.nf --mode genome_wide --species maize
```


