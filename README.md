# MINI-AC

MINI-AC stands for Motif Informed Network Inference based on Accessible Chromatin, a method that combines accessible chromatin data from bulk or single-cell experiments with transcription factor binding site enrichment to learn gene regulatory networks (GRNs) in plants.
  
MINI-AC uses a ? license.

Reference: ?.

## **Pipeline summary**
1. Generation of background model for input accessible chromatin regions (ACRs).
2. Use of background ACR to compare with real ACR and compute motif enrichment statistics.
3. Inference of GRN based on motif enrichment results.
4. Generation of functional GRN by gene ontology (GO) enrichment.
5. Integration of data to generate final user-friendly output files.

It can be run for two different species, Arabidopsis and maize. Additionally, it can be run on two different modes depeding on the non-coding genomic space considered for motif mapping: genome-wide, where the whole non-coding genome is considered, and locus-based, where the neighboring sequences of each locus and introns are considered.

## **Inputs**
* **MINI-AC mode**: genome-wide or locus-based.
* **Species**: Arabidopsis or maize.
* **ACR files**: BED files containing genomic coordinates corresponding to accessible chromatin regions (minimal format of 3 columns: chromosome, start, stop).
* **Output folder**: Path where the results will be stored.
* **(Optional) DEGs file**: Tab-separated txt file with differential expression data associated with the input ACRs. First column must be gene ID.
* **(Optional) Expressed genes file**: One-column txt file with gene IDs for genes expressed in the biological context of the input ACRs, to filter the infered GRNs.

The pipeline will run in parallel for multiple ACR BED input files. The two optional input files can be provided individually for all the input ACR datasets, meaning that one single "DEG file" and one "Expressed genes file" can be provided for all the ACRs, or there can multiple and they need to be paired to each ACR file. This is done through the naming of the files. For more details consult [here](example).  

## **Outputs**
* **enrichment_stats folder** contains tab-separated files with table that summarize the motif enrichment results and associated statistics.
* **networks folder** contains tab-separated files with the accessibility based GRNs in an edge list format (minimal format of 2 columns with source or TF and target or TG node). If specified in parameters setting, there will be, per ACR file, one extra GRN filtered for the expressed genes
* **GO_enrichment folder** contains tab-separated files with GO enrichment of the previously generated GRN regulons. If filtering for expressed genes was specified, the GO enrichment is done on the filtered network
* **integrative_outputs folder** contains 3 xlsx Excel files and 2 tab-separated txt network files that integrate and summarize all the analysis: (1) integrative analysis per TF, (2) integrative analysis per motif, (3) integrative analysis of GO enrichment, (4) functional network edge list file, (5) nodes attributes files of the functional network

##   
A detailed overview on necessary input files and expected output files can be found [here](example).
## 
Requirements:

* [Nextflow version 21.10.6](https://www.nextflow.io/)
* [Singularity version 3.8.7-1.el7](https://sylabs.io/guides/3.0/user-guide/index.html)
* [Wget version 1.14](https://www.gnu.org/software/wget/)
* Download motif mapping files executing the following commands:

  ```shell
  wget
  wget
  wget
  wget
  ```


## 
Example on how to run it:

* Define paths, and the parameters setting in the [config file](docs/configuration_pipeline.md).

```shell
nextflow -C mini_ac.config run mini_ac.nf --mode genome_wide --species maize
```


