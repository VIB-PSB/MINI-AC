# MINI-AC

MINI-AC stands for Motif-Informed Network Inference based on Accessible Chromatin, a method that combines accessible chromatin data from bulk or single-cell experiments with transcription factor binding site enrichment to learn gene regulatory networks (GRNs) in plants. The present README contains a tutorial of how to use MINI-AC, and how to modify its parameters to adapt it to user-desired settings.
  
MINI-AC uses a dual license to offer the distribution of the software under a proprietary model as well as an open source model.

## **Pipeline summary**
1. Generation of background model for input accessible chromatin regions (ACRs).
2. Use of background ACRs to compare with real ACRs and compute motif enrichment statistics.
3. Inference of a GRN based on motif enrichment results.
4. Generation of a functional GRN by gene ontology (GO) enrichment of the regulons.
5. Integration of data to generate informative, user-friendly output files.

Currently, two species are supported by MINI-AC: *Arabidopsis thaliana* and two maize genome versions (B73 RefGen_v4 and B73 RefGen_v5). Additionally, it can be run on two different modes depending on the non-coding genomic space considered for motif mapping:
* **genome-wide**: strategy where the whole non-coding genome is considered for motif mappings. It captures all the ACRs of the input dataset for the GRN prediction, which is adviced when working with species with long intergenic regions and distal regulatory elements, like maize for example.
* **locus-based**: strategy where the neighboring sequences within a pre-defined window of each locus, and introns are considered for motif mapping. It only captures the proximal ACRs of the input dataset within the pre-defined window, which can lead to missing distal ACRs in species with long intergenic regions and distal regulatory elements. However, it has the advantage of having a higher density of TFBS, which are mostly located close to the genes. The locus-based mode uses a "medium" non-coding genomic space, which corresponds, for each locus in the genome, to the 5kb upstream of the translation start site, the 1kb downstream of the translation end site, and the introns. However, for maize (but not for Arabidopsis; see publication), we generated two additional motif mapping files for the locus-based mode, that cover "large" (15kb upstream of the translation start site, the 2.5kb downstream of the translation end site, and the introns), and "small" (1kb upstream of the translation start site, the 1kb downstream of the translation end site, and the introns) non-coding genomic spaces. To use these files, check  the instructions [here](docs/pipeline_configuration.md).

A detailed overview of the necessary input files and expected output files can be found in this [example](example), done on **maize V4 with the genome-wide mode**, and using as input a single-cell-derived ACR dataset of mesophyll and bundle sheath.

## **Inputs**
* **MINI-AC mode**: genome-wide or locus-based.
* **Species**: Arabidopsis or maize (maize genome version 4 or 5).
* **ACR files**: BED files containing genomic coordinates corresponding to accessible chromatin regions (minimal format of 3 columns: chromosome, start, stop). The ACR files' coordinates **must** correspond to the genome versions of Araport11 for Arabidopsis and B73 RefGen_v4 and B73 RefGen_v5 for maize.
* **Output folder**: Path where the results will be stored.
* **(Optional) DEGs file**: Tab-separated txt file with differential expression data associated with the input ACRs. The only format requirements are that the first row has to be the header (column names), and the first column has to contain gene IDs. There is no requirement for the number of columns or content, although it should contain statistics associated to a DE analysis.
* **(Optional) Expressed genes file**: One-column txt file with gene IDs for genes expressed in the biological context of the input ACRs, to filter the inferred GRNs.

The pipeline will run in parallel for multiple ACR BED input files. The two optional input files can be provided individually for all the input ACR datasets, meaning that one single "DEG file" and one "Expressed genes file" can be provided for all the ACRs, or there can be multiple files which need to be paired to each ACR file. This is done through the naming of the files. For further details consult [here](example).  

## **Outputs**
* **enrichment_stats folder** contains tab-separated files with tables summarizing the motif enrichment statistics and numbers.
* **networks folder** contains tab-separated files with the accessibility-based GRNs in an edge list format (minimal format of 2 columns with source node or transcription factor (TF) being column 1 and the target node or target gene (TG) being column 2). If specified in parameters setting, there will be, per ACR file, one extra GRN file filtered for the expressed genes.
* **GO_enrichment folder** contains tab-separated files with GO enrichment results of the previously generated GRN regulons. If filtering for expressed genes was specified, the GO enrichment is done on the filtered network.
* **integrative_outputs folder** contains three xlsx Excel files and two tab-separated txt network files that integrate and summarize all the analysis: (1) integrative analysis per TF, (2) integrative analysis per motif, (3) integrative analysis of GO enrichment, (4) functional network edge list file, (5) nodes attributes files of the functional network. The last two files are meant to be used for network visualization using [Cytoscape](https://cytoscape.org/download.html).

## Requirements

* [Nextflow](https://www.nextflow.io/)
* [Singularity](https://sylabs.io/guides/3.0/user-guide/index.html)
* [wget](https://www.gnu.org/software/wget/)
* Motif mapping files. They need to be downloaded by executing the following commands on the **top-level directory of the repository**:
  
  For Arabidopsis

  ```
  wget https://zenodo.org/record/7974527/files/ath_genome_wide_motif_mappings.bed?download=1 -O data/ath/ath_genome_wide_motif_mappings.bed
  wget https://zenodo.org/record/7974527/files/ath_locus_based_motif_mappings_5kbup_1kbdown.bed?download=1 -O data/ath/ath_locus_based_motif_mappings_5kbup_1kbdown.bed
  ```
  For maize RefGen_v4
  ```
  wget https://zenodo.org/record/7974527/files/zma_genome_wide_motif_mappings.bed?download=1 -O data/zma_v4/zma_v4_genome_wide_motif_mappings.bed
  wget https://zenodo.org/record/7974527/files/zma_locus_based_motif_mappings_5kbup_1kbdown.bed?download=1 -O data/zma_v4/zma_v4_locus_based_motif_mappings_5kbup_1kbdown.bed
  ```
  For maize RefGen_v5
  ```
  wget https://zenodo.org/record/8386283/files/zma_v5_genome_wide_motif_mappings.bed?download=1 -O data/zma_v5/zma_v5_genome_wide_motif_mappings.bed
  wget https://zenodo.org/record/8386283/files/zma_v5_locus_based_motif_mappings_5kbup_1kbdown.bed?download=1 data/zma_v5/ -O data/zma_v5/zma_v5_locus_based_motif_mappings_5kbup_1kbdown.bed
  ```
 
NOTE: MINI-AC was developed using the following versions: Nextflow version 21.10.6, Singularity version 3.8.7-1.el7 and in a Sun Grid Engine (SGE) computer cluster.

## Usage

Define the paths with the input files and the desired parameters setting in the [configuration file](docs/pipeline_configuration.md), and run it executing the following Nextflow command:

```shell
nextflow -C mini_ac.config run mini_ac.nf --mode <genome_wide|locus_based> --species <arabidopsis|maize_v4|maize_v5>
```
 
Having problems running MINI-AC? Check the [FAQ](docs/FAQ.md).

## iCREs-based MINI-AC

Given the amount of resources available to profile regulatory DNA in maize, we curated a collection of integrated cis-regulatory elements (iCREs) by combining and comparing different CRE-profiling methods (details to be published).

We implemented a new framework in which it is possible to run MINI-AC given a list of maize genes. It works by retrieving the genomic coordinates of the iCREs associated with genes of interest, and submitting them to motif enrichment and GRN inference using the genome-wide mode of MINI-AC. iCREs-based MINI-AC can only be run for maize, and not for Arabidopsis. In addition, we offer different sets of iCREs that are used in the run: the "maxF1" (`maxf1`) set or the "all" (`all`) set. The first uses a set of putative CREs that is smaller but more precise (less false positives), while the second uses a more comprehensive and complete collection of maize putative CREs.

To download files with the genomic coordinates of the iCREs, the following commands should be executed on the **top-level directory of the repository**:

  For maize RefGen_v4
  ```
  wget https://zenodo.org/records/13143829/files/maxf1_icres_zma_v4.bed?download=1 -O data/icres/maxf1_icres_zma_v4.bed
  wget https://zenodo.org/records/13143829/files/all_icres_zma_v4.bed?download=1 -O data/icres/all_icres_zma_v4.bed
  ```
  For maize RefGen_v5
  ```
  wget https://zenodo.org/records/11192739/files/maxf1_icres_zma_v5.bed?download=1 -O data/icres/maxf1_icres_zma_v5.bed
  wget https://zenodo.org/records/11192739/files/all_icres_zma_v5.bed?download=1 -O data/icres/all_icres_zma_v5.bed
  ```
 

To run iCREs-based MINI-AC, the [configuration file](./mini_ac_icres.config) should be prepared as explained [here](./docs/pipeline_configuration.md). Only two parameters change in comparison to the regular MINI-AC runs. Instead of providing a BED file with ACR genomic coordinates, a list of gene IDs from the maize genome version V4 or V5 should be provided, as exemplified [here](./example/inputs/gene_set_files/UP_gene_set.txt). In addition, an iCREs set should be specified (`maxf1` or `all`). Next, the following Nextflow command should be executed:

```shell
nextflow -C mini_ac_icres.config run mini_ac_icres.nf --icres_set <all|maxf1> --species <maize_v4|maize_v5>
```

## Support

Should you encounter a bug or have any questions or suggestions, please [open an issue](https://github.com/VIB-PSB/MINI-AC/issues).

## Citation

When publishing results generated using MINI-AC, please cite:

Nicolás Manosalva Pérez, Camilla Ferrari, Julia Engelhorn, Thomas Depuydt, Hilde Nelissen, Thomas Hartwig, and Klaas Vandepoele. “MINI-AC: Inference of Plant Gene Regulatory Networks Using Bulk or Single-Cell Accessible Chromatin Profiles.” The Plant Journal 117, no. 1 (2024): 280–301. https://doi.org/10.1111/tpj.16483.

## Contact

If you have any questions, encounter issues, or want to contribute to this project, please feel free to reach out to Klaas Vandepoele (klaas.vandepoele@psb.vib-ugent.be).
