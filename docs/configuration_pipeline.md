# MINI-AC parameters setting and pipeline configuration

This document contains an overview of the pipeline parameters. They should be defined using the Nextflow configuration file (see below for details). 

## Input paths

MINI-AC has 4 main inputs that need to be given as paths or folder names, with two them being optional:

* **ACR files**: Path of the folder with the BED files containing genomic coordinates corresponding to accessible chromatin regions (minimal format of 3 columns: chromosome, start, stop). This path should be given to the parameter ```ACR_dir```.
* **Output folder**: Path where the results will be stored. This path should be given to the parameter ```OutDir```.
* **(Optional) DEGs file**: Path of folder with tab-separated txt files with differential expression data associated with the input ACRs. First column must be gene ID. It can be one DEGs file per input ACR file, or paired DEGs files-ACR files. For more details see [inputs format example](../example/). This path should be given to the parameter ```DE_genes_dir```.
* **(Optional) Expressed genes file**: Path of folder with one-column txt files with gene IDs for genes expressed in the biological context of the input ACRs, to filter the inferred GRNs. It can be one Expression file per input ACR file, or paired Expressed genes files-ACR files. For more details see [inputs format example](../example/). This path should be given to the parameter ```Set_genes_dir```.

## Input parameters

MINI-AC has several optional parameters that affect the output and some aspects of the network inference process:

### DEGs and Expressed genes files parameters

* **Species**: ```--species arabidopsis``` (command line) or ```species = "arabidopsis"``` (configuration file) for Arabidopsis, ```--species maize_v4``` (command line) or ```species = "maize_v4"``` (configuration file) for maize genome version 4, and ```--species maize_v5``` (command line) or ```species = "maize_v5"``` (configuration file) for maize genome version 5.

* **MINI-AC mode**: ```--mode genome_wide``` (command line) or ```mode = "genome_wide"``` (configuration file) for the genome-wide mode, and ```--mode locus_based``` (command line) or ```mode = "locus_based"``` (configuration file) for the locus-based mode.

* **DEGs parameters**: Since providing DEGs files is optional, it needs to be specified if the path with the DEGs files is available with the parameter ```DE_genes``` set to ```DE_genes = true``` or ```DE_genes = false```. Additionally, if there is only one DEG file for all the input ACRs, you need set the parameter ```One_DE_set``` to ```One_DE_set = true```, and to ```One_DE_set = false``` if otherwise.


* **Expressed genes files parameters**: Since providing Expressed genes files is optional, it needs to be specified if the path with the Expressed genes files is available with the parameter ```Filter_set_genes``` set to ```Filter_set_genes = true``` or ```Filter_set_genes = false```. Additionally, if there is only one Expression file for all the input ACRs, you need set the parameter ```Filter_set_genes``` to ```Filter_set_genes = true```, and to ```Filter_set_genes = false``` if otherwise.

### GRN inference parameters

* **Motif enrichment p-value cut-off**: This is the p-value cut-off that determines which motifs are enriched and used for GRN building. We do not recommend changing this parameter. It has been internally pre-defined for each MINI-AC mode based on the p-value cut-offs with a false discovery rate of 0 (see [publication](https://doi.org/10.1111/tpj.16483)). If wished, however, this p-value can be overwritten in the configuration file by setting the parameter ```P_val``` to whatever value (see below) or in the command line options. For example: ```nextflow -C mini_ac.config run mini_ac.nf --mode genome_wide --species maize_v4 --P_val 0.05```

* **Overlap criteria parameter**: By default, MINI-AC computes motif enrichment counting the motif matches within ACRs. This, however, is difficult if the ACRs are shorter than or of similar size to the motifs, which is the case of footprints. In this case, we observed that counting the absolute base-pair overlap is useful. Therefore, in case of using footprints or short ACRs (high resolution), we recommend setting the parameter ```Bps_intersect = true```. Otherwise it should be kept ```Bps_intersect = false```.

* **Annotation of second closest gene in genome-wide mode**: The parameters ```Second_gene_annot``` and ```Second_gene_dist``` are only taken into account by the genome-wide mode. In the genome-wide mode the motif matches are annotated to the closest gene, but in genomes like maize, there are very distal regulatory elements that regulate non-neighboring genes. Although we showed in the original publication that this does not improve results, we give the possibility of annotating the second closest genes that are within a certain distance from the motif match. To activate this option the parameter ```Second_gene_annot``` should be set to ```Second_gene_annot = true```. If so, the parameter ```Second_gene_dist``` should be used to set the specific distance cut-off (in absolute base-pairs) at which the second-closest gene has to be from the motif match in order to be assigned as target gene.

## **Nextflow configuration file**
The configuration or "config" file, is a file that Nextflow uses to manage and specify the inputs and parameters settings of a pipeline. For more details read [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html). Here we review the main aspects of the configuration file when running MINI-AC.

### Give parameters values to configuration file

To set the above-mentioned parameters of the pipeline in the configuration file, here's a code snippet with the default and the recommended settings:

```nextflow
params {

    //// Output folder
    OutDir = "/absolute/path/to/output/directory"

    //// Required input
    ACR_dir = "/absolute/path/to/acrs/directory"

    //// Optional input
    // Differential expression data
    DE_genes = true
    DE_genes_dir = "/absolute/path/to/degs_genes/directory"
    One_DE_set = true
    // Expression data
    Filter_set_genes = true
    Set_genes_dir = "/absolute/path/to/ex_files/directory"
    One_filtering_set = true


    //// Prediction parameters
    Bps_intersect = false
    // P_val = 0.01 // This is commented because we do not recommend changing it.


    //// Prediction parameters only genome-wide
    Second_gene_annot = false
    Second_gene_dist = 500
}
```
### Setting Work directory of the pipeline
	
The Work directory is where the temporary files are created and where Nextflow stores the files of the different processes. By default this directory is created in the folder where the pipeline is executed, but we recommend to set it to a scratch or tmp folder. To set it in the configuration file, the following code line should be added and edited:

```nextflow
workDir = '/absolute/path/to/work/dir'
```
### Enabling Singularity container 
[Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) is a container platform. It allows to run the Nextflow pipeline with a pre-determined environment that ensures reproducibility. We created a Docker image with the necessary dependencies for MINI-AC to run, without need for the user to install any of them. We strongly recommend to always run MINI-AC using the specified singularity container. For that, the following code needs to be included in the configuration file:

```nextflow
process.container = "vibpsb/mini-ac:latest"
singularity {
    enabled = true
    cacheDir = "singularity_cache"
    autoMounts = true
}
```
Sometimes the temporary directory used by Singularity is not in the same root path as in the pipeline, which can cause Singularity to struggle to find it. In this case, add the ```runOptions``` line below with the absolute path to the tmp folder. To know the absolute path to the tmp folder in linux execute in the command line ```echo $TMPDIR```. Then add it as shown below.

```nextflow
process.container = "vibpsb/mini-ac:latest"
singularity {
    enabled = true
    cacheDir = "singularity_cache"
    autoMounts = true
    runOptions = "--bind /absolute/path/to/tmp/folder"
}
```
### Pipeline executor

[The executor of the pipeline](https://www.nextflow.io/docs/latest/executor.html) is the system where the pipeline processes run and supervises its execution. It can be a computer, a cluster resource manager, or the cloud. In the configuration file it can be specified what is the executor of the pipeline. To execute it in a normal computer locally, the code below should be added in the configuration file:

```nextflow
executor {
    name = 'local'
}
```

MINI-AC was developed in an SGE computer cluster, for which we used the configuration below. This was used to run the genome-wide mode on maize using an input dataset of ~600,000 MOA-seq peaks. For smaller datasets, the memory values can be further reduced. Addionally, for Arabidopsis, a species with a smaller genome, less memory can also be used.

```nextflow
executor {
    name = 'sge'
    queueSize = 25
}


process {
    withName: get_ACR_shufflings {
        clusterOptions = '-l h_vmem=4G'
    }
    withName: getStats {
        clusterOptions = '-l h_vmem=10G'
    }
    withName: getStats_bps {
        clusterOptions = '-l h_vmem=50G'
    }
    withName: getNetwork {
        clusterOptions = '-l h_vmem=20G'
    }
    withName: filterSetOfGenes {
        clusterOptions = '-l h_vmem=5G'
    }
    withName: GOenrichment {
        clusterOptions = '-l h_vmem=5G'
    }
    withName: getIntegrativeOutputs {
        clusterOptions = '-l h_vmem=3G'
    }
    }
```

## **Modification of the internal files of the pipeline**

### Priority of configuration sources

The MINI-AC Nextflow pipeline contains a set of pre-defined parameter files specified within the main pipeline script (mini_ac.nf). This is because they are fixed data files for each MINI-AC mode and specie's genome version. However, there are cases where some of this files might want to be changed by the user. Nextflow allows to easily change these parameter files, either through the command line options or in the configuration file, thanks to a hierarchical prioritization of the configuration sources:

1. Parameters specified on the command line (--something value)
3. Config file specified using the -C mini_ac.config option
4. The config file named nextflow.config in the current directory
5. The config file named nextflow.config in the workflow project directory
7. Values defined within the pipeline script itself (e.g. main.nf)

Therefore, if the user wishes to change any of these parameters, it possible either through the command line options, or in the config file.

There are mainly two cases in which the user might want to alter the internal MINI-AC files, which are explained below.

### Modification of the motif mapping file for the locus-based mode of maize

By default, the maize MINI-AC locus-based mode (for both genome versions) runs on the "medium" non-coding genomic space, which corresponds, for each locus in the genome, to the 5kb upstream of the translation start site, the 1kb downstream of the translation end site, and the introns. However, we generated two additional motif mapping files for the locus-based mode of maize, that cover "large" (15kb upstream of the translation start site, the 2.5kb downstream of the translation end site, and the introns), and "small" (1kb upstream of the translation start site, the 1kb downstream of the translation end site, and the introns) non-coding genomic spaces. For Arabidopsis only the "medium" non-coding genomic space motif mapping file was generated because it already covers 73.5% of the whole non-coding genomic psace (see publication). To use these files, first they need to be downloaded, and then, the corresponding parameters of the motif mapping file (```MotMapsFile_lb```) and the non-coding genomic space coordinates file (```Promoter_file```) should be modified either on the command line or in the configuration file.

To download the maize "large" motif mapping file and coordinates of the "large" non-coding genomic space:

  For maize RefGen_v4 large locus-based mode files
```
wget https://zenodo.org/record/7974527/files/zma_locus_based_motif_mappings_15kbup_2.5kbdown.bed?download=1 -O data/zma_v4/zma_v4_locus_based_motif_mappings_15kbup_2.5kbdown.bed
wget https://zenodo.org/record/7974527/files/zma_promoter_15kbup_2.5kbdown_sorted.bed?download=1 -O data/zma_v4/zma_v4_promoter_15kbup_2.5kbdown_sorted.bed
```
  For maize RefGen_v5 large locus-based mode files
```
wget https://zenodo.org/record/8386283/files/zma_v5_locus_based_motif_mappings_15kbup_2.5kbdown_sorted.bed?download=1 -O data/zma_v5/zma_v5_promoter_15kbup_2.5kbdown_sorted.bed
wget https://zenodo.org/record/8386283/files/zma_v5_promoter_15kbup_2.5kbdown_sorted.bed?download=1 -O data/zma_v5/zma_v5_promoter_15kbup_2.5kbdown_sorted.bed
```

To download the maize "small" motif mapping file and coordinates of the "small" non-coding genomic space:

  For maize RefGen_v4 small locus-based mode files
```
wget https://zenodo.org/record/7974527/files/zma_locus_based_motif_mappings_1kbup_1kbdown.bed?download=1 -O data/zma_v4/zma_v4_locus_based_motif_mappings_1kbup_1kbdown.bed
wget https://zenodo.org/record/7974527/files/zma_promoter_1kbup_1kbdown_sorted.bed?download=1 -O data/zma_v4/zma_v4_promoter_1kbup_1kbdown_sorted.bed
```
  For maize RefGen_v5 small locus-based mode files
```
wget https://zenodo.org/record/8386283/files/zma_v5_locus_based_motif_mappings_1kbup_1kbdown_sorted.bed?download=1 -O data/zma_v5/zma_v5_locus_based_motif_mappings_1kbup_1kbdown.bed
wget https://zenodo.org/record/8386283/files/zma_v5_promoter_1kbup_1kbdown_sorted.bed?download=1 -O data/zma_v5/zma_v5_promoter_1kbup_1kbdown_sorted.bed
```
Then (using the "small" definition as example), change the parameters on the command line:

```
nextflow -C mini_ac.config run mini_ac.nf --mode locus_based --species maize_v4 --MotMapsFile_lb data/zma_v4/zma_v4_locus_based_motif_mappings_1kbup_1kbdown.bed --Promoter_file data/zma_v4/zma_v4_promoter_1kbup_1kbdown_sorted.bed
```
or add them to the configuration file, along with the other parameters:

```nextflow
params {
    /// [Other parameters...]
    MotMapsFile_lb = "$projectDir/data/zma_v4/zma_v4_locus_based_motif_mappings_1kbup_1kbdown.bed"
    Promoter_file = "$projectDir/data/zma_v4/zma_v4_promoter_1kbup_1kbdown_sorted.bed"
    /// [Other parameters...]
}
```

### Providing custom gene ontology (GO)-gene annotation file

To perform the functional network analysis, an internal gene-GO annotation file is used for each species. They were obtained as described in [here](../data/ath/README_MINI_ath_2021.1_motif_mappings.txt) for Arabidopsis and in [here](../data/zma_v4/README_MINI_zma_v4_2021.1_motifsMapping.txt) and [here](../data/zma_v5/README_MINI_zma_v5_2023.1_motifsMapping.txt) for maize. However, if the user wants to use a custom GO-gene file, the following parameters should be medified either on the command line or in the configuration file.

```
nextflow -C mini_ac.config run mini_ac.nf --mode locus_based --species maize_v4 --Feature_file custom_go_gene.txt
```

```nextflow
params {
    /// [Other parameters...]
    Feature_file = "custom_go_gene.txt"
    /// [Other parameters...]
}
```
It is important, however, to make sure that the format is correct. The GO terms should be extended for parental terms, and this file should contain two tab-separated columns (no header),  where the first column is the GO ID, and the second column is the gene ID, as shown [here](../data/zma_v4/zma_v4_go_gene_file.txt). It is vital that the gene IDs are either on Araport11 or AGPv4/NAM5.0.

This same principle can also be applied to other parameters that the user wants to change.
