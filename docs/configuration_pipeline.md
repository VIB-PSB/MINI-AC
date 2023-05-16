# MINI-AC parameters setting and pipeline configuration

This document contains an overview of the pipeline parameters. They should be defined using the Nextflow configuration file (see below for details). 

## Input paths

MINI-AC has 4 main inputs that need to be given as paths or folder names, with two them being optional:

* **ACR files**: Path of the folder with the BED files containing genomic coordinates corresponding to accessible chromatin regions (minimal format of 3 columns: chromosome, start, stop). This path should be given to the parameter ```ACR_dir```.
* **Output folder**: Path where the results will be stored. This path should be given to the parameter ```OutDir```.
* **(Optional) DEGs file**: Path of folder with tab-separated txt files with differential expression data associated with the input ACRs. First column must be gene ID. It can be one DEGs file per input ACR file, or paired DEGs files-ACR files. For more details see [inputs format example](../example/). This path should be given to the parameter ```DE_genes_dir```.
* **(Optional) Expressed genes file**: Path of folder with one-column txt files with gene IDs for genes expressed in the biological context of the input ACRs, to filter the infered GRNs. It can be one Expression file per input ACR file, or paired Expressed genes files-ACR files. For more details see [inputs format example](../example/). This path should be given to the parameter ```Set_genes_dir```.

## Input parameters

MINI-AC has several optional parameters that affect the output and some aspects of the network inference process:

### DEGs and Expressed genes files parameters

* **DEGs parameters**: Since providing DEGs files is optional, it needs to be specified if the path with the DEGs files is available with the parameter ```DE_genes``` set to ```DE_genes = true``` or ```DE_genes = false```. Additionally, if there is only one DEG file for all the input ACRs, you need set the parameter ```One_DE_set``` to ```One_DE_set = true```, and to ```One_DE_set = false``` if otherwise.


* **Expressed genes files parameters**: Since providing Expressed genes files is optional, it needs to be specified if the path with the Expressed genes files is available with the parameter ```Filter_set_genes``` set to ```Filter_set_genes = true``` or ```Filter_set_genes = false```. Additionally, if there is only one Expression file for all the input ACRs, you need set the parameter ```Filter_set_genes``` to ```Filter_set_genes = true```, and to ```Filter_set_genes = false``` if otherwise.

### GRN inference parameters

* **Motif enrichment p-value cut-off**: This is the p-value cut-off that determines which motifs are enriched and used for GRN building. We do not recommend changing this parameter. It has been internally pre-defined for each MINI-AC mode based on the p-value cut-offs with a false disocvery rate of 0 (see publication). If wished, however, this p-value can be overwritten in the configuration file by setting the parameter ```P_val``` to whatever value (see below) or in the command line options. For example: ```nextflow -C mini_ac.config run mini_ac.nf --mode genome_wide --species maize --P_val 0.05```

* **Overlap criteria parameter**: By default, MINI-AC computes motif enrichment counting the motif matches within ACRs. This, however, is difficult if the ACRs are shorter than or of similar size to the motifs, which is the case of footprints. In this case, we observed that counting the absolute base-pair overlap is useful. Therefore, in case of using footprints or short ACRs (high resolution), we recommend setting the parameter ```Bps_intersect = true```. Otherwise it should be kept ```Bps_intersect = false```.

* **Annotation of second closest gene in genome-wide mode**: The parameters ```Second_gene_annot``` and ```Second_gene_dist``` are only taken into account by the genome-wide mode. In the genome-wide mode the motif matches are annotated to the closest gene, but in genomes like maize, there are very distal regulatory elements that regulate non-neigboring genes. Although we showed in the original publication that this does not improve results, we give the possibility of annotating the second closest genes that are within a certain distance from the motif match. To activate this option the parameter ```Second_gene_annot``` should be set to ```Second_gene_annot = true```. If so, the parameter ```Second_gene_dist``` should be used to set the specific distance cut-off (in absolute base-pairs) at which the second-closest gene has to be from the motif match in order to be assigned as target gene.

## **Nextflow configuration file**
The configuration or "config" file, is a file that Nextflow uses to manage and specify the inputs and parameters settings of a pipeline. For more details read [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html). Here we review the main aspcects of the Configuration file when running MINI-AC.

### Give parameters values to Configuration file

To set the above-mentioned parameters of the pipeline in the Configuration file, here's a code snippet with the default and the recommended settings:

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
	
The Work directory is where the temporary files are created and where Nextflow stores the files of the different processes. By default this directory is created in the folder where the pipeline is executed, but we recommend to set it to a scratch or tmp folder. To set it in the Configuration file, the following code line should be added and edited:

```nextflow
workDir = '/absolute/path/to/work/dir'
```
### Enabling Singularity container 
[Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) is a container platform. It allows to run the nextflow pipeline with a pre-determined environment that ensures reproducibility. We created a Docker image with the necessary dependencies for MINI-AC to run, without need for the user to install any of them. We strongly recommend to always run MINI-AC using the specified singularity container. For that, the following code needs to be included in the Configuration file:

```nextflow
process.container = "vibpsb/mini-ac:latest"
singularity {
    enabled = true
    cacheDir = "singularity_cache"
    autoMounts = true
}
```
Sometimes the temporary file used by Singularity is not in the same root path as in the pipeline, which can cause Singularity to struggle to find it. In this case, add the ```runOptions``` line below with the absolute path to the tmp folder. To know the absolute path to the tmp folder in linux execute in the command line ```echo $TMPDIR```. Then add it as shown below.

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

[The executor of the pipeline](https://www.nextflow.io/docs/latest/executor.html) is the system where the pipeline processes run and supervises its execution. It can be a computer, a cluster resource manager, or the cloud. In the Configuration file it can be specified what is the executor of the pipeline. To execute it in a normal computer locally, the code below should be added the configuration file:

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
