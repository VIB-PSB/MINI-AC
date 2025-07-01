# FAQ

## Q: MINI-AC failed, how can I fix it?
A: 
* Check the [config file](/docs/pipeline_configuration.md):
  * Did you specify the correct [executor](https://www.nextflow.io/docs/latest/executor.html) (e.g. SGE, SLURM, ...)? Cluster-related options (i.e., all the lines starting with `clusterOptions`) should also be adapted to match the options of the selected executor.
  * Did you [specify to Singularity the path to the temporary directory](https://docs.sylabs.io/guides/3.5/user-guide/bind_paths_and_mounts.html)? It can be done by adjusting the parameter ```runOptions``` of singularity in Nextflow to ```--bind /absolute/path/to/tmp/folder```. To know the absolute path to the tmp folder in linux execute in the command line ```echo $TMPDIR```
 
     ```nextflow
    process.container = "vibpsb/mini-ac:latest"
    singularity {
        enabled = true
        cacheDir = "singularity_cache"
        autoMounts = true
        runOptions = "--bind /absolute/path/to/tmp/folder"
    }
    ```
* Check the Nextflow output to see which specific process failed. Try to increase the memory allocated to that process in the config file. For example, if the process "getStats" failed, change `-l h_vmem = 10G` to `-l h_vmem = 20G` for the getStats process in the config file.
