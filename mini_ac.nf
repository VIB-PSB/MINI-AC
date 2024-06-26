nextflow.enable.dsl = 2
include { genome_wide_miniac } from './workflows/miniac_gw'
include { locus_based_miniac } from './workflows/miniac_lb'

workflow MINIAC {

    params.OBO_file = "$projectDir/data/ontologies/go.obo"
    params.Shuffle_count = 1000
    params.Shuffle_seed = -1
    params.Csv_output = false

    // define species id used for data subfolder and data file prefix
    def species
    switch(params.species) {
        case "arabidopsis":
            species = "ath"
            break
        case "maize_v4":
            species = "zma_v4"
            break
        case "maize_v5":
            species = "zma_v5"
            break
        default:
            exit 1, "MINI-AC can only be run for the species 'arabidopsis', 'maize_v4' and 'maize_v5'. Instead it got '${params.species}'."
    }

    // set input data parameters shared between genome-wide and locus-based modes
    params.Faix_file = "$projectDir/data/${species}/${species}.fasta.fai"
    params.Motif_tf_file = "$projectDir/data/${species}/${species}_motif_TF_file.txt"
    params.Feature_file = "$projectDir/data/${species}/${species}_go_gene_file.txt"
    params.TF_fam_file = "$projectDir/data/${species}/${species}_TF_family_file.txt"
    params.Genes_metadata = "$projectDir/data/${species}/${species}_gene_metadata_file.txt"
    ACR_files = Channel.fromPath("${params.ACR_dir}/*.bed").ifEmpty { error "No *.bed files could be found in the specified ACR directory ${params.ACR_dir}" }

    if (params.mode == "genome_wide") {
        
        params.MotMapsFile = "$projectDir/data/${species}/${species}_genome_wide_motif_mappings.bed"
        params.Non_cod_genome = "$projectDir/data/${species}/${species}_noncod_merged.bed"
        params.Genes_coords = "$projectDir/data/${species}/${species}_genes_coords_sorted.bed"

        params.P_val = 0.1

        genome_wide_miniac(params, ACR_files)
    
    } else if (params.mode == "locus_based") {

        params.MotMapsFile = "$projectDir/data/${species}/${species}_locus_based_motif_mappings_5kbup_1kbdown.bed"
        params.Promoter_file = "$projectDir/data/${species}/${species}_promoter_5kbup_1kbdown_sorted.bed"

        params.P_val = 0.01
        
        locus_based_miniac(params, ACR_files)
    
    } else {
        exit 1, "MINI-AC can only be run using the modes 'genome_wide' or 'locus_based'. Instead it got '${params.mode}'."
    }
}

workflow {
    MINIAC()
}