nextflow.enable.dsl = 2
include { genome_wide_miniac } from './workflows/miniac_gw'
include { GENE_LIST_TO_ICRES_BED } from './modules/gene_list_to_iCREs_BED/main.nf'

workflow MINIAC_ICRES {

    params.OBO_file = "$projectDir/data/ontologies/go.obo"
    params.Shuffle_count = 1000
    params.Shuffle_seed = -1
    params.Csv_output = false

    // define species id used for data subfolder and data file prefix
    def species
    switch(params.species) {
        case "maize_v4":
            species = "zma_v4"
            break
        case "maize_v5":
            species = "zma_v5"
            break
        default:
            exit 1, "MINI-AC iCREs can only be run for the species 'maize_v4' and 'maize_v5' (not Arabidopsis). Instead it got '${params.species}'."
    }

    // define iCREs set used for iCREs coordinates file retrieval
    def icres_set
    switch(params.icres_set) {
        case "all":
            icres_set = "all_icres"
            break
        case "maxf1":
            icres_set = "maxf1_icres"
            break
        default:
            exit 1, "MINI-AC iCREs can only be run for two sets: all and maxf1. Instead it got '${params.icres_set}'."
    }

    // set input data parameters for the exection of MINI-AC genome-wide mode
    params.Faix_file = "$projectDir/data/${species}/${species}.fasta.fai"
    params.Motif_tf_file = "$projectDir/data/${species}/${species}_motif_TF_file.txt"
    params.Feature_file = "$projectDir/data/${species}/${species}_go_gene_file.txt"
    params.TF_fam_file = "$projectDir/data/${species}/${species}_TF_family_file.txt"
    params.Genes_metadata = "$projectDir/data/${species}/${species}_gene_metadata_file.txt"

    params.icres_set_file = "$projectDir/data/icres/${icres_set}_${species}.bed"

    Gene_sets = Channel.fromPath("${params.Gene_list_dir}/*.txt").ifEmpty { error "No *.txt files could be found in the specified gene sets directory ${params.Gene_list_dir}" }

    script_get_icres_bed = "${projectDir}/bin/geneList2iCREs.py"


    ACR_files = GENE_LIST_TO_ICRES_BED(params.icres_set_file, Gene_sets, script_get_icres_bed)
        
    params.MotMapsFile = "$projectDir/data/${species}/${species}_genome_wide_motif_mappings.bed"
    params.Non_cod_genome = "$projectDir/data/${species}/${species}_noncod_merged.bed"
    params.Genes_coords = "$projectDir/data/${species}/${species}_genes_coords_sorted.bed"

    params.P_val = 0.1

    genome_wide_miniac(params, ACR_files)

}

workflow {
    MINIAC_ICRES()
}