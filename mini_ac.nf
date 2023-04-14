nextflow.enable.dsl = 2

include { genome_wide_miniac } from './workflows/miniac_gw'
include { locus_based_miniac } from './workflows/miniac_lb'

workflow MINIAC {
    if (params.mode == "genome_wide" || params.mode == "locus_based" ) {

        params.OBO_file = "$projectDir/data/ontologies/go.obo"

        if (params.mode == "genome_wide") {

            if (params.species == "maize" || params.species == "arabidopsis" ) {

                if (params.species == "maize") {

                    params.MotMapsFile_gw = "$projectDir/data/zma/zma_genome_wide_motif_mappings.bed"
                    params.Non_cod_genome = "$projectDir/data/zma/zma_noncod_merged.bed"
                    params.Faix_file = "$projectDir/data/zma/zma.fasta.fai"
                    params.Motif_tf_file = "$projectDir/data/zma/zma_motif_TF_file.txt"
                    params.Genes_coords = "$projectDir/data/zma/zma_genes_coords_sorted.bed"
                    params.Feature_file = "$projectDir/data/zma/zma_go_gene_file.txt"
                    params.TF_fam_file = "$projectDir/data/zma/zma_TF_family_file.txt"
                    params.Genes_metadata = "$projectDir/data/zma/maize_gene_metadata_file.txt"
                    params.P_val = 0.1
                    
                    genome_wide_miniac(params.OutDir, params.ACR_dir, params.Filter_set_genes, params.Set_genes_dir, params.One_filtering_set, params.DE_genes, params.DE_genes_dir, params.One_DE_set, params.P_val, params.Bps_intersect, params.Second_gene_annot, params.Second_gene_dist, params.MotMapsFile_gw, params.Non_cod_genome, params.Faix_file, params.Motif_tf_file, params.Genes_coords, params.Feature_file, params.OBO_file, params.TF_fam_file, params.Genes_metadata)

                }
                else {

                    params.MotMapsFile_gw = "$projectDir/data/ath/ath_genome_wide_motif_mappings.bed"
                    params.Non_cod_genome = "$projectDir/data/ath/ath_noncod_merged.bed"
                    params.Faix_file = "$projectDir/data/ath/ath.fasta.fai"
                    params.Motif_tf_file = "$projectDir/data/ath/ath_motif_TF_file.txt"
                    params.Genes_coords = "$projectDir/data/ath/ath_genes_coords_sorted.bed"
                    params.Feature_file = "$projectDir/data/ath/ath_go_gene_file.txt"
                    params.TF_fam_file = "$projectDir/data/ath/ath_TF_family_file.txt"
                    params.Genes_metadata = "$projectDir/data/ath/arabidopsis_gene_metadata_file.txt"
                    params.P_val = 0.1

                    genome_wide_miniac(params.OutDir, params.ACR_dir, params.Filter_set_genes, params.Set_genes_dir, params.One_filtering_set, params.DE_genes, params.DE_genes_dir, params.One_DE_set, params.P_val, params.Bps_intersect, params.Second_gene_annot, params.Second_gene_dist, params.MotMapsFile_gw, params.Non_cod_genome, params.Faix_file, params.Motif_tf_file, params.Genes_coords, params.Feature_file, params.OBO_file, params.TF_fam_file, params.Genes_metadata)
                }
            }
            else {
                exit 1, "MINI-AC is only implemented for 2 species: maize and arabidopsis. Instead it got '${params.species}' "
                }

        }
        else {

            if (params.species == "maize" || params.species == "arabidopsis" ) {

                if (params.species == "maize") {

                    params.MotMapsFile_lb = "$projectDir/data/zma/zma_locus_based_motif_mappings_5kbup_1kbdown.bed"
                    params.Promoter_file = "$projectDir/data/zma/zma_promoter_5kbup_1kbdown_sorted.bed"
                    params.Faix_file = "$projectDir/data/zma/zma.fasta.fai"
                    params.Motif_tf_file = "$projectDir/data/zma/zma_motif_TF_file.txt"
                    params.Feature_file = "$projectDir/data/zma/zma_go_gene_file.txt"
                    params.TF_fam_file = "$projectDir/data/zma/zma_TF_family_file.txt"
                    params.Genes_metadata = "$projectDir/data/zma/maize_gene_metadata_file.txt"
                    params.P_val = 0.01
                    
                    locus_based_miniac(params.OutDir, params.ACR_dir, params.Filter_set_genes, params.Set_genes_dir, params.One_filtering_set, params.DE_genes, params.DE_genes_dir, params.One_DE_set, params.P_val, params.Bps_intersect, params.MotMapsFile_lb, params.Promoter_file, params.Faix_file, params.Motif_tf_file, params.Feature_file, params.OBO_file, params.TF_fam_file, params.Genes_metadata)


                }
                else {

                    params.MotMapsFile_lb = "$projectDir/data/ath/ath_locus_based_motif_mappings_5kbup_1kbdown.bed"
                    params.Promoter_file = "$projectDir/data/ath/ath_promoter_5kbup_1kbdown_sorted.bed"
                    params.Faix_file = "$projectDir/data/ath/ath.fasta.fai"
                    params.Motif_tf_file = "$projectDir/data/ath/ath_motif_TF_file.txt"
                    params.Feature_file = "$projectDir/data/ath/ath_go_gene_file.txt"
                    params.TF_fam_file = "$projectDir/data/ath/ath_TF_family_file.txt"
                    params.Genes_metadata = "$projectDir/data/ath/arabidopsis_gene_metadata_file.txt"
                    params.P_val = 0.01

                    locus_based_miniac(params.OutDir, params.ACR_dir, params.Filter_set_genes, params.Set_genes_dir, params.One_filtering_set, params.DE_genes, params.DE_genes_dir, params.One_DE_set, params.P_val, params.Bps_intersect, params.MotMapsFile_lb, params.Promoter_file, params.Faix_file, params.Motif_tf_file, params.Feature_file, params.OBO_file, params.TF_fam_file, params.Genes_metadata)


                }
            }
            else {
                exit 1, "MINI-AC is only implemented for two species: maize and arabidopsis. Instead it got '${params.species}' "
                }

        }
    }
    else {
        exit 1, "MINI-AC can only be run on two modes: genome_wide or locus_based. Instead got '${params.mode}' "
        }
}


workflow {
    MINIAC()
}