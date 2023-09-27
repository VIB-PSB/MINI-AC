nextflow.enable.dsl = 2

include { genome_wide_miniac } from './workflows/miniac_gw'
include { locus_based_miniac } from './workflows/miniac_lb'

workflow MINIAC {

    params.OBO_file = "$projectDir/data/ontologies/go.obo"
    params.Shuffle_count = 1000
    params.Shuffle_seed = -1
    params.Csv_output = false

    if (params.mode == "genome_wide" && params.species == "maize_v4") {

        params.MotMapsFile_gw = "$projectDir/data/zma_v4/zma_v4_genome_wide_motif_mappings.bed"
        params.Non_cod_genome = "$projectDir/data/zma_v4/zma_v4_noncod_merged.bed"
        params.Faix_file = "$projectDir/data/zma_v4/zma_v4.fasta.fai"
        params.Motif_tf_file = "$projectDir/data/zma_v4/zma_v4_motif_TF_file.txt"
        params.Genes_coords = "$projectDir/data/zma_v4/zma_v4_genes_coords_sorted.bed"
        params.Feature_file = "$projectDir/data/zma_v4/zma_v4_go_gene_file.txt"
        params.TF_fam_file = "$projectDir/data/zma_v4/zma_v4_TF_family_file.txt"
        params.Genes_metadata = "$projectDir/data/zma_v4/maize_v4_gene_metadata_file.txt"
        params.P_val = 0.1
        
        genome_wide_miniac(params.OutDir, params.ACR_dir, params.Filter_set_genes, params.Set_genes_dir,
            params.One_filtering_set, params.DE_genes, params.DE_genes_dir, params.One_DE_set, params.P_val,
            params.Bps_intersect, params.Second_gene_annot, params.Second_gene_dist, params.MotMapsFile_gw,
            params.Non_cod_genome, params.Faix_file, params.Motif_tf_file, params.Genes_coords, params.Feature_file,
            params.OBO_file, params.TF_fam_file, params.Genes_metadata, params.Shuffle_count, params.Shuffle_seed,
            params.Csv_output)
    }

    else if (params.mode == "genome_wide" && params.species == "maize_v5") {

        params.MotMapsFile_gw = "$projectDir/data/zma_v5/zma_v5_genome_wide_motif_mappings.bed"
        params.Non_cod_genome = "$projectDir/data/zma_v5/zma_v5_noncod_merged.bed"
        params.Faix_file = "$projectDir/data/zma_v5/zma_v5.fasta.fai"
        params.Motif_tf_file = "$projectDir/data/zma_v5/zma_v5_motif_TF_file.txt"
        params.Genes_coords = "$projectDir/data/zma_v5/zma_v5_genes_coords_sorted.bed"
        params.Feature_file = "$projectDir/data/zma_v5/zma_v5_go_gene_file.txt"
        params.TF_fam_file = "$projectDir/data/zma_v5/zma_v5_TF_family_file.txt"
        params.Genes_metadata = "$projectDir/data/zma_v5/maize_v5_gene_metadata_file.txt"
        params.P_val = 0.1
        
        genome_wide_miniac(params.OutDir, params.ACR_dir, params.Filter_set_genes, params.Set_genes_dir,
            params.One_filtering_set, params.DE_genes, params.DE_genes_dir, params.One_DE_set, params.P_val,
            params.Bps_intersect, params.Second_gene_annot, params.Second_gene_dist, params.MotMapsFile_gw,
            params.Non_cod_genome, params.Faix_file, params.Motif_tf_file, params.Genes_coords, params.Feature_file,
            params.OBO_file, params.TF_fam_file, params.Genes_metadata, params.Shuffle_count, params.Shuffle_seed,
            params.Csv_output)
    }

    else if (params.mode == "genome_wide" && params.species == "arabidopsis") {

        params.MotMapsFile_gw = "$projectDir/data/ath/ath_genome_wide_motif_mappings.bed"
        params.Non_cod_genome = "$projectDir/data/ath/ath_noncod_merged.bed"
        params.Faix_file = "$projectDir/data/ath/ath.fasta.fai"
        params.Motif_tf_file = "$projectDir/data/ath/ath_motif_TF_file.txt"
        params.Genes_coords = "$projectDir/data/ath/ath_genes_coords_sorted.bed"
        params.Feature_file = "$projectDir/data/ath/ath_go_gene_file.txt"
        params.TF_fam_file = "$projectDir/data/ath/ath_TF_family_file.txt"
        params.Genes_metadata = "$projectDir/data/ath/arabidopsis_gene_metadata_file.txt"
        params.P_val = 0.1

        genome_wide_miniac(params.OutDir, params.ACR_dir, params.Filter_set_genes, params.Set_genes_dir,
            params.One_filtering_set, params.DE_genes, params.DE_genes_dir, params.One_DE_set, params.P_val,
            params.Bps_intersect, params.Second_gene_annot, params.Second_gene_dist, params.MotMapsFile_gw,
            params.Non_cod_genome, params.Faix_file, params.Motif_tf_file, params.Genes_coords, params.Feature_file,
            params.OBO_file, params.TF_fam_file, params.Genes_metadata, params.Shuffle_count, params.Shuffle_seed,
            params.Csv_output)

    }

    else if (params.mode == "locus_based" && params.species == "maize_v4") {

        params.MotMapsFile_lb = "$projectDir/data/zma_v4/zma_v4_locus_based_motif_mappings_5kbup_1kbdown.bed"
        params.Promoter_file = "$projectDir/data/zma_v4/zma_v4_promoter_5kbup_1kbdown_sorted.bed"
        params.Faix_file = "$projectDir/data/zma_v4/zma_v4.fasta.fai"
        params.Motif_tf_file = "$projectDir/data/zma_v4/zma_v4_motif_TF_file.txt"
        params.Feature_file = "$projectDir/data/zma_v4/zma_v4_go_gene_file.txt"
        params.TF_fam_file = "$projectDir/data/zma_v4/zma_v4_TF_family_file.txt"
        params.Genes_metadata = "$projectDir/data/zma_v4/maize_v4_gene_metadata_file.txt"
        params.P_val = 0.01
        
        locus_based_miniac(params.OutDir, params.ACR_dir, params.Filter_set_genes, params.Set_genes_dir,
            params.One_filtering_set, params.DE_genes, params.DE_genes_dir, params.One_DE_set, params.P_val,
            params.Bps_intersect, params.MotMapsFile_lb, params.Promoter_file, params.Faix_file, params.Motif_tf_file,
            params.Feature_file, params.OBO_file, params.TF_fam_file, params.Genes_metadata, params.Shuffle_count, params.Shuffle_seed,
            params.Csv_output)

    }

    else if (params.mode == "locus_based" && params.species == "maize_v5") {

        params.MotMapsFile_lb = "$projectDir/data/zma_v5/zma_v5_locus_based_motif_mappings_5kbup_1kbdown.bed"
        params.Promoter_file = "$projectDir/data/zma_v5/zma_v5_promoter_5kbup_1kbdown_sorted.bed"
        params.Faix_file = "$projectDir/data/zma_v5/zma_v5.fasta.fai"
        params.Motif_tf_file = "$projectDir/data/zma_v5/zma_v5_motif_TF_file.txt"
        params.Feature_file = "$projectDir/data/zma_v5/zma_v5_go_gene_file.txt"
        params.TF_fam_file = "$projectDir/data/zma_v5/zma_v5_TF_family_file.txt"
        params.Genes_metadata = "$projectDir/data/zma_v5/maize_v5_gene_metadata_file.txt"
        params.P_val = 0.01
        
        locus_based_miniac(params.OutDir, params.ACR_dir, params.Filter_set_genes, params.Set_genes_dir,
            params.One_filtering_set, params.DE_genes, params.DE_genes_dir, params.One_DE_set, params.P_val,
            params.Bps_intersect, params.MotMapsFile_lb, params.Promoter_file, params.Faix_file, params.Motif_tf_file,
            params.Feature_file, params.OBO_file, params.TF_fam_file, params.Genes_metadata, params.Shuffle_count, params.Shuffle_seed,
            params.Csv_output)

    }

    else if (params.mode == "locus_based" && params.species == "arabidopsis") {

        params.MotMapsFile_lb = "$projectDir/data/ath/ath_locus_based_motif_mappings_5kbup_1kbdown.bed"
        params.Promoter_file = "$projectDir/data/ath/ath_promoter_5kbup_1kbdown_sorted.bed"
        params.Faix_file = "$projectDir/data/ath/ath.fasta.fai"
        params.Motif_tf_file = "$projectDir/data/ath/ath_motif_TF_file.txt"
        params.Feature_file = "$projectDir/data/ath/ath_go_gene_file.txt"
        params.TF_fam_file = "$projectDir/data/ath/ath_TF_family_file.txt"
        params.Genes_metadata = "$projectDir/data/ath/arabidopsis_gene_metadata_file.txt"
        params.P_val = 0.01

        locus_based_miniac(params.OutDir, params.ACR_dir, params.Filter_set_genes, params.Set_genes_dir,
            params.One_filtering_set, params.DE_genes, params.DE_genes_dir, params.One_DE_set, params.P_val,
            params.Bps_intersect, params.MotMapsFile_lb, params.Promoter_file, params.Faix_file, params.Motif_tf_file,
            params.Feature_file, params.OBO_file, params.TF_fam_file, params.Genes_metadata, params.Shuffle_count, params.Shuffle_seed,
            params.Csv_output)
    }

    else {
        exit 1, "MINI-AC can only be run using the modes 'genome_wide' and 'locus_based', and with the species 'arabidopsis', 'maize_v4' and 'maize_v5'. Instead it got '${params.species}' and '${params.mode}' "
        }
}


workflow {
    MINIAC()
}