nextflow.enable.dsl=2

process get_ACR_shufflings{
    tag "${acr_input_file.baseName}"

    input:
    path acr_input_file
    path faix
    path non_cod_genome
    val count
    val seed

    output:
    path "${acr_input_file.baseName}_allshuff_sorted.bed", emit: shufflings
    path "${acr_input_file.baseName}_3col.bed", emit: acr_input
    

    script:

    seed_arg = ''
    if (seed) {
        seed_arg = "-seed $seed"
    }

    """
    cut -f1,2,3 ${acr_input_file} | sort-bed - | bedops -m - > ${acr_input_file.baseName}_3col.bed
    
    sed "s/\$/\treal_ints/" ${acr_input_file.baseName}_3col.bed > ${acr_input_file.baseName}_allshuff.bed

    bedtools subtract -a $non_cod_genome -b ${acr_input_file.baseName}_3col.bed > shuffling_space.bed

    for i in {1..$count};
    do
       bedtools shuffle $seed_arg -i ${acr_input_file.baseName}_3col.bed -g $faix -incl shuffling_space.bed | sed "s/\$/\t\$i/" >> ${acr_input_file.baseName}_allshuff.bed
    done
   
    sort-bed --max-mem 2G ${acr_input_file.baseName}_allshuff.bed > ${acr_input_file.baseName}_allshuff_sorted.bed

    """
}

process getStats {
   publishDir "${outDir}/enrichment_stats", pattern: "*_miniac_stats.txt", mode: 'copy'
   tag "${shuffled_real_acr.baseName}"

   input:
   tuple path(shuffled_real_acr), path(gw_motmaps)
   path get_stats_script
   path mot_tf
   val outDir
   val shuffle_count

   output:
   path "${shuffled_real_acr.baseName}_motif_int_data.txt", emit: raw_stats
   path "${shuffled_real_acr.baseName}_miniac_stats.txt", emit: proc_stats

   script:

   """
   #bedmap --delim '\t' --fraction-ref 0.5 --echo --echo-map-id --skip-unmapped $gw_motmaps $shuffled_real_acr | awk '{n=split(\$0,a,";"); split(a[1],b,"\t"); print b[4], b[5]} n>1 {for (i=2;i<=n;i++) print b[4], a[i]}' | sort | uniq -c > ${shuffled_real_acr.baseName}_motif_int_data.txt

   #bedmap --delim '\t' --fraction-either 0.5 --echo --echo-map-id --skip-unmapped $gw_motmaps $shuffled_real_acr | awk '{n=split(\$0,a,";"); split(a[1],b,"\t"); print b[4], b[5]} n>1 {for (i=2;i<=n;i++) print b[4], a[i]}' | sort | uniq -c > ${shuffled_real_acr.baseName}_motif_int_data.txt

   #bedmap --delim '\t' --bp-ovr 5 --echo --echo-map-id --skip-unmapped $gw_motmaps $shuffled_real_acr | awk '{n=split(\$0,a,";"); split(a[1],b,"\t"); print b[4], b[5]} n>1 {for (i=2;i<=n;i++) print b[4], a[i]}' | sort | uniq -c > ${shuffled_real_acr.baseName}_motif_int_data.txt

   bedmap --delim '\t' --echo --echo-map-id --skip-unmapped $gw_motmaps $shuffled_real_acr | awk '{n=split(\$0,a,";"); split(a[1],b,"\t"); print b[4], b[5]} n>1 {for (i=2;i<=n;i++) print b[4], a[i]}' | sort | uniq -c > ${shuffled_real_acr.baseName}_motif_int_data.txt

   num_peaks=\$(grep -c "real_ints" $shuffled_real_acr)

   OMP_NUM_THREADS=1 python3 $get_stats_script ${shuffled_real_acr.baseName}_motif_int_data.txt \$num_peaks $mot_tf ${shuffled_real_acr.baseName}_miniac_stats.txt $shuffle_count

   """

}

process getStats_bps {
   publishDir "${outDir}/enrichment_stats", pattern: "*_miniac_stats.txt", mode: 'copy'
   tag "${shuffled_real_acr.baseName}"

   input:
   tuple path(shuffled_real_acr), path(gw_motmaps)
   path get_stats_script_bps
   path mot_tf
   val outDir
   val shuffle_count

   output:
   path "${shuffled_real_acr.baseName}_int_bps_data.bed", emit: raw_stats
   path "${shuffled_real_acr.baseName}_miniac_stats.txt", emit: proc_stats

   script:

   """
   #bedmap --delim '\t' --echo --echo-overlap-size --echo-map-id --skip-unmapped $gw_motmaps $shuffled_real_acr > ${shuffled_real_acr.baseName}_int_bps_data.bed

   #bedmap --delim '\t' --fraction-either 0.5 --echo --echo-overlap-size --echo-map-id --skip-unmapped $gw_motmaps $shuffled_real_acr > ${shuffled_real_acr.baseName}_int_bps_data.bed

   #bedmap --delim '\t' --bp-ovr 5 --echo --echo-overlap-size --echo-map-id --skip-unmapped  $gw_motmaps $shuffled_real_acr > ${shuffled_real_acr.baseName}_int_bps_data.bed

   bedmap --delim '\t' --echo --echo-overlap-size --echo-map-id --skip-unmapped  $gw_motmaps $shuffled_real_acr > ${shuffled_real_acr.baseName}_int_bps_data.bed

   num_peaks=\$(grep -c "real_ints" $shuffled_real_acr)

   OMP_NUM_THREADS=1 python3 $get_stats_script_bps ${shuffled_real_acr.baseName}_int_bps_data.bed \$num_peaks $mot_tf ${shuffled_real_acr.baseName}_miniac_stats.txt $shuffle_count
   """

}

process getNetwork {
   publishDir "${outDir}/networks", mode: 'copy'
   tag "${acr_file_name}"


   input:
   tuple val(acr_file_name), path(enr_stats), path(acr_file), path(motmaps)
   path zma_genes_coords
   path script_get_net
   path mot_tf
   val second_gene_dist
   val second_gene_annot
   val pval
   val outDir

   output:
   path "${acr_file_name}_network.txt"

   script:
   ann_second_gene = second_gene_dist == true ? "--second_gene_dist ${second_gene_dist}" : ''
   num_genes_annot = second_gene_annot == true ? "2" : '1'

   """
   #bedops -e 50% $motmaps $acr_file | bedtools closest -a stdin -b $zma_genes_coords -t all -d -D b -k $num_genes_annot > ${acr_file_name}_peak_ann.bed

   #bedmap --delim '\t' --fraction-either 0.5 --echo --skip-unmapped $motmaps $acr_file | bedtools closest -a stdin -b $zma_genes_coords -t all -d -D b -k $num_genes_annot > ${acr_file_name}_peak_ann.bed

   #bedmap --delim '\t' --bp-ovr 5 --echo --skip-unmapped $motmaps $acr_file | bedtools closest -a stdin -b $zma_genes_coords -t all -d -D b -k $num_genes_annot > ${acr_file_name}_peak_ann.bed

   bedmap --delim '\t' --echo --skip-unmapped $motmaps $acr_file | bedtools closest -a stdin -b $zma_genes_coords -t all -d -D b -k $num_genes_annot > ${acr_file_name}_peak_ann.bed

   OMP_NUM_THREADS=1 python3 $script_get_net ${acr_file_name}_peak_ann.bed $enr_stats $mot_tf $pval ${acr_file_name}_network_dups.txt $ann_second_gene

   sort -u ${acr_file_name}_network_dups.txt > ${acr_file_name}_network.txt

   """

}

process filterSetOfGenes{
    publishDir "${outDir}/networks", mode: 'copy'
    tag "${network.baseName}"

    input:
    path filterScript
    tuple path(network), path(gene_set)
    val outDir

    output:
    path "${network.baseName}_filtered.txt"

    script:

    """
    python3 $filterScript $network $gene_set ${network.baseName}_filtered.txt
    
    """
}


process GOenrichment{
    publishDir "${outDir}/GO_enrichment", mode: 'copy'
    tag "${acr_file_name}"

    input:
    tuple val(acr_file_name), path(set_file)
    path feature_file
    path enricher
    path reduce_go
    path add_go_names
    path obo_file
    path filter_reduced_script
    val outDir

    output:
    path "${acr_file_name}_GO_enrichment.txt"

    script:

    """
    cut -f1,2 $set_file | tail -n +2 > ${acr_file_name}_2col.txt

    ./$enricher $feature_file ${acr_file_name}_2col.txt -p > ${acr_file_name}_GO_enrichment_rawOut.txt

    count=\$(wc -l ${acr_file_name}_GO_enrichment_rawOut.txt | cut -d" " -f1)

    if [[ \$count -gt 7 ]]
    
    then

        grep -v '#' ${acr_file_name}_GO_enrichment_rawOut.txt > ${acr_file_name}_GO_enrichment_noheader.txt

        python $add_go_names ${acr_file_name}_GO_enrichment_noheader.txt $obo_file | sort -g -k4,4 > ${acr_file_name}_GO_enrichment_go_names.txt

        awk 'BEGIN {FS="\t"; OFS="\t"}; {print \$2, \$1}' ${acr_file_name}_GO_enrichment_go_names.txt > go_gene.txt

         python $reduce_go go_gene.txt $obo_file > go_gene_reduced.txt

        python3 $filter_reduced_script go_gene_reduced.txt ${acr_file_name}_GO_enrichment_go_names.txt ${acr_file_name}_GO_enrichment_reduced.txt

        sort -g -k4,4 ${acr_file_name}_GO_enrichment_reduced.txt > ${acr_file_name}_GO_enrichment.txt

    else

        echo "### The network corresponding to this file did not yield GO enrichment" > ${acr_file_name}_GO_enrichment.txt
        echo "### The network corresponding to this file did not yield GO enrichment" > ${acr_file_name}_GO_enrichment_go_names.txt

    fi

    mv ${acr_file_name}_GO_enrichment_go_names.txt ${acr_file_name}_GO_enrichment_nonred.txt

    """
}

process getIntegrativeOutputs {
   publishDir "${outDir}/integrative_outputs", mode: 'copy'
   tag "${acr_file_name}"

   input:
   tuple val(acr_file_name), path(enr_stats), path(exp_genes), path(network), path(go_enr), path(de_genes)
   path mot_tf
   path tf_fam
   path info_file
   val pval
   val exp_genes_bool
   val de_genes_bool
   path script_tfs
   path script_motifs
   path script_go
   path script_net
   val outDir
   val csvOutput            // set to true for CSV instead of Excel output


   output:
   path "${acr_file_name}_TF_centric.xlsx"
   path "${acr_file_name}_motif_centric.xlsx"
   path "${acr_file_name}_GO_enrichment.xlsx"
   path "${acr_file_name}_functional_network.txt"
   path "${acr_file_name}_node_attributes.txt"


   script:
   de_genes_table = de_genes_bool == true ? "-de $de_genes" : ''
   exp_genes_table = exp_genes_bool == true ? "-ex $exp_genes" : ''
   csv_output_flag = csvOutput ? '-csv' : ''

   """
   OMP_NUM_THREADS=1 python3 $script_tfs $enr_stats $network $go_enr $mot_tf $tf_fam $info_file $pval ${acr_file_name}_TF_centric.xlsx $de_genes_table $exp_genes_table $csv_output_flag

   OMP_NUM_THREADS=1 python3 $script_motifs $enr_stats $mot_tf $tf_fam $info_file $pval ${acr_file_name}_motif_centric.xlsx $de_genes_table $exp_genes_table $csv_output_flag

   OMP_NUM_THREADS=1 python3 $script_go $go_enr $tf_fam $info_file ${acr_file_name}_GO_enrichment.xlsx $de_genes_table $exp_genes_table $csv_output_flag

   OMP_NUM_THREADS=1 python3 $script_net $enr_stats $network $go_enr $mot_tf $tf_fam $info_file $pval ${acr_file_name}_functional_network.txt ${acr_file_name}_node_attributes.txt $de_genes_table $exp_genes_table

   """

}

workflow genome_wide_miniac {
    take:
    OutDir
    ACR_dir
    Filter_set_genes
    Set_genes_dir
    One_filtering_set
    DE_genes
    DE_genes_dir
    One_DE_set
    P_val
    Bps_intersect
    Second_gene_annot
    Second_gene_dist
    MotMapsFile_gw
    Non_cod_genome
    Faix_file
    Motif_tf_file
    Genes_coords
    Feature_file
    OBO_file
    TF_fam_file
    Genes_metadata
    Shuffle_count
    Shuffle_seed
    Csv_output

    main:
    
    if (!file(MotMapsFile_gw).exists()) { error "Please make sure that you downloaded the motif mapping files as described in the documentation." }
    
    ACR_files = Channel.fromPath("${ACR_dir}/*.bed").ifEmpty { error "No *.bed files could be found in the specified ACR directory ${ACR_dir}" }
    
    get_ACR_shufflings(ACR_files, Faix_file, Non_cod_genome, Shuffle_count, Shuffle_seed)    

    acr_shufflings_ch = get_ACR_shufflings.out.shufflings

    parsed_acr = get_ACR_shufflings.out.acr_input
                                        .map {n -> [n.baseName.split("_")[0..-2].join("_"), n]}

    motmaps_ch = Channel.fromPath(MotMapsFile_gw)

    input_stats = acr_shufflings_ch.combine(motmaps_ch)

    if (Bps_intersect == false) {

        script_proc_stats = "${projectDir}/bin/processStats_gw.py"

        getStats(input_stats, script_proc_stats, Motif_tf_file, OutDir, Shuffle_count)

        stats_ch = getStats.out.proc_stats
                    .map { n -> [n.BaseName.split("_")[0..-5].join("_"), n]}

                                            }

    else {
        
        script_proc_stats_bps = "${projectDir}/bin/processStats_bps_gw.py"

        getStats_bps(input_stats, script_proc_stats_bps, Motif_tf_file, OutDir, Shuffle_count)

        stats_ch = getStats_bps.out.proc_stats
                    .map { n -> [n.BaseName.split("_")[0..-5].join("_"), n]}

                                        }

    stats_acr_motmaps_ch = stats_ch.join(parsed_acr).combine(motmaps_ch)


    script_getnetwork = "${projectDir}/bin/getNetwork_gw.py"

    getNetwork(stats_acr_motmaps_ch, Genes_coords, script_getnetwork, Motif_tf_file, Second_gene_dist, Second_gene_annot, P_val, OutDir)

    networks = getNetwork.out

    if (Filter_set_genes == true) {

        filteringScript = "${projectDir}/bin/filterNetwork_gw.py"

        filt_set_files = Channel.fromPath("${Set_genes_dir}/*.txt")
                .ifEmpty { error "Cannot find any directory: ${Set_genes_dir}" }

        if (One_filtering_set == true) {

            networks_gene_set = networks.combine(filt_set_files)

            int_input = stats_ch.combine(filt_set_files)

                                                }

        else {
            
            filt_set_files_tup = filt_set_files.flatten()
                                .map { n -> [n.baseName.split("_")[0..-3].join("_"), n] }

            networks_tup = networks.flatten()
                                .map { n -> [n.baseName.split("_")[0..-2].join("_"), n] }

            networks_gene_set = networks_tup.join(filt_set_files_tup)
                                .map {n -> [ n[1], n[2] ]}
            
            networks_gene_set.ifEmpty { error "Expression and ACR files name don't match" }

            int_input = stats_ch.join(filt_set_files_tup)

                                            }

        filterSetOfGenes(filteringScript, networks_gene_set, OutDir)

        net_tuple = filterSetOfGenes.out
                    .map { n -> [n.BaseName.split("_")[0..-3].join("_"), n]}

                                        }
    else {
        
        int_input = stats_ch.combine(Channel.fromPath('NO_FILE'))

        net_tuple = getNetwork.out
                    .map { n -> [n.BaseName.split("_")[0..-2].join("_"), n]}

                                        }

    int_input = int_input.join(net_tuple)

    script_enricher = "${projectDir}/bin/enricherv2.4"
    script_reduce_go = "${projectDir}/bin/reduce_go.py"
    script_add_go_names = "${projectDir}/bin/add_go_names.py"
    script_filter_reduced = "${projectDir}/bin/FilterReduceGO.py"

    GOenrichment(net_tuple, Feature_file,
                script_enricher, script_reduce_go,
                script_add_go_names, OBO_file, script_filter_reduced, OutDir)

    go_tuple = GOenrichment.out
            .map { n -> [n.BaseName.split("_")[0..-3].join("_"), n]}

    int_input = int_input.join(go_tuple)

    if (DE_genes == true) {

        de_files = Channel.fromPath("${DE_genes_dir}/*.txt")
                                .ifEmpty { error "Cannot find any directory: ${Set_genes_dir}" }

        if (One_DE_set == true) {

            int_input = int_input.combine(de_files)

                                            }

        else {

            de_files_tup = de_files.flatten()
                                    .map { n -> [n.baseName.split("_")[0..-3].join("_"), n] }

            int_input = int_input.join(de_files_tup)
            
                    }
            }
    else {

        int_input = int_input.combine(Channel.fromPath('NO_FILE2'))

                                    }

    script_tf_file = "${projectDir}/bin/getTFCentricOutput_gw.py"
    script_motif_file = "${projectDir}/bin/getMotifCentricOutput_gw.py"
    script_go_file = "${projectDir}/bin/getGO_xlsx_gw.py"
    script_net_files = "${projectDir}/bin/getNetVisualizationOutput_gw.py"

    getIntegrativeOutputs(int_input, Motif_tf_file, TF_fam_file, Genes_metadata, P_val, Filter_set_genes, DE_genes,
                          script_tf_file, script_motif_file, script_go_file, script_net_files, OutDir, Csv_output)

        }
