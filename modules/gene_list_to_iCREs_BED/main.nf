process GENE_LIST_TO_ICRES_BED {
   tag "${gene_list.baseName}"


   input:
   path annotated_icres
   path gene_list
   path script_geneslist2bed

   output:
   path "${gene_list.baseName}_icres.bed"

   script:

   """
   python3 $script_geneslist2bed $annotated_icres $gene_list ${gene_list.baseName}_icres.bed

   """

}