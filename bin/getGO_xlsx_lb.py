### Importing modules ###

import pandas as pd
import sys
import argparse

### Function to read input arguments ###

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to get GO enrichment integrative MINI-AC output', 
                        conflict_handler = 'resolve')

    parser.add_argument('go_enrichment_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'MINI-AC GO enrichment raw output file')

    parser.add_argument('tf_fam_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'TF-family file')

    parser.add_argument('info_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'genes information file')

    parser.add_argument('output_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'Name of the output file')

    parser.add_argument('-de', '--de_table_file', nargs = 1, type = str,
                        default = None, help = '',
                        metavar = 'Table with differential expression data')

    parser.add_argument('-ex', '--expressed_genes_file', nargs = 1, type = str,
                        default = None, help = '',
                        metavar = 'List of genes expressed in biological context of experiment')
    
    parser.add_argument('-out', '--output_type', nargs = 1, type = str,
                        default = 'excel', help = '',
                        metavar = 'output type (excel or csv, defaults to excel)')
    
    args = parser.parse_args()

    return args

### Reading inputs ###

args = parseArgs()

go_enrichment_file = args.go_enrichment_file[0]
tf_fam_file = args.tf_fam_file[0]
info_file = args.info_file[0]
output_file = args.output_file[0]
DE_table_file = args.de_table_file
expressed_genes_file = args.expressed_genes_file
output_type = args.output_type

info_df = pd.read_csv(info_file)
tf_fam = pd.read_csv(tf_fam_file)
info_df_columns = list(info_df.columns)[1:]

### Reading expression data tables ###

if expressed_genes_file:
    exp_genes = set()
    with open(expressed_genes_file[0], 'r') as f:
        for line in f:
            rec = line.strip().split("\t")
            gene = rec[0]
            exp_genes.add(gene)

if DE_table_file:
    DE_table = pd.read_csv(DE_table_file[0], sep = "\t")

    for col_num in range(len(DE_table.columns.values)):
        if col_num == 0:
            DE_table.columns.values[col_num] = "gene_id"
        else:
            DE_table.columns.values[col_num] = DE_table.columns.values[col_num] + " DE"

    DE_genes = set(DE_table.gene_id)
    DE_table_columns = list(DE_table.columns)[1:]

### Reading and processing GO enrichment data ###

GO_info = {}

with open(go_enrichment_file, 'r') as f:
    for line in f:
        rec = line.strip().split("\t")
        if line.startswith("###"):
            GO_info = None
            break
        if line.startswith("#set_id"):
            continue
        tf, go_term = rec[0], rec[-1]
        if (tf, go_term) not in GO_info:
            GO_info[(tf, go_term)] = rec[1:-1]
            
if not GO_info:
    empty_table = pd.DataFrame(["### This dataset did not yield any GO enrichment"])
    with pd.ExcelWriter(output_file) as writer:
        empty_table.to_excel(writer, index = False, header = False)
    sys.exit()
    
### Integrating data ###

if GO_info:                
    go_df = pd.DataFrame().from_dict(GO_info, orient = 'index')
    go_df.index = pd.MultiIndex.from_tuples(go_df.index, names=['gene_id', 'go_term'])
    go_df = go_df.reset_index()
    go_df.columns = ['Gene ID', 'GO term', 'GO term ID', 'p-value', 'q-value', 'Enrichment fold', 'Number of target genes', 'Number of genes annotated to GO', 'Number of target genes annotated with GO', 'Target genes']
    go_df = go_df.merge(tf_fam, how = 'left', left_on = 'Gene ID', right_on = 'gene_id').merge(info_df, how = 'left', left_on = 'Gene ID', right_on = 'gene_id')

    if DE_table_file:
    
        if expressed_genes_file:
            go_df['Expressed gene'] = go_df.loc[:, 'Gene ID'].isin(exp_genes)
            info_df_columns = info_df_columns + ['Expressed gene']

        go_df = go_df.merge(DE_table, how = 'left', left_on = 'Gene ID', right_on = 'gene_id')
        
        go_df = go_df[['Gene ID', 'GO term', 'q-value', 'family'] + info_df_columns + DE_table_columns + ['GO term ID', 'p-value', 'Enrichment fold', 'Number of target genes', 'Number of genes annotated to GO', 'Number of target genes annotated with GO', 'Target genes']]
        
    if not DE_table_file:
        
        if expressed_genes_file:
            go_df['Expressed gene'] = go_df.loc[:, 'Gene ID'].isin(exp_genes)
            info_df_columns = info_df_columns + ['Expressed gene']
        
        go_df = go_df[['Gene ID', 'GO term', 'q-value', 'family'] + info_df_columns + ['GO term ID', 'p-value', 'Enrichment fold', 'Number of target genes', 'Number of genes annotated to GO', 'Number of target genes annotated with GO', 'Target genes']]
        

go_df = go_df.rename(columns = {'gene_name': 'Gene name', 'description': 'Gene description', 'family': 'TF family', 'ath_gene_id': 'Arabidopsis gene ID', 'ath_alias': 'Arabidopsis gene name', 'gene_name_figs': 'Gene name (+ alias)'})

### Writing output file ###

if (output_type == 'csv'):
    go_df.to_csv(output_file, index = False)
else:
    with pd.ExcelWriter(output_file) as writer:
        go_df.to_excel(writer, index = False)
