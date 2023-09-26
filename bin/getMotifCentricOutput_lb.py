### Importing modules ###

import pandas as pd
import numpy as np
import math
import sys
import argparse

### Function to read input arguments ###

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to get motif-centric integrative MINI-AC output', 
                        conflict_handler = 'resolve')

    parser.add_argument('motif_enrichment_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'MINI-AC motif enrichment raw output file')

    parser.add_argument('mot_tf_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'Motif-TF file')

    parser.add_argument('tf_fam_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'TF-family file')

    parser.add_argument('info_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'genes information file')

    parser.add_argument('p_val', nargs = '?', type = float,
                        default = 0.01, help = '',
                        metavar = 'p-value for network inference')

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

motif_enrichment_file = args.motif_enrichment_file[0]
mot_tf_file = args.mot_tf_file[0]
tf_fam_file = args.tf_fam_file[0]
info_file = args.info_file[0]
pval = args.p_val
output_file = args.output_file[0]
DE_table_file = args.de_table_file
expressed_genes_file = args.expressed_genes_file
output_type = args.output_type

info_df = pd.read_csv(info_file)
mot_tf = pd.read_csv(mot_tf_file)
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

    lst = list(DE_table)
    DE_table[lst] = DE_table[lst].astype(str)

    DE_genes = set(DE_table.gene_id)
    DE_table_columns = list(DE_table.columns)
    DE_table_columns = DE_table_columns[1:]

### Reading and processing motif enrichment data ###

enr_stats = pd.read_csv(motif_enrichment_file, sep = "\t")
enr_stats.insert(7, 'pi_value', enr_stats.enr_fold * (-np.log10(enr_stats.adj_pval)))
enr_stats = enr_stats.sort_values(by = 'pi_value', ascending = False)
enr_stats["rank_pi_val"] = enr_stats.groupby("dataset")["pi_value"].rank("first", ascending = False)
enr_stats = enr_stats.merge(mot_tf, how = 'right', left_on = 'motif', right_on = 'motif_id')
enr_stats = enr_stats[enr_stats.adj_pval <= pval]

if enr_stats.empty:
    empty_table = pd.DataFrame(["### This dataset did not yield any motif enrichment"])
    with pd.ExcelWriter(output_file) as writer:
        empty_table.to_excel(writer, index = False, header = False)
    sys.exit()

for col in enr_stats.select_dtypes(include = ['float']).columns:
    if col == 'rank_pi_val':
        continue
    enr_stats[col] = enr_stats[col].apply(lambda x: x if(math.isnan(x)) else round(x,3))
enr_stats = enr_stats.astype({'real_int':'int'})
enr_stats = enr_stats.astype({'shuffled_int':'int'})
enr_stats = enr_stats.astype({'rank_pi_val':'int'})
enr_stats = enr_stats.merge(tf_fam, how = 'left', on = 'gene_id').drop('motif_id', axis = 1)

### Integrating everything ###

if expressed_genes_file:
    enr_stats['Any expressed gene'] = enr_stats.gene_id.isin(exp_genes)

    enr_stats = enr_stats.groupby(['dataset', 'input_total_peaks', 'peaks_in_promoter', 'motif','real_int', 'shuffled_int', 'p_val', 'enr_fold', 'adj_pval', 'pi_value', 'rank_pi_val']).agg({'gene_id': ','.join, 'family': lambda x: ','.join(list(set(x))), 'Any expressed gene': any}).reset_index().sort_values(by = 'rank_pi_val').drop('gene_id', axis = 1)
    
if not expressed_genes_file:

    enr_stats = enr_stats.groupby(['dataset', 'input_total_peaks', 'peaks_in_promoter', 'motif', 'real_int', 'shuffled_int', 'p_val', 'enr_fold', 'adj_pval', 'pi_value', 'rank_pi_val']).agg({'gene_id': ','.join, 'family': lambda x: ','.join(list(set(x)))}).reset_index().sort_values(by = 'rank_pi_val').drop('gene_id', axis = 1)

enr_stats = enr_stats.merge(mot_tf, how = 'right', left_on = 'motif', right_on = 'motif_id').drop('motif_id', axis = 1)

enr_stats = enr_stats.merge(info_df, how = 'left', on = 'gene_id')

if DE_table_file:
    enr_stats = enr_stats.merge(DE_table, how = 'left', on = 'gene_id')

if expressed_genes_file:
    enr_stats = enr_stats.groupby(['dataset', 'input_total_peaks', 'peaks_in_promoter', 'motif', 'real_int', 'shuffled_int', 'p_val', 'enr_fold', 'adj_pval', 'pi_value', 'rank_pi_val', 'family', 'Any expressed gene']).agg(lambda x: ','.join((x.dropna()))).reset_index().sort_values(by = 'rank_pi_val')
    
if not expressed_genes_file:
    enr_stats = enr_stats.groupby(['dataset', 'input_total_peaks', 'peaks_in_promoter', 'motif', 'real_int', 'shuffled_int', 'p_val', 'enr_fold', 'adj_pval', 'pi_value', 'rank_pi_val', 'family']).agg(lambda x: ','.join((x.dropna()))).reset_index().sort_values(by = 'rank_pi_val')


enr_stats = enr_stats.rename(columns = {'dataset': 'Dataset name', 'input_total_peaks': 'Total number of input peaks', 'peaks_in_promoter': 'Number of peaks within promoter definition', 'gene_id': 'Gene ID', 'TF_rank': 'TF rank', 'gene_name': 'Gene name', 'description': 'Gene description', 'family': 'TF family', 'ath_gene_id': 'Arabidopsis gene ID', 'ath_alias': 'Arabidopsis gene name', 'gene_name_figs': 'Gene name (+ alias)', 'motif': 'Motif IDs', 'enr_fold': 'Motif enrichement fold', 'p_val': 'Motif enrichment p-value', 'adj_pval': 'Motif enrichment q-value', 'pi_value': 'Motif enrichment pi-value', 'rank_pi_val': 'Motif rank (pi-value)', 'real_int': 'Real motif overlaps', 'shuffled_int': 'Shuffled motif overlaps'})

### Writing output file ###

if (output_type == 'csv'):
    enr_stats.to_csv(output_file, index = False)
else:
    with pd.ExcelWriter(output_file) as writer:
        enr_stats.to_excel(writer, index = False)
