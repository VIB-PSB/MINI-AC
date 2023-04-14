### Importing modules ###

import pandas as pd
import numpy as np
import math
import sys
import argparse

### Function to read input arguments ###

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to get TF-centric integrative MINI-AC output', 
                        conflict_handler = 'resolve')

    parser.add_argument('motif_enrichment_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'MINI-AC motif enrichment raw output file')

    parser.add_argument('network_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'MINI-AC network raw output file')

    parser.add_argument('go_enrichment_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'MINI-AC GO enrichment raw output file')

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

    args = parser.parse_args()

    return args

### Reading inputs ###

args = parseArgs()

motif_enrichment_file = args.motif_enrichment_file[0]
network_file = args.network_file[0]
go_enrichment_file = args.go_enrichment_file[0]
mot_tf_file = args.mot_tf_file[0]
tf_fam_file = args.tf_fam_file[0]
info_file = args.info_file[0]
pval = args.p_val
output_file = args.output_file[0]
DE_table_file = args.de_table_file
expressed_genes_file = args.expressed_genes_file

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
for col in enr_stats.select_dtypes(include = ['float']).columns:
    if col == 'rank_pi_val':
        continue
    enr_stats[col] = enr_stats[col].apply(lambda x: x if(math.isnan(x)) else round(x,3))

enr_stats = enr_stats.astype({'real_int':'int'})
enr_stats = enr_stats.astype({'shuffled_int':'int'})
enr_stats = enr_stats.astype({'rank_pi_val':'int'})

enr_stats['min_motif_rank'] = enr_stats.rank_pi_val
enr_stats['motif_min_rank'] = enr_stats.motif
enr_stats['min_qval'] = enr_stats.adj_pval
enr_stats = enr_stats.sort_values(by = 'rank_pi_val')
enr_stats = enr_stats.astype({'min_motif_rank':'int'})
enr_stats = enr_stats.astype({'min_qval':'float'})
enr_stats = enr_stats.astype({'real_int':'string', 'shuffled_int':'string',
                              'enr_fold':'string', 'adj_pval':'string',
                              'pi_value':'string', 'rank_pi_val':'string'})

enr_stats = enr_stats.groupby('gene_id').agg({'dataset': 'first', 'input_total_peaks': 'first', 'peaks_in_promoter': 'first', 'motif_min_rank': 'first', 'min_motif_rank': min, 'min_qval': min, 'motif': ','.join,'enr_fold': ','.join, 'adj_pval': ','.join,'rank_pi_val': ','.join}).sort_values(by = 'min_motif_rank').reset_index()
enr_stats = enr_stats.merge(tf_fam, how = 'left', on = 'gene_id').merge(info_df, how = 'left', on = 'gene_id')

if enr_stats.empty:
    empty_table = pd.DataFrame(["### This dataset did not yield any motif enrichment"])
    with pd.ExcelWriter(output_file) as writer:
        empty_table.to_excel(writer, index = False, header = False)
    sys.exit()

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
            
if GO_info:                
    go_df = pd.DataFrame().from_dict(GO_info, orient = 'index')
    go_df.index = pd.MultiIndex.from_tuples(go_df.index, names=['gene_id', 'go_term'])
    go_df = go_df.reset_index()
    go_df.columns = ['gene_id', 'go_term', 'go_id', "p_value",
                          "q_value_go", "enr_fold", 'total_genes_with_GO',
                          'total_target_genes', 'number_target_genes_with_GO', "target_genes_with_GO"]

    go_df = go_df.astype({'number_target_genes_with_GO': int})

    go_df = go_df[go_df.number_target_genes_with_GO > 1]

    go_df = go_df[['gene_id', 'go_term', 'q_value_go', 'enr_fold', 'total_genes_with_GO', 'total_target_genes', 'number_target_genes_with_GO', 'target_genes_with_GO']]


    go_df = go_df.astype({'q_value_go': 'float'})
    go_df.q_value_go = go_df.q_value_go.apply(lambda x: '{:0.3e}'.format(x))
    go_df = go_df.astype({'q_value_go': 'str'})

    go_df = go_df.groupby('gene_id').agg({'go_term' : ','.join, 'q_value_go' : ','.join}).reset_index()
    go_df_cols = list(go_df.columns[1:])

### Reading and processing network data ###

net_dict = {}
with open(network_file, 'r') as f:
    for line in f:
        if line.startswith("###"):
            net_dict = None
            break
        rec = line.strip().split("\t")
        net_dict.setdefault(rec[0], set()).add(rec[1])
        
if net_dict:
        
    if DE_table_file:
        net_dict_de_tgs_perc = {}
        for tf in net_dict:
            tgs = net_dict[tf].intersection(DE_genes)
            perc = len(tgs) / len(net_dict[tf])
            net_dict_de_tgs_perc.setdefault(tf, [len(net_dict[tf]), len(tgs), round(perc, 2)])

        df_tgs_info = pd.DataFrame().from_dict(net_dict_de_tgs_perc).T.reset_index()
        df_tgs_info.columns = ['gene_id', 'total TGs', 'DE TGs', '% DE TGs']
        df_tgs_info = df_tgs_info.astype({'DE TGs':'int'})

    elif not DE_table_file:
        net_dict_de_tgs_perc = {}
        for tf in net_dict:
            net_dict_de_tgs_perc.setdefault(tf, [len(net_dict[tf])])

        df_tgs_info = pd.DataFrame().from_dict(net_dict_de_tgs_perc).T.reset_index()
        df_tgs_info.columns = ['gene_id', 'total TGs']

    net_df_cols = list(df_tgs_info.columns[1:])

### Integrating everything ###

if DE_table_file:
    
    if expressed_genes_file:
        enr_stats['Expressed gene'] = enr_stats.gene_id.isin(exp_genes)
        info_df_columns = info_df_columns + ['Expressed gene']
        
    enr_stats = enr_stats.merge(DE_table, how = 'left', on = 'gene_id')
    
    net_cols = []
    
    if net_dict:
        enr_stats = enr_stats.merge(df_tgs_info, how = 'left', on = 'gene_id')
        
        net_cols = net_df_cols
    
    if GO_info:
        enr_stats = enr_stats.merge(go_df, how = 'left', on = 'gene_id')
        
        net_cols = net_cols + go_df_cols
    
    enr_stats = enr_stats[['dataset', 'input_total_peaks', 'peaks_in_promoter', 'gene_id',  'motif_min_rank', 'min_motif_rank', 'min_qval', 'family'] + info_df_columns + DE_table_columns + net_cols + ['motif', 'enr_fold', 'adj_pval', 'rank_pi_val']]

elif not DE_table_file:
    
    if expressed_genes_file:
        enr_stats['Expressed gene'] = enr_stats.gene_id.isin(exp_genes)
        info_df_columns = info_df_columns + ['Expressed gene']
        
    net_cols = []
    
    if net_dict:
        enr_stats = enr_stats.merge(df_tgs_info, how = 'left', on = 'gene_id')
        
        net_cols = net_df_cols
    
    if GO_info:
        enr_stats = enr_stats.merge(go_df, how = 'left', on = 'gene_id')
        
        net_cols = net_cols + go_df_cols
        
    enr_stats = enr_stats[['dataset', 'input_total_peaks', 'peaks_in_promoter', 'gene_id',  'motif_min_rank', 'min_motif_rank', 'min_qval', 'family'] + info_df_columns + net_cols + ['motif', 'enr_fold', 'adj_pval', 'rank_pi_val']]
    
enr_stats.insert(4, 'TF_rank', enr_stats.min_motif_rank.rank(axis = 0, method = 'dense'))
enr_stats = enr_stats.astype({'TF_rank':'int'})

enr_stats = enr_stats.rename(columns = {'dataset': 'Dataset name', 'input_total_peaks': 'Total number of input peaks', 'peaks_in_promoter': 'Number of peaks within promoter definition', 'gene_id': 'Gene ID', 'TF_rank': 'TF rank', 'motif_min_rank': 'Motif with min rank', 'min_motif_rank': 'Rank of motif', 'min_qval': 'q-value of motif with min rank', 'gene_name': 'Gene name', 'description': 'Gene description', 'family': 'TF family', 'ath_gene_id': 'Arabidopsis gene ID', 'ath_alias': 'Arabidopsis gene name', 'gene_name_figs': 'Gene name (+ alias)', 'go_term': 'Enriched GO terms', 'q_value_go': 'GO enrichment q-value', 'motif': 'TF motifs IDs', 'enr_fold': 'Motif enrichment folds', 'adj_pval': 'Motif enrichment q-values', 'rank_pi_val': 'Motif ranks (pi-value)'})

### Writing output file ###

with pd.ExcelWriter(output_file) as writer:
    enr_stats.to_excel(writer, index = False)