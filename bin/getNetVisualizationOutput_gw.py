### Importing modules ###

import pandas as pd
import numpy as np
import sys
import argparse

### Function to read input arguments ###

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to get network visualization files of MINI-AC output ready to use with Cytoscape', 
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
                        default = 0.1, help = '',
                        metavar = 'p-value for network inference')

    parser.add_argument('output_file_func_net', nargs = 1, type = str,
                        help = '',
                        metavar = 'Name of the output file')

    parser.add_argument('output_file_nodes_atts', nargs = 1, type = str,
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
output_file_func_net = args.output_file_func_net[0]
output_file_nodes_atts = args.output_file_nodes_atts[0]
DE_table_file = args.de_table_file
expressed_genes_file = args.expressed_genes_file

info_df = pd.read_csv(info_file).drop_duplicates(subset = "gene_id")
mot_tf = pd.read_csv(mot_tf_file)
tf_fam = pd.read_csv(tf_fam_file)

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
            DE_table.columns.values[col_num] = DE_table.columns.values[col_num] + "_DE"

    lst = list(DE_table)
    DE_table[lst] = DE_table[lst].astype(str)

    DE_genes = set(DE_table.gene_id)
    DE_table_columns = list(DE_table.columns)
    DE_table_columns = DE_table_columns[1:]

### Reading and processing network data ###

net_list = []
with open(network_file, 'r') as f:
    for line in f:
        if line.startswith("###"):
            net_list = None
            break
        if line.startswith("#TF"):
            continue
        rec = line.strip().split("\t")
        net_list.append([rec[0], rec[1]])
        
if net_list:
    net_df = pd.DataFrame(net_list)
    net_df.columns = ['TF', 'TG']
    
if not net_list:
    empty_table = pd.DataFrame(["### This dataset did not yield any motif enrichment or network"])
    empty_table.to_csv(output_file_func_net, sep = "\t", index = None, header = None)
    empty_table.to_csv(output_file_nodes_atts, sep = "\t", index = None, header = None)
    sys.exit()

### Reading and processing GO enrichment data ###

GO_net = []

with open(go_enrichment_file, 'r') as f:
    for line in f:
        rec = line.strip().split("\t")
        if line.startswith("###"):
            GO_net = None
            break
        if line.startswith("#set_id"):
            continue
        tf, go_term, q_val, enr_fold = rec[0], rec[-1], rec[3], rec[4]
        tgs = rec[-2].split(",")
        for tg in tgs:
            GO_net.append([tf, tg, go_term, q_val, enr_fold])

if GO_net:
    go_df = pd.DataFrame(GO_net)
    go_df.columns = ['TF', 'TG', 'GO_term', 'q-value', 'enrichment_fold']
    go_df = go_df.astype({'q-value': 'float'})
    go_df.loc[:, 'q-value'] = go_df.loc[:, 'q-value'].apply(lambda x: '{:0.3e}'.format(x))
    func_net = net_df.merge(go_df, how = 'left', on = ['TF', 'TG'])
if not GO_net:
    func_net = net_df

### Reading and processing motif enrichment data ###

enr_stats = pd.read_csv(motif_enrichment_file, sep = "\t")
enr_stats.insert(7, 'pi_value', enr_stats.enr_fold * (-np.log10(enr_stats.adj_pval)))
enr_stats = enr_stats.sort_values(by = 'pi_value', ascending = False)
enr_stats["rank_pi_val"] = enr_stats.groupby("dataset")["pi_value"].rank("first", ascending = False)
enr_stats = enr_stats.merge(mot_tf, how = 'right', left_on = 'motif', right_on = 'motif_id')
enr_stats = enr_stats[enr_stats.adj_pval <= pval]
enr_stats = enr_stats.groupby('gene_id').agg({'motif': 'first', 'rank_pi_val': min, 'adj_pval': min, 'enr_fold': max, 'pi_value': max}).sort_values(by = 'rank_pi_val').reset_index()
enr_stats.columns = ['gene_id', 'motif_min_rank', 'min_motif_rank', 'min_q_val', 'min_enr_fold', 'min_pi_value']
enr_stats = enr_stats.astype({'min_motif_rank':'int'})

if enr_stats.empty:
    empty_table = pd.DataFrame(["### This dataset did not yield any motif enrichment or network"])
    empty_table.to_csv(output_file_func_net, sep = "\t", index = None, header = None)
    empty_table.to_csv(output_file_nodes_atts, sep = "\t", index = None, header = None)
    sys.exit()

### Integrating everything ###

tf_df = pd.DataFrame(func_net.TF)
tf_df['type'] = 'TF'
tf_df.columns = ['gene_id', 'type']
tg_df = pd.DataFrame(func_net.TG)
tg_df['type'] = 'TG'
tg_df.columns = ['gene_id', 'type']
nodes_table = pd.concat([tf_df, tg_df]).drop_duplicates(subset = "gene_id")
nodes_table = nodes_table.merge(enr_stats, how = 'left', on = 'gene_id').merge(info_df, how = 'left', on = 'gene_id').drop_duplicates(subset = "gene_id")

if expressed_genes_file:
    nodes_table['Expressed gene'] = nodes_table.gene_id.isin(exp_genes)

if DE_table_file:
    nodes_table = nodes_table.merge(DE_table, how = 'left', on = 'gene_id').drop_duplicates(subset = "gene_id")

### Writing output file ###

func_net.to_csv(output_file_func_net, sep = "\t", index = None)
nodes_table.to_csv(output_file_nodes_atts, sep = "\t", index = None)