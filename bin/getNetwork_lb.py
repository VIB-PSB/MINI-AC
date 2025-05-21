import pandas as pd
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to process MINIAC enrichment ' + \
                                            'results and obtain ACR based network',
                        conflict_handler='resolve')

    parser.add_argument('mot_tg', nargs = 1, type = str,
                        help = '',
                        metavar = 'motif-target gene file')

    parser.add_argument('stats_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'File with the enrichment statistics of the ACR file')

    parser.add_argument('mot_tf_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'Motif-TF file')

    parser.add_argument('p_val', nargs = '?', type = float,
                        default = 0.01, help = '',
                        metavar = 'p-value for network inference')

    parser.add_argument('output', nargs = 1, type = str,
                        help = '',
                        metavar = 'Name of the output file')

    args = parser.parse_args()

    return args

args = parseArgs()

mot_tg_file = args.mot_tg[0]
stats_input = args.stats_file[0]
motif_tf_file = args.mot_tf_file[0]
p_value = args.p_val
out_file = args.output[0]

mot_tf_dic = {}

with open(motif_tf_file, 'r') as f:
    for line in f:
        if line.startswith('motif_id'):
            continue
        rec = line.strip().split(",")
        motif, tf = rec[0], rec[1]
        mot_tf_dic.setdefault(motif, set()).add(tf)


stats_df = pd.read_csv(stats_input, sep = "\t")

enr_mots = list(stats_df[stats_df.adj_pval <= p_value].motif)


if not enr_mots:
    with open(out_file, 'w') as fout:
        fout.write("### It was not possible to predict a network for this dataset")
        
else:

    with open(out_file, 'w') as fout:
        fout.write("\t".join(["#TF", "TG"]))
        fout.write("\n")
        with open(mot_tg_file, 'r') as fin:
            for line in fin:
                rec = line.strip().split("\t")
                mot, tg = rec[0], rec[1]
                if mot in enr_mots:
                    for tf in mot_tf_dic[mot]:
                        fout.write("\t".join([tf, tg]))
                        fout.write("\n")
