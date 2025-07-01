import pandas as pd
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to process MINIAC enrichment ' + \
                                            'results and obtain ACR based network',
                        conflict_handler='resolve')

    parser.add_argument('peak_ann', nargs = 1, type = str,
                        help = '',
                        metavar = 'motif-target gene file')

    parser.add_argument('stats_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'File with the enrichment statistics of the ACR file')

    parser.add_argument('mot_tf_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'Motif-TF file')

    parser.add_argument('-sg', '--second_gene_dist', nargs = '?', type = int,
                        help = '', default = 0,
                        metavar = 'Flag for annotating the 2nd closest '+\
                            'gene, and what distance difference with '+\
                            'the 1st annotated gene is the threshold '+\
                            'to annotate such gene')

    parser.add_argument('p_val', nargs = '?', type = float,
                        default = 0.1, help = '',
                        metavar = 'p-value for network inference')

    parser.add_argument('output', nargs = 1, type = str,
                        help = '',
                        metavar = 'Name of the output file')

    args = parser.parse_args()

    return args

args = parseArgs()

peak_ann_file = args.peak_ann[0]
stats_input = args.stats_file[0]
motif_tf_file = args.mot_tf_file[0]
chosen_dist = args.second_gene_dist
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

    peak_ann_d = {}
    with open(peak_ann_file, 'r') as f:
        for line in f:
            rec = line.strip().split("\t")
            peak_ann_d.setdefault(tuple(rec[0:4]), list()).append([rec[-1], rec[-2]])

    if chosen_dist == 0:

        with open(out_file, 'w') as fout:

            fout.write("\t".join(["#TF", "TG"]))
            fout.write("\n")
            
            for coords in peak_ann_d:
                for i in peak_ann_d[coords]:
                    mot, tg = coords[3], i[1]
                    if tg == ".":
                        continue
                    if mot in enr_mots:
                        for tf in mot_tf_dic[mot]:
                            fout.write("\t".join([tf, tg]))
                            fout.write("\n")

    if chosen_dist != 0:

        with open(out_file, 'w') as fout:

            fout.write("\t".join(["#TF", "TG"]))
            fout.write("\n")
            
            for coords in peak_ann_d:
                if len(peak_ann_d[coords]) >= 2:
                    for i in peak_ann_d[coords]:
                        if i == peak_ann_d[coords][0]:
                            mot, tg = coords[3], peak_ann_d[coords][0][1]
                            if tg == ".":
                                continue
                            if mot in enr_mots:
                                for tf in mot_tf_dic[mot]:
                                    fout.write("\t".join([tf, tg]))
                                    fout.write("\n")
                        else:
                            other_dist = abs(int(i[0]))
                            if other_dist <= chosen_dist:
                                mot, tg = coords[3], i[1]
                                if tg == ".":
                                    continue
                                if mot in enr_mots:
                                    for tf in mot_tf_dic[mot]:
                                        fout.write("\t".join([tf, tg]))
                                        fout.write("\n")

                if len(peak_ann_d[coords]) == 1:
                    mot, tg = coords[3], peak_ann_d[coords][0][1]
                    if tg == ".":
                        continue
                    if mot in enr_mots:
                        for tf in mot_tf_dic[mot]:
                            fout.write("\t".join([tf, tg]))
                            fout.write("\n")


