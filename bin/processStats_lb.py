import numpy as np
from statsmodels.stats.multitest import multipletests
import pandas as pd
import argparse
import warnings
warnings.filterwarnings('ignore')

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to process MINIAC raw numbers ' + \
                                            'into statistics',
                        conflict_handler='resolve')

    parser.add_argument('stats_input_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'File with raw intersection numbers')

    parser.add_argument('num_peaks', nargs = 1, type = int,
                        help = '',
                        metavar = 'Number of peaks in input ACR real file')

    parser.add_argument('peaks_in_prom', nargs = 1, type = int,
                        help = '',
                        metavar = 'Number of peaks in input ' +\
                            'ACR inside promoter definition')

    parser.add_argument('mot_tf', nargs = 1, type = str,
                        help = '',
                        metavar = 'Motif_TF_file')

    parser.add_argument('output', nargs = 1, type = str,
                    help = '',
                    metavar = 'output_file')

    args = parser.parse_args()

    return args

def correct_pval(df):
    FDR = multipletests(df['p_val'], method = 'fdr_bh', alpha = 0.05)
    df.insert(6, 'adj_pval', FDR[1])
    return df

args = parseArgs()

raw_file = args.stats_input_file[0]
out_file = args.output[0]
mot_tf_file = args.mot_tf[0]
total_peaks = args.num_peaks[0]
peaksInProm = args.peaks_in_prom[0]

file_name = "_".join(raw_file.split("/")[-1].split("_")[0:-5])

data_real = {}
data_shuff = {}

with open(raw_file, 'r') as f:
    for line in f:
        rec = line.strip().split(" ")
        motif = rec[1]
        ints = int(rec[0])
        it_num = rec[2]
        if it_num == "real_ints":
            data_real[motif] = ints
        else:
            data_shuff.setdefault(motif, list()).append(ints)

mots = pd.read_csv(mot_tf_file).motif_id

stats = {}

for motif in mots:
    try:
        shuff_motmaps = np.array(data_shuff[motif])
    except KeyError:
        shuff_motmaps = np.zeros(1000, dtype = int)
    try:
        real_motmaps = data_real[motif]
    except KeyError:
        real_motmaps = 0
    if len(shuff_motmaps) < 1000:
        zero_pad = np.zeros(1000 - len(shuff_motmaps), dtype = int)
        shuff_motmaps = np.concatenate([shuff_motmaps, zero_pad])
    else:
        pass
    p_val_1000 = len(shuff_motmaps[shuff_motmaps >= real_motmaps])
    if p_val_1000 == 0:
        p_val_1000 = 0.9
    median = np.median(shuff_motmaps)
    if real_motmaps == 0 or median == 0:
        enrichment_fold = 0
    else:
        enrichment_fold = real_motmaps / median
    
    stats[(file_name, motif)] = [int(total_peaks), int(peaksInProm), real_motmaps, median, (p_val_1000/1000), enrichment_fold]


data_df = pd.DataFrame.from_dict(stats).T
data_df = data_df.reset_index()
data_df.columns = ['dataset', 'motif', 'input_total_peaks', 'peaks_in_promoter', 'real_int', 'shuffled_int', 'p_val', 'enr_fold']
data_df = data_df[['dataset', 'input_total_peaks', 'peaks_in_promoter', 'motif', 'real_int', 'shuffled_int', 'p_val', 'enr_fold']]
FDR = multipletests(data_df['p_val'], method = 'fdr_bh', alpha = 0.05)
data_df.insert(7, 'adj_pval', FDR[1])

data_df.to_csv(out_file, sep = "\t", index = None, na_rep = "nan")

