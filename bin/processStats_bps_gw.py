import warnings
warnings.filterwarnings('ignore')
import numpy as np
from statsmodels.stats.multitest import multipletests
import pandas as pd
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to process enrichment raw numbers ' + \
                                            'into statistics',
                        conflict_handler='resolve')

    parser.add_argument('stats_input_file', nargs = 1, type = str,
                        help = '',
                        metavar = 'File with raw intersection numbers')

    parser.add_argument('num_peaks', nargs = 1, type = int,
                        help = '',
                        metavar = 'Number of peaks in input ACR real file')

    parser.add_argument('cns_sets_list', nargs = 1, type = str,
                        help = '',
                        metavar = 'File with list of CNS sets')

    parser.add_argument('output', nargs = 1, type = str,
                    help = '',
                    metavar = 'output_file')

    args = parser.parse_args()

    return args

args = parseArgs()

raw_file = args.stats_input_file[0]
out_file = args.output[0]
cns_list_file = args.cns_sets_list[0]
total_peaks = args.num_peaks[0]

file_name = "_".join(raw_file.split("/")[-1].split("_")[0:-5])

cns_list = set()
with open(cns_list_file, 'r') as f:
    for line in f:
        if line.startswith("motif_id"):
            continue
        rec = line.strip().split(",")
        cns_set = rec[0]
        cns_list.add(cns_set)

shuff_dict = {}
real_dict = {}
with open(raw_file, 'r') as fin:
    for line in fin:
        rec = line.strip().split("\t")
        cns_set = rec[3]
        gs_shuff_event = rec[5].split(";")
        int_bps = rec[4].split(";")
        for i in range(len(gs_shuff_event)):
            if gs_shuff_event[i] == "real_ints":
                real_dict.setdefault(cns_set, []).append(int(int_bps[i]))
            else:
                shuff_dict.setdefault(cns_set, {}).setdefault(gs_shuff_event[i], []).append(int(int_bps[i]))

for i in real_dict:
    real_dict.update({i: sum(real_dict[i])})

for cns_set in shuff_dict:
    summed_bps = []
    for shuff_event in shuff_dict[cns_set]:
        summed_bps.append(sum(shuff_dict[cns_set][shuff_event]))
    shuff_dict.update({cns_set: summed_bps})

stats = {}
for cns_set in cns_list:
    try:
        shuff_overlap = np.array(shuff_dict[cns_set])
    except KeyError:
        shuff_overlap = np.zeros(1000, dtype = int)
    try:
        real_overlap = real_dict[cns_set]
    except KeyError:
        real_overlap = 0
    if len(shuff_overlap) < 1000:
        zero_pad = np.zeros(1000 - len(shuff_overlap), dtype = int)
        shuff_overlap = np.concatenate([shuff_overlap, zero_pad])
    else:
        pass
    p_val_1000 = len(shuff_overlap[shuff_overlap >= real_overlap])
    if p_val_1000 == 0:
        p_val_1000 = 0.9
    median = np.median(shuff_overlap)
    if real_overlap == 0 or median == 0:
        enrichment_fold = 0
    else:
        enrichment_fold = real_overlap / median
    
    stats[(file_name, cns_set)] = [int(total_peaks), real_overlap, median, (p_val_1000/1000), enrichment_fold]

data_df = pd.DataFrame.from_dict(stats).T
data_df = data_df.reset_index()
data_df.columns = ['dataset', 'motif', 'input_total_peaks', 'real_int', 'shuffled_int', 'p_val', 'enr_fold']
data_df = data_df[['dataset', 'input_total_peaks', 'motif', 'real_int', 'shuffled_int', 'p_val', 'enr_fold']]
FDR = multipletests(data_df['p_val'], method = 'fdr_bh', alpha = 0.05)
data_df.insert(6, 'adj_pval', FDR[1])

data_df.to_csv(out_file, sep = "\t", index = None, na_rep = "nan")

