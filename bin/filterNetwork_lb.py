import os
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to filter a network file for genes ' + \
                                            'only present in a specified gene set file',
                        conflict_handler='resolve')

    parser.add_argument('input_network', nargs = 1, type = str,
                        help = '',
                        metavar = 'Input network file')
    
    parser.add_argument('gene_set', nargs = 1, type = str,
                        help = '',
                        metavar = 'One column file containing gene IDs '+ \
                                  'of the set for filtering')

    parser.add_argument('output_network', nargs = 1, type = str,
                        help = '',
                        metavar = 'Output network file filtered for genes in set file')

    args = parser.parse_args()

    return args

args = parseArgs()

input_file = args.input_network[0]
filtering_set = args.gene_set[0]
output_file = args.output_network[0]

filt_set = []
with open(filtering_set, 'r') as f:
    for line in f:
        gene = line.strip().split(",")[0]
        filt_set.append(gene)
        
filt_set = set(filt_set)

def wcl(file):
    num_lines = sum(1 for line in open(file))
    return num_lines

if wcl(input_file) == 1:
    with open(output_file, 'w') as fout:
        fout.write("### It was not possible to predict a network for this dataset")

with open(output_file, 'w') as fout:
    with open(input_file, 'r') as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            rec = line.strip().split("\t")
            TF = rec[0]
            TG = rec[1]
            if TF in filt_set and TG in filt_set:
                fout.write("\t".join([TF, TG]))
                fout.write("\n")
            

