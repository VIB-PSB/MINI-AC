import os
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to filter a network file for genes ' + \
                                            'only present in a specified gene set file',
                        conflict_handler='resolve')

    parser.add_argument('reduced_go_gene', nargs = 1, type = str,
                        help = '',
                        metavar = 'Reduced GO Gene File')
    
    parser.add_argument('GO_enrichemnt', nargs = 1, type = str,
                        help = '',
                        metavar = 'Non reduced GO enrichment file')

    parser.add_argument('output', nargs = 1, type = str,
                        help = '',
                        metavar = 'Reduced GO enrichment file')

    args = parser.parse_args()

    return args

args = parseArgs()

reduced_go_gene = args.reduced_go_gene[0]
go_enr = args.GO_enrichemnt[0]
output = args.output[0]

go_genes = set()

with open(reduced_go_gene, 'r') as f:
    for file in f:
        rec = file.strip().split("\t")
        go = rec[0]
        gene = rec[1]
        go_genes.add(tuple([go, gene]))
        
with open(output, 'w') as fout:
    fout.write("\t".join(["#set_id", "ftr_id", "p-val", "q-val", "enr_fold", "set_size", "ftr_size", "n_hits", "hits", "go_term"]))
    fout.write("\n")
    with open(go_enr, 'r') as fin:
        for line in fin:
            rec = line.strip().split("\t")
            go = rec[1]
            gene = rec[0]
            if (go, gene) in go_genes:
                fout.write(line)