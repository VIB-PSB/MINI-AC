# %%
import argparse

def parseArgs():

    parser = argparse.ArgumentParser(prog = 'Script to get a BED file with iCREs ' + \
                                            'coordinates given a list of genes',
                        conflict_handler='resolve')

    parser.add_argument('annotated_icres', nargs = 1, type = str,
                        help = '',
                        metavar = 'BED file with 4th column being ' +\
                                    'an annotated gene ID')
    
    parser.add_argument('gene_list', nargs = 1, type = str,
                        help = '',
                        metavar = 'One column file containing gene IDs '+ \
                                'of interest')

    parser.add_argument('bed_of_genes_icres', nargs = 1, type = str,
                        help = '',
                        metavar = 'Output BED file with coordinates '+\
                            'of iCREs associated with genes of interest')

    args = parser.parse_args()

    return args

args = parseArgs()

annot_icres = args.annotated_icres[0]
genes_oi_file = args.gene_list[0]
output_file = args.bed_of_genes_icres[0]

# %%
genes_oi = set()

with open(genes_oi_file, "r") as fin:
    for line in fin:
        rec = line.strip().split("\t")
        gene_id = rec[0]
        genes_oi.add(gene_id)

with open(output_file, "w") as fout:
    with open(annot_icres, "r") as fin:
        for line in fin:
            rec = line.strip().split("\t")
            gene_id = rec[3]
            if gene_id in genes_oi:
                fout.write("\t".join(rec[0:3]))
                fout.write("\n")