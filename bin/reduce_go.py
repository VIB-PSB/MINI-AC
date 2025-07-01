#!/usr/bin/python
import go_manipulations
from os import path
from sys import argv

gene_go_file = argv[1]
ontology_file = argv[2]

go_tree = go_manipulations.GOtree(ontology_file)

gene_go = go_manipulations.GOgenes(gene_go_file, go_tree)
go_tree.reduce(gene_go)
gene_go.write()
