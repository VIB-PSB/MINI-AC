#!/usr/bin/python
import go_manipulations
from os import path
from sys import argv

gene_go_file = argv[1]

go_tree = go_manipulations.GOtree(path.join(path.dirname(path.dirname(argv[0])), "ontologies", "go.obo"))

gene_go = go_manipulations.GOgenes(gene_go_file, go_tree)
go_tree.reduce(gene_go)
gene_go.write()
