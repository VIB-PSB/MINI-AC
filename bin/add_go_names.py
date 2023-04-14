#!/usr/bin/python
from sys import argv
from os import path
import go_manipulations

gene_go_file = argv[1]

go_tree = go_manipulations.GOtree(path.join(path.dirname(path.dirname(argv[0])), "ontologies", "go.obo"))

go_tree.add_descriptions(gene_go_file)
