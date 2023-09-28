#!/usr/bin/python
from sys import argv
from os import path
import go_manipulations

gene_go_file = argv[1]
ontology_file = argv[2]

go_tree = go_manipulations.GOtree(ontology_file)

go_tree.add_descriptions(gene_go_file)
