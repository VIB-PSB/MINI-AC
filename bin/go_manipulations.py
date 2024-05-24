#!/usr/bin/python
from collections import deque
from sys import argv, stderr

parental_links = set(
    ["negatively_regulates", "positively_regulates", "regulates", "part_of", "is_a"])

exp_evidence = ["EXP", "IMP", "IDA", "IPI", "IGI", "IEP"]
aut_evidence = ["TAS", "NAS", "IC"]
com_evidence = ["ISO", "ISS", "ISA", "IBA", "IBD",
    "ISM", "IKR", "IRD", "IGC", "RCA", "IEA"]
non_evidence = ["ND", "NA", "extended"]
all_evidence = exp_evidence + aut_evidence + com_evidence + non_evidence
evidence_priority = dict()
counter = 0
for evidence_code in all_evidence:
	evidence_priority[evidence_code] = counter
	counter += 1


class GOterm:
	def __init__(self, go_id):
		self.go_id = go_id
		self.go_name = ""
		self.go_namespace = ""
		self.go_subsets = set()
		self.go_parents = set()
		self.go_children = set()


class GOtree:
	# Constructor
	def __init__(self, obo_file):
		self.go_terms = dict()
		self.alt_ids = dict()
		with open(obo_file, 'r') as reader:
			line = ""
			in_term_field = False
			term_id = ""
			for line in reader:
				split = line.strip().split(" ")
				in_term_field = split[0] == "[Term]" or (
				    in_term_field and not split[0] == "")
				if not in_term_field: continue

				if split[0] == "id:":
					term_id = split[1]
					if not term_id in self.go_terms:
						self.go_terms[term_id] = GOterm(term_id)
				elif split[0] == "name:":
					self.go_terms[term_id].go_name = " ".join(split[1:len(split)])
				elif split[0] == "is_a:":
					self.go_terms[term_id].go_parents.add(split[1])
					if not split[1] in self.go_terms:
						self.go_terms[split[1]] = GOterm(split[1])
					self.go_terms[split[1]].go_children.add(term_id)
				elif split[0] == "namespace:":
					self.go_terms[term_id].go_namespace = split[1]
				elif split[0] == "relationship:" and split[1] in parental_links:
					self.go_terms[term_id].go_parents.add(split[2])
					if not split[2] in self.go_terms:
						self.go_terms[split[2]] = GOterm(split[2])
					self.go_terms[split[2]].go_children.add(term_id)
				elif split[0] == "alt_id:":
					self.alt_ids[split[1]] = term_id

	# Getters
	def size(self):
		return len(self.go_terms)

	def get_namespace(self, go_id):
		go_id = self.get_consensus_id(go_id)
		return self.go_terms[go_id].go_namespace

	def get_parents(self, go_id):
		go_id = self.get_consensus_id(go_id)
		return self.go_terms[go_id].go_parents

	def get_children(self, go_id):
		go_id = self.get_consensus_id(go_id)
		return self.go_terms[go_id].go_children

	def get_ancestors(self, go_id):
		go_id = self.get_consensus_id(go_id)
		ancestors = set(self.get_parents(go_id))
		todo = deque(ancestors)
		while todo:
			next_parent = todo.pop()
			for ancestor in self.get_parents(next_parent):
				if not ancestor in ancestors:
					ancestors.add(ancestor)
					todo.append(ancestor)
		return ancestors

	def get_descendants(self, go_id):
		go_id = self.get_consensus_id(go_id)
		descendants = set(self.get_children(go_id))
		todo = deque(descendants)
		while todo:
			next_child = todo.pop()
			for descendant in self.get_children(next_child):
				if not descendant in descendants:
					descendants.add(descendant)
					todo.append(descendant)
		return descendants

	# Extend/Reduce
	def extend(self, go_genes):
		for gene_id in go_genes.gene_go:
			known_gos = list(go_genes.gene_go[gene_id])
			for known_go in known_gos:
				ancestors = self.get_ancestors(known_go)
				for ancestor in ancestors:
					go_genes.add(gene_id, ancestor, "extended")

	#def extend_downstream(self, go_genes):
	#	for gene_id in go_genes.gene_go:
	#		known_gos = list(go_genes.gene_go[gene_id])
	#		for known_go in known_gos:
	#			descendants = self.get_descendants(known_go)
	#			for descendant in descendants:
	#				go_genes.add(gene_id, descendant, "extended")

	def reduce(self, go_genes):
		for gene_id in go_genes.gene_go:
			known_gos = list(go_genes.gene_go[gene_id])
			for known_go in known_gos:
				ancestors = self.get_ancestors(known_go)
				for ancestor in ancestors:
					go_genes.remove(gene_id, ancestor)

	def filter_category(self, go_genes, select):
		if not select in ["molecular_function", "biological_process", "cellular_component", "external"]:
			stderr.write("[ERROR] Category filter not recognized: " + select + "\n")
			return
		all_genes = go_genes.gene_go.keys()
		for gene_id in all_genes:
			known_gos = list(go_genes.gene_go[gene_id])
			for known_go in known_gos:
				if not self.go_terms[known_go].go_namespace == select:
					go_genes.remove(gene_id, known_go)

	def filter_evidence(self, go_genes, select):
		if not select in ["exp", "com", "cur", "non"]: return
		filter_list = all_evidence
		if select == "exp": filter_list = exp_evidence
		if select == "com": filter_list = com_evidence
		if select == "cur": filter_list = aut_evidence
		if select == "non": filter_list = non_evidence
		all_genes = go_genes.gene_go.keys()
		for gene_id in all_genes:
			known_gos = list(go_genes.gene_go[gene_id])
			for known_go in known_gos:
				if not go_genes.go_evidence[gene_id + "-" + known_go] in filter_list:
					go_genes.remove(gene_id, known_go)

	def filter_content(self, go_genes):
		all_genes = go_genes.gene_go.keys()
		for gene_id in all_genes:
			known_gos = list(go_genes.gene_go[gene_id])
			for known_go in known_gos:
				if not known_go in self.go_terms:
					go_genes.remove(gene_id, known_go)

	def get_consensus_id(self, go_id):
		if go_id in self.alt_ids:
			return self.alt_ids[go_id]
		elif not go_id in self.go_terms:
			stderr.write("[WARN] No such GO term: "+go_id+"\n")
			return go_id
		else:
			return go_id

	# Add descriptions
	def add_descriptions_to_file(self, in_file, out_file):
		with open(out_file, 'r') as writer:
			with open(in_file, 'r') as reader:
				for line in reader:
					if line[0] == "#":
						print(line.strip())
						continue
					description = "NA"
					split = line.strip().split("\t")
					for col in split:
						if col in self.go_terms:
							description = self.go_terms[col].go_name
					writer.write(line.strip()+"\t"+description+"\n")
	def add_descriptions(self, in_file):
		with open(in_file, 'r') as reader:
			for line in reader:
				if line[0] == "#":
					print(line.strip())
					continue
				description = "NA"
				split = line.strip().split("\t")
				for col in split:
					if col in self.go_terms:
						description = self.go_terms[col].go_name
				print(line.strip() + "\t" + description)

class GOgenes:
	# Constructor
	def __init__(self, go_file, go_tree):
		self.gene_go = dict() 		# { gene_id => set( go_id ) }
		self.go_gene = dict() 		# { go_id => set( gene_id ) }
		self.go_evidence = dict() 	# { "gene_id"+"go_id" => evidence }
		with open(go_file, 'r') as reader:
			for line in reader:
				if line[0] == "#": continue
				# Get go_id, gene_id and evidence_code
				split = line.strip().split("\t")
				gene_id = split[1]
				go_id = go_tree.get_consensus_id(split[0])
				if gene_id[0:3] == "GO:": (go_id, gene_id) = (gene_id, go_id)
				evidence_code = "NA"
				if len(split) > 2 and split[2] in all_evidence: evidence_code = split[2]
				# Add the data
				self.add(gene_id, go_id, evidence_code)
	# Adding/Removing annotations
	def add(self, gene_id, go_id, evidence_code):
		# Check if the go_id and gene_id were already known
		if not gene_id in self.gene_go: self.gene_go[gene_id] = set()
		if not go_id   in self.go_gene: self.go_gene[go_id]   = set()
		# Store gene-go and go-gene info
		self.gene_go[gene_id].add(go_id)
		self.go_gene[go_id].add(gene_id)
		# Store the evidence code if it's priority is higher than already known evidence codes
		gene_go_string = gene_id+"-"+go_id
		if not gene_go_string in self.go_evidence or evidence_priority[self.go_evidence[gene_go_string]] > evidence_priority[evidence_code]:
			self.go_evidence[gene_go_string] = evidence_code
	def remove(self, gene_id, go_id):
		# Remove the gene-go annotation
		if gene_id in self.gene_go:
			if go_id   in self.gene_go[gene_id]: self.gene_go[gene_id].remove(go_id)
			if not self.gene_go[gene_id]: del self.gene_go[gene_id]
		if go_id   in self.go_gene:
			if gene_id in self.go_gene[go_id]:   self.go_gene[go_id].remove(gene_id)
			if not self.go_gene[go_id]: del self.go_gene[go_id]
		# Remove the evidence code
		if gene_id+"-"+go_id in self.go_evidence: del self.go_evidence[gene_id+"-"+go_id]
	# Getters
	def get_n_gen(self):
		return len(self.gene_go)
	def get_n_go(self):
		return len(self.go_gene)
	# Write
	def write_to_file(self, out_file):
		with open(out_file, 'w') as writer:
			for gene_id in sorted(self.gene_go):
				for go_id in sorted(self.gene_go[gene_id]):
					writer.write(go_id+"\t"+gene_id+"\t"+self.go_evidence[gene_id+"-"+go_id]+"\n")
	def write(self):
		for gene_id in sorted(self.gene_go):
			for go_id in sorted(self.gene_go[gene_id]):
				if self.go_evidence[gene_id+"-"+go_id] != "NA":
					print(go_id+"\t"+gene_id+"\t"+self.go_evidence[gene_id+"-"+go_id])
				else:
					print(go_id+"\t"+gene_id)
