"""
Parsimony Tree Improvement 2017
By Sachandhan Ganesh

General Usage:
	parsimony.py -a <alignment> -n <newick_tree> -r <number_of_random_rearrangements> -p <probability>
"""

from Bio import AlignIO, Phylo # Format reader for clustal alignments and newick trees
from ete3 import Tree # Newick tree reader

class ParsimonyTree(object):
	def __init__(self, msa_filepath, tree_filepath):
		# Read files and unpack into iterable objects
		self._msa = AlignIO.read(msa_filepath, "clustal") # multiple sequene alignment object
		self._tree = Phylo.read(tree_filepath, "newick") # tree object of alignment

	def get_parsimony_score(self):
		total_score = 0

		# for each nucleotide position
		for i in range(self._msa.get_alignment_length()):
			score = 0 # score for tree given nucleotide states

			# Assign states for leaves
			clade_states = self.get_terminal_states(i)

			# postorder traversal of internal nodes to resolve commonalities
			for clade in self._tree.get_nonterminals(order="postorder"):
				children = clade.clades

				# get left and right children of internal node
				alpha = clade_states[children[0]]
				beta = clade_states[children[1]]

				# find common nucleotide states between children
				common_states = alpha.intersection(beta)

				# if there are common states
				if common_states:
					# use common state
					clade_states[clade] = common_states
				# otherwise
				else:
					# take all states and increment parsimony score
					clade_states[clade] = alpha.union(beta)
					score += 1

			total_score += score

		return total_score

	def find_alignment_index_from_clade(self, clade):
		clade = str(clade)

		for i in range(len(self._msa)):
			if self._msa[i].id == clade:
				return i

	def get_terminal_states(self, nt_pos):
		terminals = self._tree.get_terminals()
		nts = []

		for clade in terminals:
			i = self.find_alignment_index_from_clade(clade)
			nts.append(self._msa[i, nt_pos])

		states = {}
		for clade, nt in zip(terminals, nts):
			states[clade] = set(nt)

		return states


	def visialize_tree(self):
		Phylo.draw(self._tree)


def main():
	pt = ParsimonyTree("./test_data/test_msa.txt", "./test_data/test_tree.txt")
	score = pt.get_parsimony_score()
	print(score)

	pt.visialize_tree()

if __name__ == "__main__":
	main()
