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
		print(self._msa)

		for i in range(self._msa.get_alignment_length()):
			score = 0

			for clade, nt in self.pair_msa_tree_leaves(i):
				print(clade, nt)

			break

	def find_alignment_index_from_clade(self, clade):
		clade = str(clade)

		for i in range(len(self._msa)):
			if self._msa[i].id == clade:
				return i

	def pair_msa_tree_leaves(self, nt_pos):
		terminals = self._tree.get_terminals()
		nts = []

		for clade in terminals:
			i = self.find_alignment_index_from_clade(clade)
			nts.append(self._msa[i, nt_pos])

		return list(zip(terminals, nts))


def main():
	pt = ParsimonyTree("./test_data/test_msa.txt", "./test_data/test_tree.txt")
	pt.get_parsimony_score()

if __name__ == "__main__":
	main()
