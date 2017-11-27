"""
Parsimony Tree Improvement 2017
By Sachandhan Ganesh

General Usage:
	parsimony.py -a <alignment> -n <newick_tree> -r <number_of_random_rearrangements> -p <probability>
"""

from Bio import AlignIO, Phylo # Format reader for clustal alignments and newick trees
from ete3 import Tree # Newick tree reader

class ParsimonyTree(object):
	self._msa = None # multiple sequene alignment object
	self._tree = None # tree object of alignment

	def _init_(msa_filepath, tree_filepath):
		# Read files and unpack into iterable objects
		self._msa = AlignIO.read(msa_filepath, "clustal")
		self._tree = Phylo.read(tree_filepath, "newick")

	def parsimony_score(tree):
		print("score")
