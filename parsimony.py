"""
Parsimony Tree Improvement 2017
By Sachandhan Ganesh

General Usage:
	parsimony.py -a <alignment> -n <newick_tree> -r <number_of_random_rearrangements> -p <probability>
"""

from Bio import AlignIO, Phylo # Format reader for clustal alignments and newick trees
from Bio.Phylo.TreeConstruction import *
from ete3 import Tree # Newick tree reader

import numpy as np
from copy import deepcopy


class ParsimonyTree(object):
	def __init__(self, msa_filepath, tree_filepath):
		# Read files and unpack into iterable objects
		self._msa = AlignIO.read(msa_filepath, "clustal") # multiple sequene alignment object
		self._tree = Phylo.read(tree_filepath, "newick") # tree object of alignment

		self._num_rearrangements = 1000


	def get_parsimony_score(self):
		total_score = 0

		# for each nucleotide position
		for i in range(self._msa.get_alignment_length()):
			score = 0 # score for tree given nucleotide states

			# Assign states for leaves
			clade_states = self._get_terminal_states(i)

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

			# add to parsimony score
			total_score += score

		return total_score


	def _find_alignment_index_from_clade(self, clade):
		# get clade id
		clade = str(clade)

		# search for clade id in msa and return index
		for i in range(len(self._msa)):
			if self._msa[i].id == clade:
				return i


	def _get_terminal_states(self, nt_pos):
		terminals = self._tree.get_terminals()
		nts = []

		# get the nucleotide at a position for each leaf
		for clade in terminals:
			i = self._find_alignment_index_from_clade(clade)
			nts.append(self._msa[i, nt_pos])

		# combine both nucleotide and leaf data into dict
		states = {}
		for clade, nt in zip(terminals, nts):
			states[clade] = set(nt)

		return states


	def _get_parents(self, tree):
		parents = {}

		for clade in tree.find_clades():
			path = tree.get_path(clade)

			if len(path) == 0:
				parents[clade] = None
			elif len(path) == 1:
				parents[clade] = tree.root
			else:
				parents[clade] = path[-2]

		return parents


	def _get_neighbors(self):
		neighbors = []
		skippable = []
		node_parent = self._get_parents(self._tree)

		for clade in self._tree.get_nonterminals():
			if clade is not self._tree.root and clade not in skippable:
				parent = node_parent[clade]
				left_child = clade.clades[0]
				right_child = clade.clades[1]

				if clade == parent.clades[0]:
					sibling = parent.clades[1]

					parent.clades.remove(parent.clades[1])
					clade.clades.remove(clade.clades[0])

					parent.clades.append(left_child)
					clade.clades.append(sibling)

					neighbors.append(deepcopy(self._tree))

					parent.clades.remove(parent.clades[1])
					clade.clades.remove(clade.clades[0])

					parent.clades.append(right_child)
					clade.clades.append(left_child)

					neighbors.append(deepcopy(self._tree))

					# revert back
					parent.clades.remove(parent.clades[1])
					clade.clades.remove(clade.clades[0])

					parent.clades.append(sibling)
					clade.clades.append(right_child)
				else:
					sibling = parent.clades[0]

					clade.clades.remove(clade.clades[0])
					parent.clades.remove(parent.clades[0])

					clade.clades.append(sibling)
					parent.clades.insert(0, left_child)

					neighbors.append(deepcopy(self._tree))

					clade.clades.remove(clade.clades[0])
					parent.clades.remove(parent.clades[0])

					clade.clades.append(left_child)
					parent.clades.insert(0, right_child)

					neighbors.append(deepcopy(self._tree))

					# revert back
					clade.clades.remove(clade.clades[0])
					parent.clades.remove(parent.clades[0])

					clade.clades.append(right_child)
					parent.clades.insert(0, sibling)
			elif clade in skippable:
				# do nothing
				pass
			else:
				left_child = clade.clades[0]
				right_child = clade.clades[1]

				skippable.extend([left_child, right_child])

				if not left_child.is_terminal() and not right_child.is_terminal():
					left_child_right = left_child.clades[1]
					right_child_left = right_child.clades[0]
					right_child_right = right_child.clades[1]

					left_child.clades.remove(left_child.clades[1])
					right_child.clades.remove(right_child.clades[0])

					left_child.clades.append(right_child_left)
					right_child.clades.append(left_child_right)

					neighbors.append(deepcopy(self._tree))

					left_child.clades.remove(left_child.clades[1])
					right_child.clades.remove(right_child.clades[0])

					left_child.clades.append(right_child_right)
					right_child.clades.append(right_child_left)

					neighbors.append(deepcopy(self._tree))

					# revert back
					left_child.clades.remove(left_child.clades[1])
					right_child.clades.remove(right_child.clades[0])

					left_child.clades.append(left_child_right)
					right_child.clades.append(right_child_right)

		return neighbors


		def


	def visialize_tree(self):
		Phylo.draw(self._tree)


def main():
	pt = ParsimonyTree("./test_data/test_msa.txt", "./test_data/test_tree.txt")
	score = pt.get_parsimony_score()
	print(score)

	print("Original")
	Phylo.draw_ascii(pt._tree)

	neighbors = pt.get_neighbors()
	print(len(neighbors))

if __name__ == "__main__":
	main()
