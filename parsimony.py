"""
Parsimony Tree Improvement 2017
By Sachandhan Ganesh

General Usage:
	parsimony.py -a <alignment> -n <newick_tree> -r <number_of_random_rearrangements> -p <probability>
"""

from Bio import AlignIO, Phylo # Format reader for clustal alignments and newick trees
from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.BaseTree import *
from ete3 import Tree # Newick tree reader
from copy import copy, deepcopy

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns


class ParsimonyTree(object):
	@staticmethod
	def read_msa(msa_filepath):
		return AlignIO.read(msa_filepath, "clustal") # multiple sequene alignment object

	@staticmethod
	def read_tree(tree_filepath):
		return Phylo.read(tree_filepath, "newick") # tree object of alignment

	@staticmethod
	def _find_alignment_index_from_clade(msa, clade):
		# get clade id
		clade = str(clade)

		# search for clade id in msa and return index
		for i in range(len(msa)):
			if msa[i].id == clade:
				return i

	@staticmethod
	def _get_terminal_states(msa, tree, nt_pos):
		terminals = tree.get_terminals()
		nts = []

		# get the nucleotide at a position for each leaf
		for clade in terminals:
			i = ParsimonyTree._find_alignment_index_from_clade(msa, clade)
			nts.append(msa[i, nt_pos])

		# combine both nucleotide and leaf data into dict
		states = {}
		for clade, nt in zip(terminals, nts):
			states[clade] = set(nt)

		return states

	@staticmethod
	def _get_parents(tree):
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

	@staticmethod
	def _remove_from_tree(tree, removable):
		clade_parent = ParsimonyTree._get_parents(tree)
		rem_parent = clade_parent[removable]

		rem_ind = rem_parent.clades.index(removable)
		rem_parent.clades.remove(removable)

		if rem_parent is not tree.root and len(rem_parent.clades) == 1:
			rem_sibling = rem_parent.clades[0]
			extra_parent = clade_parent[rem_parent]

			extra_ind = extra_parent.clades.index(rem_parent)
			extra_parent.clades.remove(rem_parent)

			extra_parent.clades.insert(extra_ind, rem_sibling)

			return (rem_ind, removable), (extra_ind, rem_parent)

		return (rem_ind, removable), None

	def _undo_remove_from_tree(tree, parent, clade_dt, pruned_dt):
		if pruned_dt is None:


	@staticmethod
	def _insert_bifurcation(parent, clade_pos, new_clade_attr, insertable):
		c = Clade(**new_clade_attr)

		c.clades.append(insertable)
		c.clades.append(parent.clades.pop(clade_pos))
		parent.clades.insert(clade_pos, c)

	@staticmethod
	def _remove_bifurcation(clade, removable):
		clade_parent = ParsimonyTree._get_parents(clade)
		parent = clade_parent[removable]

		if parent.clades[0] == removable:
			ind = 1
			keep = parent.clades[1]
		else:
			ind = 0
			keep = parent.clades[0]

		clade.clades.remove(parent)
		clade.clades.insert(ind, keep)
		del parent

	@staticmethod
	def _has_subtree(clade):
		return not all([c.is_terminal() for c in clade.clades])

	@staticmethod
	def get_parsimony_score(msa, tree):
		total_score = 0

		# for each nucleotide position
		for i in range(msa.get_alignment_length()):
			score = 0 # score for tree given nucleotide states

			# Assign states for leaves
			clade_states = ParsimonyTree._get_terminal_states(msa, tree, i)

			# postorder traversal of internal nodes to resolve commonalities
			for clade in tree.get_nonterminals(order="postorder"):
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

	@staticmethod
	def get_nni_neighbors(tree):
		neighbors = []
		skippable = []
		node_parent = ParsimonyTree._get_parents(tree)

		for clade in tree.get_nonterminals():
			if clade is not tree.root and clade not in skippable:
				parent = node_parent[clade]
				left_child = clade.clades[0]
				right_child = clade.clades[1]

				if clade == parent.clades[0]:
					sibling = parent.clades[1]

					parent.clades.remove(parent.clades[1])
					clade.clades.remove(clade.clades[0])

					parent.clades.append(left_child)
					clade.clades.append(sibling)

					neighbors.append(deepcopy(tree))

					parent.clades.remove(parent.clades[1])
					clade.clades.remove(clade.clades[0])

					parent.clades.append(right_child)
					clade.clades.append(left_child)

					neighbors.append(deepcopy(tree))

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

					neighbors.append(deepcopy(tree))

					clade.clades.remove(clade.clades[0])
					parent.clades.remove(parent.clades[0])

					clade.clades.append(left_child)
					parent.clades.insert(0, right_child)

					neighbors.append(deepcopy(tree))

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

					neighbors.append(deepcopy(tree))

					left_child.clades.remove(left_child.clades[1])
					right_child.clades.remove(right_child.clades[0])

					left_child.clades.append(right_child_right)
					right_child.clades.append(right_child_left)

					neighbors.append(deepcopy(tree))

					# revert back
					left_child.clades.remove(left_child.clades[1])
					right_child.clades.remove(right_child.clades[0])

					left_child.clades.append(left_child_right)
					right_child.clades.append(right_child_right)

		return neighbors

	@staticmethod
	def get_spr_neighbors(msa, tree):
		print(tree)
		neighbors = []
		clade_parent = ParsimonyTree._get_parents(tree)

		for clade in tree.get_nonterminals():
			if clade is not tree.root and ParsimonyTree._has_subtree(clade):
				parent = clade_parent[clade]
				siblings = copy(parent.clades)

				clade_dt, pruned_dt = ParsimonyTree._remove_from_tree(tree, clade)
				print("pruned tree:")
				print(tree)
				for other_clade in tree.find_clades():
					if other_clade is not tree.root and other_clade not in siblings:
						ind = clade_parent[other_clade].clades.index(other_clade)
						ParsimonyTree._insert_bifurcation(clade_parent[other_clade], ind, {}, clade)

						neighbors.append(tree)
						Phylo.draw_ascii(tree)
						print("neighbor:")
						print(tree)
						# print(ParsimonyTree.get_parsimony_score(msa, tree))

						ParsimonyTree._remove_bifurcation(clade_parent[other_clade], clade)
						exit(1)

				parent.clades.insert(clade_dt[0], clade)

		return neighbors

	@staticmethod
	def get_tbr_neighbors(tree):
		pass

	@staticmethod
	def visualize_tree(tree):
		Phylo.draw(tree)


class MonteCarlo(object):
	def __init__(self, msa, init_tree, neighbor_fn, num_iterations, accept_prob):
		self._msa = msa
		self._states = [init_tree]
		self._past_scores = []
		self._neighbor_fn = neighbor_fn
		self._num_iterations = num_iterations
		self._accept_prob = accept_prob

	def _update(self, state):
		self._states.append(state)

	def _get_next_state(self):
		state = self._states[-1]
		current_score = ParsimonyTree.get_parsimony_score(self._msa, state)
		self._past_scores.append(current_score)

		neighbors = self._neighbor_fn(state)
		scores = []
		n_high = 0

		if len(neighbors) > 1:
			for n in neighbors:
				s = ParsimonyTree.get_parsimony_score(self._msa, n)
				scores.append(s)

				if s > current_score:
					n_high += 1

			if n_high == 0:
				high_w = 0
			elif n_high == len(scores):
				high_w = 1 / n_high
			else:
				high_w = self._accept_prob / n_high

			unif_w = (1 - high_w * n_high) / (len(neighbors) - n_high)
			w = [high_w if s > current_score else unif_w for s in scores]

			self._past_scores.extend(scores)
			return np.random.choice(neighbors, p=w)
		else:
			return None

	def get_tree(self):
		next_state = self._get_next_state()

		for i in range(self._num_iterations):
			if next_state is not None:
				self._update(next_state)
				next_state = self._get_next_state()
			else:
				return self._states[-1]

		return next_state

	def get_scores(self):
		return self._past_scores


def main():
	msa = ParsimonyTree.read_msa("./test_data/test_msa.txt")
	i_tree = ParsimonyTree.read_tree("./test_data/test_tree.txt")
	Phylo.draw_ascii(i_tree)

	ParsimonyTree.get_spr_neighbors(msa, i_tree)

	# mcmc = MonteCarlo(msa, i_tree, ParsimonyTree.get_spr_neighbors, 20, 0.1)
	# f_tree = mcmc.get_tree()
	#
	# Phylo.draw_ascii(i_tree)
	# print("inital score:", ParsimonyTree.get_parsimony_score(msa, i_tree))
	# Phylo.draw_ascii(f_tree)
	# print("final score:", ParsimonyTree.get_parsimony_score(msa, f_tree))
	#
	# sns.distplot(mcmc.get_scores())
	# plt.show()


if __name__ == "__main__":
	main()
