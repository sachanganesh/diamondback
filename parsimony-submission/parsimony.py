"""
Parsimony Tree Improvement 2017
By Sachandhan Ganesh

General Usage:
	parsimony.py -a <alignment> -n <newick_tree> -r <number_of_random_rearrangements> -p <probability>
"""

import argparse
from Bio import AlignIO, Phylo # Format reader for clustal alignments and newick trees
from Bio.Phylo.TreeConstruction import *
from Bio.Phylo.BaseTree import *
from copy import copy, deepcopy

import numpy as np
from matplotlib import pyplot as plt


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
		if clade is not "Clade":
			# search for clade id in msa and return index
			for i in range(len(msa)):
				if msa[i].id == clade:
					return i
		else:
			return None

	@staticmethod
	def _get_terminal_states(msa, tree, nt_pos):
		terminals = tree.get_terminals()
		nts = []

		# get the nucleotide at a position for each leaf
		for clade in terminals:
			i = ParsimonyTree._find_alignment_index_from_clade(msa, clade)
			if i is not None:
				nts.append(msa[i, nt_pos])
			else:
				nts.append(None)

		# combine both nucleotide and leaf data into dict
		states = {}
		for clade, nt in zip(terminals, nts):
			states[clade] = set([nt])

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

	@staticmethod
	def _remove_extraneous_clade(tree):
		for c in tree.find_clades():
			clade_parent = ParsimonyTree._get_parents(tree)
			p = clade_parent[c]

			if len(c.clades) is 1:
				if p is not tree.root and p is not None:
					if p.branch_length is not None and c.branch_length is not None:
						p.branch_length += c.branch_length

					p.clades.extend(c.clades)
					p.clades.remove(c)
				else:
					ch = c.clades[0]

					if c.branch_length is not None and ch.branch_length is not None:
						c.branch_length += ch.branch_length

					c.clades.extend(ch.clades)
					c.clades.remove(ch)
			elif c.is_terminal() and c.name == None:
				p.clades.remove(c)

	@staticmethod
	def _has_subtree(clade):
		return not all([c.is_terminal() for c in clade.clades])

	@staticmethod
	def _neighbor_visited(clade, other_clade, visited_pairs):
		for pair in visited_pairs:
			if pair[0] == clade and pair[1] == other_clade:
				return True

		return False

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
				common_states = set()

				for i, child in enumerate(children):
					if i is 0:
						common_states = clade_states[child]
					else:
						beta = clade_states[child]

						# find common nucleotide states between children
						common_states = common_states.intersection(beta)

				# if there are common states, use them
				if common_states:
					clade_states[clade] = common_states
				else:
					for child in children:
						common_states = common_states.union(clade_states[child])

					# take all states and increment parsimony score
					clade_states[clade] = common_states
					score += 1

			# add to parsimony score
			total_score += score

		return total_score

	@staticmethod
	def get_nni_neighbors(tree):
		neighbors = []
		skippable = []
		node_parent = ParsimonyTree._get_parents(tree)

		# Assuming a binary tree:
		# for every internal node
		for clade in tree.get_nonterminals():
			# if the clade is an internal node that's not the root or its direct children
			if clade is not tree.root and clade not in skippable:
				parent = node_parent[clade]
				left_child = clade.clades[0]
				right_child = clade.clades[1]

				# depending on if the clade is a left or right node
				# switch the node with the sibling branch
				# then alternate
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

				# as long as the root's children aren't leaves
				# switch all the children of each branch
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
	def get_spr_neighbors(tree):
		neighbors = []
		visited_pairs = []
		clade_parent = ParsimonyTree._get_parents(tree)

		# for every node in the tree
		for clade in tree.find_clades():
			# as long as it's not the root
			if clade is not tree.root:
				# get the parent and siblings of the node
				parent = clade_parent[clade]
				siblings = copy(parent.clades)

				# get the original index and remove the node from the tree
				clade_ind = parent.clades.index(clade)
				parent.clades.remove(clade)

				# for every other node in the tree
				for other_clade in tree.find_clades():
					# as long as it hasn't been visited before and it's not the root or sibling
					if not ParsimonyTree._neighbor_visited(clade, other_clade, visited_pairs) and other_clade is not tree.root and other_clade not in siblings:
						if other_clade is parent and len(siblings) - 1 is 1:
							continue
						else:
							# get the index and save this neighbor-tree pair
							ind = clade_parent[other_clade].clades.index(other_clade)
							visited_pairs.append((clade, other_clade))

							# insert the branch into this place in the tree
							ParsimonyTree._insert_bifurcation(clade_parent[other_clade], ind, {}, clade)

							# copy the tree
							cp_tree = deepcopy(tree)

							# clean up extraneous clades
							ParsimonyTree._remove_extraneous_clade(cp_tree)

							# add tree to neighbors
							neighbors.append(cp_tree)

							# undo insertion of the branch
							ParsimonyTree._remove_bifurcation(clade_parent[other_clade], clade)

				# put the node back where it was
				parent.clades.insert(clade_ind, clade)

		return neighbors

	@staticmethod
	def get_tbr_neighbors(tree):
		neighbors = []
		visited_pairs = []
		clade_parent = ParsimonyTree._get_parents(tree)

		for clade in tree.find_clades():
			if clade is not tree.root:
				parent = clade_parent[clade]
				siblings = copy(parent.clades)

				clade_ind = parent.clades.index(clade)
				parent.clades.remove(clade)

				tmp = Tree.from_clade(deepcopy(clade), rooted=True)

				for other_clade in tree.find_clades():
					if not ParsimonyTree._neighbor_visited(tmp, other_clade, visited_pairs) and other_clade is not tree.root and other_clade not in siblings:
						if other_clade is parent and len(siblings) - 1 is 1:
							continue
						else:
							for order in tmp.find_clades():
								# tmp.root_with_outgroup(order)

								ind = clade_parent[other_clade].clades.index(other_clade)
								visited_pairs.append((tmp.root, other_clade))

								ParsimonyTree._insert_bifurcation(clade_parent[other_clade], ind, {}, tmp.root)

								cp_tree = deepcopy(tree)

								# clean up extraneous clades
								ParsimonyTree._remove_extraneous_clade(cp_tree)

								neighbors.append(cp_tree)

								ParsimonyTree._remove_bifurcation(clade_parent[other_clade], tmp.root)

				parent.clades.insert(clade_ind, clade)

		return neighbors

	@staticmethod
	def visualize_tree(tree, figure_title="Phylogenetic Tree"):
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


def parse_arguments():
	parser = argparse.ArgumentParser(description="Mossback: a parsimony tree searcher using Metropolis sampling")
	parser.add_argument("a", metavar="msa_filepath", type=str, help="Multiple sequence alignment file in Clustal format")
	parser.add_argument("n", metavar="tree_filepath", type=str, help="Tree representation file in Newick format")
	parser.add_argument("r", metavar="r", type=int, default=50, help="Number of random local tree rearrangements")
	parser.add_argument("p", metavar="p", type=float, default=0.05, help="Probability at which worse parsimony scores will be accepted")
	parser.add_argument("o", metavar="output_filepath", type=str, help="Filepath to write final tree found through search")
	parser.add_argument("-nni", "--nni", action="store_true", help="Use nearest neighbor interchange for neighbor-tree creation")
	parser.add_argument("-spr", "--spr", action="store_true", help="Use subtree pruning and regrafting for neighbor-tree creation")
	parser.add_argument("-tbr", "--tbr", action="store_true", help="Use tree bisection and reconnection for neighbor-tree creation")

	args = parser.parse_args()
	return args


def main():
	args = parse_arguments()

	msa = ParsimonyTree.read_msa(args.a)
	i_tree = ParsimonyTree.read_tree(args.n)

	nb_f = ParsimonyTree.get_nni_neighbors
	if args.spr:
		nb_f = ParsimonyTree.get_spr_neighbors
	elif args.tbr:
		nb_f = ParsimonyTree.get_tbr_neighbors

	mcmc = MonteCarlo(msa, i_tree, nb_f, args.r, args.p)
	f_tree = mcmc.get_tree()

	with open(args.o, "w") as outfile:
		Phylo.write(f_tree, outfile, "newick")

	print("\n=========================\n")
	print("Original Tree")
	print("Score:", ParsimonyTree.get_parsimony_score(msa, i_tree))
	Phylo.draw(i_tree)
	Phylo.draw_ascii(i_tree)

	print("\n=========================\n")
	print("Final Tree")
	print("Score:", ParsimonyTree.get_parsimony_score(msa, f_tree))
	Phylo.draw(f_tree)
	Phylo.draw_ascii(f_tree)

	print("\n=========================\n")
	print("Histogram of Parsimony Scores")
	plt.title("Histogram of Parsimony Scores")
	plt.hist(mcmc.get_scores())
	plt.show()


if __name__ == "__main__":
	main()
