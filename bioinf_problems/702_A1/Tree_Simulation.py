"""
Tree_Simulation.py
Author: Edwardo Reynolds
Date: 28/4/2016
"""

import Tree
import Node
import numpy as np
from scipy.linalg import expm


# Read in parameters from file
# From here, using the inputs, simulate the tree and sequences on the tree
def read_input():
    # read each line, set it to a variable
    f = open('input.txt', 'r')
    n = int(f.readline())
    model = f.readline().strip()
    par = float(f.readline())
    length = int(f.readline())
    pi = map(float, f.readline().split())
    free_vars = map(float, f.readline().split())
    u = float(f.readline())
    filename = f.readline().strip()
    do_sample = f.readline().strip()
    size = f.readline().strip()
    f.close()
    # create Q_hat, beta, and Q
    Q = generate_q(pi, free_vars)

    # call create_tree, create_seq, and sim_tree
    tree = create_tree(n, model, par)
    seq = create_sequence(pi, length)
    tree.get_root().set_sequence(seq)
    root = simulate_seq_on_tree(tree.get_root(), seq, Q, u)

    out = open(filename + ".tree", 'w')
    out.write(root.get_newick())
    out.close()

    out2 = open(filename + ".nex", 'w')
    out2.write("begin trees;\n")
    for leaf in root.get_leaves():
        out2.write(leaf.get_sequence())
        out2.write("\n")
    out2.write("end;\n")
    out2.close()

    # For simulating many 2 taxa trees
    if do_sample == "true":
        sim_single_edges(n, model, par, length, pi, Q, u, filename, int(size))


def generate_q(pi, free_vars):
    # Find Q_hat, and beta
    # Generate normalised rate matrix
    q_hat = np.matrix([[0, free_vars[0]*pi[1], free_vars[1]*pi[2], free_vars[2]*pi[3]],
                       [free_vars[0]*pi[0], 0, free_vars[3]*pi[2], free_vars[4]*pi[3]],
                       [free_vars[1]*pi[0], free_vars[3]*pi[1], 0, free_vars[5]*pi[3]],
                       [free_vars[2]*pi[0], free_vars[4]*pi[1], free_vars[5]*pi[2],0]])
    for i in range(4):
        q_hat[i, i] = 0.0 - np.sum(q_hat[i])

    beta = 1.0/(2*(free_vars[0]*pi[0]*pi[1] + free_vars[1]*pi[0]*pi[2] + free_vars[2]*pi[0]*pi[3] +
                   free_vars[3]*pi[1]*pi[2] + free_vars[4]*pi[1]*pi[3] + free_vars[5]*pi[2]*pi[3]))
    q = beta * q_hat
    return q


def get_random_exponential(k, model, par):
    # calculate tk based on model and parameter
    tk = 0.0

    if model == 'coalescent':
        tk = np.random.exponential(1/(k*(k-1.0)/par))
    elif model == 'yule':
        tk = np.random.exponential(1/(k*par))
    return tk


def create_tree(n, model, par):
    # create a list of all nodes to be coalesced.
    list_of_nodes = []
    for i in range(n):
        node = Node.Node(str(i))
        node.set_height(0.0)
        list_of_nodes.append(node)

    # create tree from this list of nodes
    k = n
    t = 0.0
    label = n
    while k >= 2:
        tk = get_random_exponential(k, model, par)
        t = t + tk
        # choose 2 nodes at random, add a parent of these two nodes, next_node
        lineages_to_coalesce = np.random.choice(list_of_nodes, 2, False)
        next_node = Node.Node(str(label))
        label += 1
        next_node.set_height(t)
        for child in lineages_to_coalesce:
            next_node.add_child(child)
            child.set_parent(next_node)
            list_of_nodes.remove(child)

        list_of_nodes.append(next_node)
        k -= 1

    tree = Tree.Tree()
    tree.set_root(list_of_nodes[0])
    return tree


def create_sequence(pi, length):
    # generate a random sequence of length, length, based on stationary dist, pi.
    seq = ""
    possible_values = [0, 1, 2, 3]
    for i in range(length):
        next_base = np.random.choice(possible_values, 1, False, pi)[0]
        seq += str(next_base)

    return seq


def mutate_sequence(seq, Q, u, t):
    # Takes seq, and mutates it according to Q, t, and u.
    next_seq = ""
    matrix_exp = expm(u*Q*t)
    possible_values = [0, 1, 2, 3]

    for base in seq:
        base = int(base)
        new_base = np.random.choice(possible_values, 1, False, matrix_exp[base])[0]
        next_seq += str(new_base)

    return next_seq


def get_base_sequence(seq):
    # Converts sequence of 0,1,2,3 to a sequence of DNA bases (A, C, G, T)
    base_seq = ""
    possible_bases = ['A', 'C', 'G', 'T']
    for i in seq:
        base_seq += possible_bases[int(i)]

    return base_seq


def simulate_seq_on_tree(node, seq, Q, u):
    # Use sequence from parent, mutate it according to Q, u, and t (edge length).
    # Continue to mutate that sequence until you reach leaves.  -- Uses DFS approach.
    if node.get_sequence() is None:
        t = node.get_parent().get_height() - node.get_height()
        next_seq = mutate_sequence(seq, Q, u, t)
        node.set_sequence(next_seq)

    new_seq = node.get_sequence()
    if not node.is_leaf():
        for child in node.get_children():
            simulate_seq_on_tree(child, new_seq, Q, u)

    base_seq = get_base_sequence(new_seq)
    node.set_sequence(base_seq)

    return node


def sim_single_edges(n, model, par, length, pi, Q, u, filename, size):
    # Simulate 'size' 'n' taxa trees as samples
    # Generate a list of edge_lengths, and a list of equivalent hamming distances
    # Write lists to file to be imported into R

    edge_lengths = []
    hamming_dist = []

    for i in range(size):
        # simulate trees as in the default, and take distances from parent to one of the children.
        tree = create_tree(n, model, par)
        seq = create_sequence(pi, length)
        tree.get_root().set_sequence(seq)
        root = simulate_seq_on_tree(tree.get_root(), seq, Q, u)

        edge_lengths.append(root.get_height() - root.get_children()[0].get_height())
        seq1 = root.get_sequence()
        seq2 = root.get_children()[0].get_sequence()
        correct_bases = 0.0

        for j in range(len(seq1)):
            if seq1[j] == seq2[j]:
                correct_bases += 1.0
        hamming_dist.append(1 - correct_bases/len(seq1))

    out = open(filename + ".sample", 'w')
    for i in range(len(edge_lengths)):
        out.write(str(hamming_dist[i]))
        out.write(" ")
        out.write(str(edge_lengths[i]))
        out.write("\n")

    out.close()

read_input()
