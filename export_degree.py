# Export the in / out degrees of 
# the vertices in a network.
# 
# DMD, 011113-13-38

import numpy
import pylab
import ipdb

from network_methods import *

# The weighted adjacency matrix, stored as a dictionary of the form:
# {from_node : [(to_node1, weight1), (to_node2, weight2)]}

# adj_matrix_out = {0 : [(1, uniform_weight), (2, uniform_weight), (3, uniform_weight), (4, uniform_weight), (5, 0.005)], 5 : [(6, uniform_weight), (7, uniform_weight), (8, uniform_weight), (0, 0.005)]}

# network_type = 'toy_transducer'
network_type = 'twitter_p1_i2_ensemble'

# adj_mat_fname = 'adj_mat_toy.txt'
adj_mat_fname = 'edge_list_3K_user_connected_directed.txt'

adj_matrix_out, adj_matrix_in, Nv = get_adjacency_to_weight(adj_mat_fname)

with open('{}_indegree.dat'.format(network_type), 'w') as wfile:
	for nodeid in range(0, len(adj_matrix_in)):
		wfile.write('{}\t{}\n'.format(nodeid, len(adj_matrix_in[nodeid])))