import collections
import numpy

def get_adjacency_weighted(fname):
	# Get the directed, weighted
	# adjacency matrix from a file. Also get
	# out the number of vertices.

	adj_matrix_out = collections.defaultdict(list)
	adj_matrix_in = collections.defaultdict(list)

	Nv = 0

	ofile = open(fname)

	ofile.readline() # To take care of the header

	line = ofile.readline()

	while line != '':
		if ',' in line:
			from_ind, to_ind, weight = line.strip().split(',')
		elif '\t' in line:
			from_ind, to_ind, weight = line.strip().split('\t')
		else:
			from_ind, to_ind, weight = line.strip().split(' ')

		# adj_matrix_out takes node_id and outputs the nodes that
		# node points to.

		# adj_matrix_in  takes node_id and outputs the nodes that
		# point *to* that node.

		adj_matrix_out[int(from_ind)].append([int(to_ind), float(weight)])
		adj_matrix_in[int(to_ind)].append([int(from_ind), float(weight)])

		line = ofile.readline()

		Nv = numpy.max((int(from_ind), int(to_ind), Nv))

	ofile.close()

	# Since we started labeling nodes at 0.

	Nv = Nv + 1

	return adj_matrix_out, adj_matrix_in, Nv

def get_adjacency_to_weight(fname):
	# Get the directed, weighted
	# adjacency matrix from a file. Also get
	# out the number of vertices.

	adj_matrix_out = collections.defaultdict(list)
	adj_matrix_in = collections.defaultdict(list)

	Nv = 0

	ofile = open(fname)

	ofile.readline() # To take care of the header

	line = ofile.readline()

	while line != '':
		if ',' in line:
			from_ind, to_ind = line.strip().split(',')
		elif '\t' in line:
			from_ind, to_ind = line.strip().split('\t')
		else:
			from_ind, to_ind = line.strip().split(' ')

		weight = 0.

		# adj_matrix_out takes node_id and outputs the nodes that
		# node points to.

		# adj_matrix_in  takes node_id and outputs the nodes that
		# point *to* that node.

		adj_matrix_out[int(from_ind)].append([int(to_ind), float(weight)])
		adj_matrix_in[int(to_ind)].append(from_ind)

		line = ofile.readline()

		Nv = numpy.max((int(from_ind), int(to_ind), Nv))

	# After we've read in the adjacency matrix, we need to weight
	# each edge appropriately.

	# For now, let k^{in}_{j} be the in-degree of node j. Then w_{ij}
	# (the weight on the directed edge *from* i *to* j) is
	# 	w_{ij} = 1/k^{in}_{j},
	# based on some sort of attention-type rule.

	for incoming_node in adj_matrix_out: # This gives me a list of the edges that incoming_node attaches to
		for ind in range(len(adj_matrix_out[incoming_node])):
			connecting_node = adj_matrix_in[adj_matrix_out[incoming_node][ind][0]]
			adj_matrix_out[incoming_node][ind][1] = 1./float(len(connecting_node))

	ofile.close()

	# Since we started labeling nodes at 0.

	Nv = Nv + 1

	return adj_matrix_out, adj_matrix_in, Nv