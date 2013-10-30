# Generate a realization from coupled
# inhomogeneous Bernoulli processes
# with a prescribed weighted, directed 
# adjacency matrix. 

import numpy
import pylab
import ipdb

from network_methods import *

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Various parameters:

# The amplitude of the base rate, from
# 0 to lam_max.

p_max = 0.1

# Loop until we've reached the appropriate amount of 
# time.

T_final = 1000

# Loop until we've run the appropriate number of trials.

num_trials = 1000

# The weighted adjacency matrix, stored as a dictionary of the form:
# {from_node : [(to_node1, weight1), (to_node2, weight2)]}

# adj_matrix_out = {0 : [(1, uniform_weight), (2, uniform_weight), (3, uniform_weight), (4, uniform_weight), (5, 0.005)], 5 : [(6, uniform_weight), (7, uniform_weight), (8, uniform_weight), (0, 0.005)]}

network_type = 'toy_transducer'

adj_mat_fname = 'adj_mat_toy.txt'
# adj_mat_fname = 'edge_list_3K_user_connected_directed.txt'

adj_matrix_out, adj_matrix_in, Nv = get_adjacency_to_weight(adj_mat_fname)

scale_weights = 1 # How much to scale the influence weights by

#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def p(t):
	return p_max*numpy.ones(t.shape)

# Generate a bunch of empty files.

for ind in range(Nv):
		ofile = open('timeseries_synthetic/{}/sample{}.dat'.format(network_type, ind), 'w')
		ofile.write('')
		ofile.close()

for trial_ind in range(num_trials):
	if (trial_ind % 10) == 0:
		print 'On trial {}...'.format(trial_ind)
	# The baseline rate, before accounting for any spikes.

	Pt = p(numpy.tile(numpy.arange(T_final), (Nv, 1)))

	# Scale the baseline rate so that each individual
	# vertex has its own base rate.

	pv = numpy.array([5, 1, 5, 5, 5, 5])

	Pt = Pt*pv[:, numpy.newaxis]

	# A placeholder for all of the spikes. It is a 
	# Nv by T_final array.

	Xt = numpy.zeros((Nv, T_final))

	# Us is all of the random draws we will need.

	Us = numpy.random.rand(Nv, T_final)

	# The skeleton to add for each observed spike.

	# impulse_skeleton = numpy.concatenate((numpy.ones(1), numpy.power(numpy.arange(1., 50.), -3)))
	impulse_skeleton = numpy.ones(1)

	for t in range(0, T_final):
		if (t % 10) == 0:
			print 'At timestep {}...'.format(t)

		# Determine which nodes were active.

		active_bool = (Us[:, t] < Pt[:, t])

		active_inds  = numpy.arange(Nv)[active_bool]

		# Update the matrix of states.

		Xt[active_inds, t] = 1

		# Update Pt based on which nodes spikes.
		# To do this, we look at the weighted
		# adjacency matrix and only add a
		# contribution to the appropriate nodes.

		for from_ind in active_inds:
			neighbors = adj_matrix_out.get(from_ind, [])

			# Account for when impulse_skeleton is longer
			# than the remaining timeseries.

			diff_T = T_final - t

			if diff_T <= len(impulse_skeleton):
				amount_forward = diff_T
			else:
				amount_forward = len(impulse_skeleton)+1


			for to_ind, to_weight in neighbors:
				Pt[to_ind, t+1:t+amount_forward] += scale_weights*to_weight*impulse_skeleton[:amount_forward-1]

	# Save the spikes to output files.

	for ind in range(Nv):
		if (ind % 10) == 0:
			print 'Saving data for vertex {}...'.format(ind)
		ofile = open('timeseries_synthetic/{}/sample{}.dat'.format(network_type, ind), 'a')
		for spike in Xt[ind, :]:
			ofile.write('{}'.format(int(spike)))
		ofile.write('\n')
		ofile.close()

	# Plot a few figures.

	# num_plot = Nv
	# num_plot = 5

	# f, axarr = pylab.subplots(num_plot, sharex = True)

	# for axind in range(num_plot):
	# 	axarr[axind].vlines(numpy.arange(T_final)[Xt[axind, :] == 1], ymin = -0.5, ymax = 0.5)
	# 	axarr[axind].yaxis.set_visible(False)

	# f, axarr = pylab.subplots(Nv, sharex = True)

	# for axind in range(num_plot):
	# 	axarr[axind].plot(numpy.arange(T_final), Pt[axind, :])
	# 	axarr[axind].yaxis.set_visible(False)

	# pylab.show()

print '\n\nBe sure that if you\'re running trials using the transducer code you convert the output to a transducer using parse_to_transducer.py.\n\n'