# 'stf' stands for Spatio-Temporal Filter, and implements
# the state merging algorithm a la Shalizi's LICORS method.
# This does estimation for both a transducer (a trivial 
# network), and for an arbitrary network, both of which
# are prescribed by their adjacency matrix.
#
# The 'ensemble' refers to the fact that we compute
# the causal states by 'averaging' over trials at a
# fixed time t rather than 'averaging' over time
# for a fixed trial. This is the more appropriate
# approach when trying to compute how the local
# statistical complexity C(t, v) varies over time.
#
# DMD 291013-13-36

import scipy.stats
import numpy
import copy
from collections import deque
import ipdb
import sys

from stf_methods import *

print '\n\n\nThis code isn\'t actually set up to work with ensembles. You need to modify generate_hist_dict to generate_hist_dict_ensemble and filter_states to filter_states_ensemble.\n\n\n\n'

sys.exit(1)

# numpy.random.seed(1) # Fix the random number generator so we get
					   # reproducible results.

# adj_file = 'edge_list_3K_user_connected_directed.txt'
adj_file = 'adj_mat_toy.txt'

ofile = open(adj_file)

weighted = False

ofile.readline()

line = ofile.readline()

# Node of interest.

noi = '1'

sources = []

if weighted:
	while line != '':
		if '\t' in line:
			source, dest, weight = line.strip().split('\t')
		elif ',' in line:
			source, dest, weight = line.strip().split(',')
		else:
			source, dest, weight = line.strip().split(' ')

		if dest == noi:
			sources.append(source)

		line = ofile.readline()
else:
	while line != '':
		if '\t' in line:
			source, dest = line.strip().split('\t')
		elif ',' in line:
			source, dest = line.strip().split(',')
		else:
			source, dest = line.strip().split(' ')

		if dest == noi:
			sources.append(source)

		line = ofile.readline()

ofile.close()

network_type = 'toy_transducer'
# network_type = 'twitter_p1_i2'

datatype = 'timeseries_synthetic/{}'.format(network_type)

# dataset = 8
# datatype = 'timeseries/NEURO-Set' + str(dataset)

# sources_ts contains all of the time series
# for the inputs *into* a particular node (in
# the case of a transducer) or *adjacent to*
# a particular node (in the case of a spatio-
# temporal random field).

# Get the dimension of the timeseries,
# i.e. (number of samples)x(number of timepoints).

sample_count = 0

with open('{}/sample{}.dat'.format(datatype, noi)) as ofile:
	for line in ofile:
		sample_count += 1

	if ';' in line: # If we have an alphabet of size greater than 10.
		symbol_count = len(line.strip()[:-1].split(';'))
	else: # If we have an alphabet of size less than or equal to 10.
		symbol_count = len(line.strip())

if len(sources) == 0:
	sources_ts = None
else:
	sources_ts = numpy.zeros((sample_count, symbol_count, len(sources)), dtype = 'int16')

	for source_index, source in enumerate(sources):
		with open('{}/sample{}.dat'.format(datatype, source)) as ofile:
			for sample_index, line in enumerate(ofile):
				if ';' in line:
					ts = line.strip()[:-1].split(';')

					sources_ts[sample_index, :, source_index] = cur_ts[sample_index, :] = map(int, ts)
				else:
					ts = line.strip()

					sources_ts[sample_index, :, source_index] = map(int, ts)

# noi_ts contains the time series for the node that
# we wish to predict.

noi_ts = numpy.zeros((sample_count, symbol_count), dtype = 'int16')

with open('{}/sample{}.dat'.format(datatype, noi)) as ofile:
	for sample_index, line in enumerate(ofile):
		if ';' in line:
			ts = line.strip()[:-1].split(';')

			noi_ts[sample_index, :] = map(int, ts)
		else:
			ts = line.strip()

			noi_ts[sample_index, :] = map(int, ts)

L = 1 # The past to consider

alphabet_size = 2 # The number of possible symbols

hist_dict, hists = generate_hist_dict(noi_ts = noi_ts, sources_ts = sources_ts, num_symbols = alphabet_size, L = L)

df = alphabet_size - 1

# Perform the hypothesis test:
#
#	   H0: Two distributions are equal
#	   H1: Two distributions are *not* equal
#
# at a level alpha


alpha = 0.001

test_type = 'chisquared'

states_counts, states_probs, hist_lookup = csmr(hists, alphabet_size, alpha = 0.001, H_test = test_type)

states_final = states_probs

print '\n\n'

for state in states_counts:
	for hist in state[0]:
		print hist

	print state[1]

	print '\n'

print '\n\n'

state_seq = filter_states(noi_ts, sources_ts, hist_lookup, L = L)

state_props = numpy.bincount(state_seq)

state_probs = state_props / float(numpy.sum(state_props))

C = 0

for prob in state_probs:
	C += prob*numpy.log2(prob)

C = -C

print 'The local statistical complexity is {}...'.format(C)