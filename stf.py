# 'stf' stands for Spatio-Temporal Filter, and implements
# the state merging algorithm a la Shalizi's LICORS method.
# This does estimation for both a transducer (a trivial 
# network), and for an arbitrary network, both of which
# are prescribed by their adjacency matrix.
#
# DMD 030713-11-40

import scipy.stats
import numpy
import copy
from collections import deque
import ipdb

from stf_methods import *

# numpy.random.seed(1) # Fix the random number generator so we get
					   # reproducible results.

adj_file = 'edge_list_3K_user_connected_directed.txt'
# adj_file = 'adj_mat_toy.txt'

ofile = open(adj_file)

weighted = False

ofile.readline()

line = ofile.readline()

# Node of interest.

noi = '17'

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

# sources_ts contains all of the time series
# for the inputs *into* a particular node (in
# the case of a transducer) or *adjacent to*
# a particular node (in the case of a spatio-
# temporal random field).

sources_ts = []

datatype = 'timeseries_synthetic/twitter_p1_i2'
# datatype = 'timeseries_synthetic/toy_transducer'

# dataset = 8
# datatype = 'timeseries/NEURO-Set' + str(dataset)

for source in sources:
	ofile = open('{}/sample{}.dat'.format(datatype, source))

	sources_ts.append(ofile.readline().rstrip('\n'))

	ofile.close()

ofile = open('{}/sample{}.dat'.format(datatype, noi))

# noi_ts contains the time series for the node that
# we wish to predict.

noi_ts = ofile.readline().rstrip('\n')

n = len(noi_ts)

ofile.close()

L = 1 # The past to consider

num_symbols = 2 # The number of possible symbols

hist_dict, hists = generate_hist_dict(noi_ts = noi_ts, sources_ts = sources_ts, num_symbols = num_symbols, L = L)

df = num_symbols - 1

# Perform the hypothesis test:
#
#	   H0: Two distributions are equal
#	   H1: Two distributions are *not* equal
#
# at a level alpha


alpha = 0.001

test_type = 'chisquared'

states_counts, states_probs, hist_lookup = csmr(hists, num_symbols, alpha = 0.001, H_test = test_type)

states_final = states_probs

print states_final

state_seq = filter_states(noi_ts, sources_ts, hist_lookup, L = L)

state_props = numpy.bincount(state_seq)

state_probs = state_props / float(numpy.sum(state_props))

C = 0

for prob in state_probs:
	C += prob*numpy.log2(prob)

C = -C

print 'The local statistical complexity is {}...'.format(C)