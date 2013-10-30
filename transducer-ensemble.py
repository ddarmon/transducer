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
# DMD 301013-13-23

import scipy.stats
import numpy
import copy
from collections import deque
import ipdb

from stf_methods import *

# numpy.random.seed(1) # Fix the random number generator so we get
					   # reproducible results.

print '\n\n\n\n\n\nWARNING!!! This does not currently handle input symbols with more than 0..9 ... \n\n\n\n\n'

adj_file = 'edge_list_3K_user_connected_directed.txt'
# adj_file = 'adj_mat_toy.txt'

ofile = open(adj_file)

weighted = False

ofile.readline()

line = ofile.readline()

# Node of interest.

noi = '1'

# sources_ts contains all of the time series
# for the inputs *into* a particular node (in
# the case of a transducer) or *adjacent to*
# a particular node (in the case of a spatio-
# temporal random field).

# datatype = 'timeseries_synthetic/twitter_p1_i2'
datatype = 'timeseries_synthetic/toy_transducer'

source = '{}/input{}'.format(datatype, noi)

sample_count = 0

with open('{}.dat'.format(source)) as ofile:
	for line in ofile:
		sample_count += 1

	if ';' in line: # If we have an alphabet of size greater than 10.
		symbol_count = len(line.strip()[:-1].split(';'))
	else: # If we have an alphabet of size less than or equal to 10.
		symbol_count = len(line.strip())

sources_ts = numpy.zeros((sample_count, symbol_count, 1), dtype = 'int8')

with open('{}.dat'.format(source)) as ofile:
	for sample_index, line in enumerate(ofile):
		if ';' in line:
			ts = line.strip()[:-1].split(';')

			sources_ts[sample_index, :, 0] = map(int, ts)
		else:
			ts = line.strip()

			sources_ts[sample_index, :, 0] = numpy.fromstring(ts, dtype = 'int8') - 48


# noi_ts contains the time series for the node that
# we wish to predict.

noi_ts = numpy.zeros((sample_count, symbol_count), dtype = 'int8')

with open('{}/sample{}.dat'.format(datatype, noi)) as ofile:
	for sample_index, line in enumerate(ofile):
		if ';' in line:
			ts = line.strip()[:-1].split(';')

			noi_ts[sample_index, :] = map(int, ts)
		else:
			ts = line.strip()

			noi_ts[sample_index, :] = numpy.fromstring(ts, dtype = 'int8') - 48

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

# print states_final

state_seq = filter_states(noi_ts, sources_ts, hist_lookup, L = L)

state_props = numpy.bincount(state_seq)

state_probs = state_props / float(numpy.sum(state_props))

C = 0

for prob in state_probs:
	C += prob*numpy.log2(prob)

C = -C

print 'The local statistical complexity is {}...'.format(C)

# print hist_dict

# for hist in hist_dict:
# 	count_0, count_1 = hist_dict[hist]

# 	print hist, count_1/float(count_0 + count_1), count_0, count_1