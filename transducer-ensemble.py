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
import pylab

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

noi = '0'

# sources_ts contains all of the time series
# for the inputs *into* a particular node (in
# the case of a transducer) or *adjacent to*
# a particular node (in the case of a spatio-
# temporal random field).

datatype = 'timeseries_synthetic/twitter_p1_i2_ensemble'
# datatype = 'timeseries_synthetic/toy_transducer'

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
		if sample_index == sample_count:
			break # We have more samples for noi than for the input
				  # so we stop.

		if ';' in line:
			ts = line.strip()[:-1].split(';')

			noi_ts[sample_index, :] = map(int, ts)
		else:
			ts = line.strip()

			noi_ts[sample_index, :] = numpy.fromstring(ts, dtype = 'int8') - 48

L = 1 # The past to consider

alphabet_size = 2 # The number of possible symbols

# hist_dict, hists = generate_hist_dict(noi_ts = noi_ts, sources_ts = sources_ts, num_symbols = alphabet_size, L = L)

ts = range(L, noi_ts.shape[1])
# ts = range(L, L+1)

Cs = numpy.zeros(len(ts))

Cs_correct_model = numpy.zeros(len(ts))

for t_ind, t in enumerate(ts):
	print 'At time index {} of {}...'.format(t, ts[-1])

	hist_dict, hists = generate_hist_dict_ensemble(noi_ts = noi_ts[:, (t-L):(t+1)], sources_ts = sources_ts[:, (t-L):t, :], num_symbols = alphabet_size)

	df = alphabet_size - 1

	# Perform the hypothesis test:
	#
	#	   H0: Two distributions are equal
	#	   H1: Two distributions are *not* equal
	#
	# at a level alpha

	alpha = 0.001

	# alpha = 0.001 # We set the size of the test, an upper bound on the
				  #	probability of rejecting that two histories are in
				  # the same equivalence class when they actually *are*.
				  # That is, we fix an allowable Type I error.
	
	# Note: Setting alpha also fixes the Type II error, but 
	# *does not* control it. That is, fixing alpha also
	# fixes the probability that we *do not* reject that
	# two states are in the same equivalence class, even
	# when they *are in two* different equivalence classes.
	#
	# Thus, we have to balance between how often we don't
	# merge when we should (which is controlled by alpha),
	# and how much we do merge when we shouldn't (which
	# is only controlled indirectly by alpha).
	#
	# For *too small* alpha, we'll end up merging too
	# much. For *too large* alpha, we won't end up 
	# merging enough.
	
	# Note: The test statistic associated with the chi-squared
	# test is is only chi-squared *asymptotically*, since the
	# sampling distribution result depends on the 
	# central limit theorem. Thus, for small
	# sample sizes, it is better to use Fisher's exact 
	# test, which uses a test statistic whose sampling 
	# distribution is known *exactly*, even for finite n.

	# test_type = 'chisquared'
	test_type = 'exact'

	states_counts, states_probs, hist_lookup = csmr(hists, alphabet_size, alpha = alpha, H_test = test_type)

	states_final = states_probs

	# print states_final

	# Computing LSC with estimated states:

	state_seq = filter_states(noi_ts, sources_ts, hist_lookup, L = L)

	state_props = numpy.bincount(state_seq)

	state_probs = state_props / float(numpy.sum(state_props))

	print '\n\n'

	for state in states_counts:
		for hist in state[0]:
			print hist

		print state[1]

		print '\n'

	print '\n\n'

	C = 0

	for prob in state_probs:
		C += prob*numpy.log2(prob)

	C = -C

	# Perform the Miller-Maddow entropy correction,
	# 	H_{MM} = H_{MLE} + (|X| - 1)/(2n)

	Cs[t_ind] = C + (len(state_probs) - 1)/(2*len(state_seq))

	print 'The local statistical complexity is {}...'.format(C)

	# Compute LSC with the correct model.

	state_seq = sources_ts[:, t, :].flatten()

	state_props = numpy.bincount(state_seq)

	state_probs = state_props / float(numpy.sum(state_props))

	C = 0

	for prob in state_probs:
		C += prob*numpy.log2(prob)

	C = -C

	# Perform the Miller-Maddow entropy correction,
	# 	H_{MM} = H_{MLE} + (|X| - 1)/(2n)

	Cs_correct_model[t_ind] = C + (len(state_probs) - 1)/(2*len(state_seq))

pylab.figure()
pylab.plot(ts, Cs)
pylab.show()

with open('Cs-twitter_MM.dat', 'w') as wfile:
	for C in Cs:
		wfile.write('{}\n'.format(C))

pylab.figure()
pylab.plot(ts, Cs_correct_model)
pylab.show()

with open('Cs-twitter_MM_true.dat', 'w') as wfile:
	for C in Cs_correct_model:
		wfile.write('{}\n'.format(C))