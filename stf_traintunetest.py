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
from traintunetest import *

# numpy.random.seed(1) # Fix the random number generator so we get
                       # reproducible results.

ofile = open('adj_mat.txt')

ofile.readline()

line = ofile.readline()

# Node of interest.

noi = '1'

sources = []

while line != '':
    source, dest, weight = line.split(',')

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

datatype = 'timeseries_synthetic/Bernoulli'

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
#       H0: Two distributions are equal
#       H1: Two distributions are *not* equal
#
# at a level alpha


alpha = 0.001

states_counts, states_probs, hist_lookup = csmr(hists, num_symbols, alpha = 0.001)

states_final = states_probs

# print states_final

state_seq = filter_states(noi_ts, sources_ts, hist_lookup, L = L)

state_props = numpy.bincount(state_seq)

state_fractions = state_props / float(numpy.sum(state_props))

prediction = predict(noi_ts, sources_ts, hist_lookup, states_probs, L = L)

# num_forward = 1000

# for ind in range(num_forward):
# 	print '{}\t{}'.format(noi_ts[L+ind], prediction[ind])

correct = 0
base_rate = 0

for ind in range(len(prediction)):
	true = noi_ts[L+ind]
	pred = prediction[ind]

	if true == pred:
		correct += 1

	if true == '1':
		base_rate += 1

print 'The accuracy rate is {}...'.format(float(correct)/len(prediction))
print 'Compared to using a biased coin {}...'.format(float(base_rate)/len(prediction))