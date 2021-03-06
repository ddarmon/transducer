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

import cssr_interface
from filter_data_methods import *

# numpy.random.seed(1) # Fix the random number generator so we get
                       # reproducible results.

adj_file = 'edge_list_3K_user_connected_directed.txt'
# adj_file = 'adj_mat_toy.txt'

ofile = open(adj_file)

ofile.readline()

line = ofile.readline()

# Node of interest.

noi = '17'

L = 1 # The past to consider

sources = []

weighted = False

if weighted:
    while line != '':
        if ',' in line:
            source, dest, weight = line.split(',')
        elif '\t' in line:
            source, dest, weight = line.split('\t')
        else:
            source, dest, weight = line.split(' ')

        if dest == noi:
            sources.append(source)

        line = ofile.readline()

    ofile.close()
else:
    while line != '':
        if ',' in line:
            source, dest = line.strip().split(',')
        elif '\t' in line:
            source, dest = line.strip().split('\t')
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
#   print '{}\t{}'.format(noi_ts[L+ind], prediction[ind])

correct = 0
base_rate = 0

if numpy.mean(numpy.fromstring(noi_ts, dtype = 'int8') - 48) < 0.5:
    majority_vote = '0'
else:
    majority_vote = '1'

for ind in range(len(prediction)):
    true = noi_ts[L+ind]
    pred = prediction[ind]

    if true == pred:
       correct += 1

    if true == majority_vote:
       base_rate += 1

    # print true, pred

print 'Accuracy rate: {}...'.format(float(correct)/len(prediction))
print 'Compared to using a biased coin {}...'.format(numpy.max((float(base_rate)/len(prediction), 1 - float(base_rate)/len(prediction))))

L_CSSR = L
L_max  = L_CSSR

metric = 'accuracy'

fname = '{}/sample{}'.format(datatype, noi)

zero_order_predict = generate_zero_order_CSM(fname)

cssr_interface.run_CSSR(filename = fname, L = L_CSSR, savefiles = True, showdot = True, is_multiline = False, showCSSRoutput = False)

CSM = get_CSM(fname = fname)

epsilon_machine = get_epsilon_machine(fname = fname)

states, L = get_equivalence_classes(fname) # A dictionary structure with the ordered pair
                                                             # (symbol sequence, state)

correct_rates = run_tests(fname = fname, CSM = CSM, zero_order_CSM = zero_order_predict, states = states, epsilon_machine = epsilon_machine, L = L_CSSR, L_max = L_max, metric = metric, print_predictions = False, print_state_series = False)

print 'Compared to using CSSR on the timeseries {}...'.format(correct_rates[0])

# print hist_dict

# for hist in hist_dict:
#   count_0, count_1 = hist_dict[hist]

#   print hist, count_1/float(count_0 + count_1), count_0, count_1