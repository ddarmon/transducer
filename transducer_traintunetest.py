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
import glob

from stf_methods import *
from traintunetest import *

import cssr_interface
from filter_data_methods import *

# numpy.random.seed(1) # Fix the random number generator so we get
                       # reproducible results.

# Node of interest.

noi = '1'

L = 1
L_CSSR = L
L_max  = L_CSSR


# sources_ts contains all of the time series
# for the inputs *into* a particular node (in
# the case of a transducer) or *adjacent to*
# a particular node (in the case of a spatio-
# temporal random field).

datatype = 'timeseries_synthetic/toy_transducer'
# datatype = 'timeseries_synthetic/twitter_p1_i2'

source = '{}/input{}'.format(datatype, noi)

try: # Some nodes have *no* input, and we must handle those.
  with open('{}.dat'.format(source)) as ofile:
    line = ofile.readline()

    if ';' in line:
      ts = line.strip()[:-1].split(';')

      sources_ts = numpy.array(map(int, ts))[numpy.newaxis, :, numpy.newaxis]
    else:
      ts = line.strip()

      sources_ts = numpy.array(map(int, ts))[numpy.newaxis, :, numpy.newaxis]
except IOError, e:
  sources_ts = None 

# noi_ts contains the time series for the node that
# we wish to predict.

with open('{}/sample{}.dat'.format(datatype, noi)) as ofile:
  line = ofile.readline()

  if ';' in line:
    ts = line.strip()[:-1].split(';')

    noi_ts = numpy.array(map(int, ts))[numpy.newaxis, :]
  else:
    ts = line.strip()

    noi_ts = numpy.array(map(int, ts))[numpy.newaxis, :]

n = noi_ts.shape[1]

num_symbols = 2 # The number of possible symbols

hist_dict, hists = generate_hist_dict(noi_ts = noi_ts, sources_ts = sources_ts, num_symbols = num_symbols, L = L)

df = num_symbols - 1

# Perform the hypothesis test:
#
#       H0: Two distributions are equal
#       H1: Two distributions are *not* equal
#
# at a level alpha


alpha = 0.05

states_counts, states_probs, hist_lookup = csmr(hists, num_symbols, alpha = alpha)

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

if numpy.mean(noi_ts) < 0.5:
    majority_vote = 0
else:
    majority_vote = 1

for ind in range(len(prediction)):
    true = noi_ts[0, L+ind]
    pred = prediction[ind]

    if true == pred:
       correct += 1

    if true == majority_vote:
       base_rate += 1

    # print true, pred

print 'Accuracy rate: {}...'.format(float(correct)/len(prediction))
print 'Compared to using a biased coin {}...'.format(numpy.max((float(base_rate)/len(prediction), 1 - float(base_rate)/len(prediction))))

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

for state_ind, state in enumerate(states_final):
  print '\nState {}\n========'.format(state_ind)
  enum_hists, enum_probs = state

  for hist in enum_hists:
    print hist

  print enum_probs

print '\n\n\n\n'

for hist in hist_dict:
  count_0, count_1 = hist_dict[hist]

  print hist, count_1/float(count_0 + count_1), count_0, count_1