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

sources_ts = []

datatype = 'timeseries_synthetic/Bernoulli'

# dataset = 8
# datatype = 'timeseries/NEURO-Set' + str(dataset)

for source in sources:
    ofile = open('{}/sample{}.dat'.format(datatype, source))

    sources_ts.append(ofile.readline().rstrip('\n'))

    ofile.close()

ofile = open('{}/sample{}.dat'.format(datatype, noi))

noi_ts = ofile.readline().rstrip('\n')

ofile.close()

n = len(noi_ts)

L = 1 # The past to consider

num_symbols = 2 # The number of possible symbols

# A dictionary of histories. It is of the form
#   {history : [# of 0s, # of 1s, ..., # of num_symbols]}

hist_dict = {}

for start in range(0, n - L):
    cur_string = ''

    for source_ind in range(len(sources)):
        cur_string += sources_ts[source_ind][start:(start+L)]

    cur_string += noi_ts[start:(start+L+1)]
    
    if cur_string[0:-1] in hist_dict: # We have found an observed sequence already in our collection
        hist_dict[cur_string[0:-1]][int(cur_string[-1])] += 1
    else: # We have not found the observed sequence in our collection yet
        hist_dict[cur_string[0:-1]] = [0 for ind in range(num_symbols)]
        hist_dict[cur_string[0:-1]][int(cur_string[-1])] += 1

hists = []

for histories in hist_dict:
    hists.append([[histories], hist_dict[histories]])

# The quantile (inverse CDF) for a standard chi-square variable

quantile = scipy.stats.chi.ppf

df = num_symbols - 1

# Perform the hypothesis test:
#
#       H0: Two distributions are equal
#       H1: Two distributions are *not* equal
#
# at a level alpha

alpha = 0.001

def merge_states(state1, state2):
    # Input: Two lists of the form [[hist1, hist2, ...], [# 0s, #1s, ..., #num_symbols]]
    #
    # Output: One list of the same form, with the histories combined and the symbol counts updated
    
    merged_state = copy.copy(state1)
    
    for hist in state2[0]:
        merged_state[0].append(hist)
    
    for ind in range(len(state2[1])):
        merged_state[1][ind] += state2[1][ind]
    
    return merged_state

def chisquared_test(state1, state2):
    # Input: Two lists of the form [[hist1, hist2, ...], [# 0s, #1s, ..., #num_symbols]]
    #
    # Output: A boolean with whether we reject the null hypothesis (and do not merge),
    # or do not reject the null hypothesis (and do merge)
    
    chisquared_statistic = 0

    S1 = numpy.sum(state1[1])
    S2 = numpy.sum(state2[1])

    K1 = numpy.sqrt(S2 / float(S1))
    K2 = 1/K1

    for ind in range(num_symbols):
        numerator = K1*state1[1][ind] - K2*state2[1][ind]
    
        denominator = state1[1][ind] + state2[1][ind]
    
        chisquared_statistic += numerator**2/float(denominator)
    
    test = (chisquared_statistic > quantile(1-alpha, df)) # If false, we reject the null hypothesis. Otherwise, we do not.
    
    return test

def exact_test(state1, state2):
    # Perform Fisher's exact test.
    
    c00, c01 = state1[1]
    c10, c11 = state2[1]

    oddsratio, pvalue = scipy.stats.fisher_exact([[c00, c01], [c10, c11]])

    test = (pvalue < alpha)

    return test

# This algorithm for merging is due to Shalizi.

numpy.random.shuffle(hists) # NOTE: numpy.random.shuffle shuffles in place

# Create a list to contain the causal states.

states = []

nhists = len(hists) # the number of histories

# Move the first past into its own state

states.append(hists.pop(0))

# Run through each of the remaining histories and see if it belongs to
# one of the existing states.

print 'Using Fisher\'s exact test...'

for historyind, history in enumerate(hists):
    nomatch = True # nomatch keeps track, for each history, whether or not we've found
                   # a matching causal state.
    
    for stateind, state in enumerate(states):
        # test = chisquared_test(history, state)
        test = exact_test(history, state)
        
        if test == False: 
            states[stateind] = merge_states(history, state)
            
            nomatch = False # We've found a match, so we indicate that in nomatch
            
            break # We've merged, so break out of the inner for loop
        
    if nomatch == True:
        states.append(history)
    
# print states

states_final = []
hist_lookup  = {}

for state_label, state in enumerate(states):
    histories, props = state

    tot = numpy.sum(props)

    probs = []

    for prop in props:
        probs.append(prop/float(tot))

    states_final.append([histories, probs])

    for history in histories:
        hist_lookup[history] = state_label

print states_final

state_seq = []

for start in range(0, n - L):
    cur_string = ''

    for source_ind in range(len(sources)):
        cur_string += sources_ts[source_ind][start:(start+L)]

    cur_string += noi_ts[start:(start+L)]
    
    state_seq.append(hist_lookup[cur_string])

state_props = numpy.bincount(state_seq)

state_probs = state_props / float(numpy.sum(state_props))

C = 0

for prob in state_probs:
    C += prob*numpy.log2(prob)

C = -C

print 'The local statistical complexity is {}...'.format(C)