import scipy.stats
import numpy
import copy
from collections import deque
import ipdb

# numpy.random.seed(1) # Fix the random number generator so we get
                     # reproducible results.

# execfile("even_process.py")
ofile = open('even.dat')

# execfile("generate_coin.py")
# ofile = open('coin.dat')

seq = ofile.readline()

n = len(seq) # The number of time points

ofile.close()

L = 3 # The past to consider

num_symbols = 2 # The number of possible symbols

# A dictionary of histories. It is of the form
#   {history : [# of 0s, # of 1s, ..., # of num_symbols]}

hist_dict = {}

for start in range(0, n - L):
    cur_string = seq[start:(start+L+1)]
    
    if cur_string[0:L] in hist_dict: # We have found an observed sequence already in our collection
        hist_dict[cur_string[0:L]][int(cur_string[L])] += 1
    else: # We have not found the observed sequence in our collection yet
        hist_dict[cur_string[0:L]] = [0 for ind in range(num_symbols)]
        hist_dict[cur_string[0:L]][int(cur_string[L])] += 1

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
# at a level alpha

alpha = 0.05

# NOTE: I'm not convinced this is working. At least,
# it doesn't seem to be working as well as I'd like it
# to be working. Which may be more of a statement about me
# than the hypothesis test.

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

# This algorithm for merging is due to Shalizi.

numpy.random.shuffle(hists) # NOTE: numpy.random.shuffle shuffles in place

# Create a list to contain the causal states.

states = []

nhists = len(hists) # the number of histories

# Move the first past into its own state

states.append(hists.pop(0))

# Run through each of the remaining histories and see if it belongs to
# one of the existing states.

for historyind, history in enumerate(hists):
    nomatch = True # nomatch keeps track, for each history, whether or not we've found
                   # a matching causal state.
    
    for stateind, state in enumerate(states):
        test = chisquared_test(history, state)
        
        if test == False: 
            states[stateind] = merge_states(history, state)
            
            nomatch = False # We've found a match, so we indicate that in nomatch
            
            break # We've merged, so break out of the inner for loop
        
    if nomatch == True:
        states.append(history)
    
print states