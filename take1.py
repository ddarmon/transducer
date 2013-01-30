import scipy.stats
import numpy
import copy
from collections import deque

execfile("even_process.py")

ofile = open('symseq.dat')

seq = ofile.readline()

n = len(seq) # The number of time points

ofile.close()

L = 2 # The past to consider

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

alpha = 0.001

# NOTE: I'm not convinced this is working. At least,
# it doesn't seem to be working as well as I'd like it
# to be working. Which may be more of a statement about me
# than the hypothesis test.

def merge_states(state1, state2):
    # Input: Two lists of the form [[hist1, hist2, ...], [# 0s, #1s, ..., #num_symbols]]
    #
    # Output: One list of the same form, with the histories combined and the symbol counts updated
    
    merged_state = state1
    
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

# for hist_ind in range(1, len(hists)):
#     test = chisquared_test(hists[0], hists[hist_ind])
#     
#     if test == False:        
#         print 'According to the chi-squared test, the two histories merge.'
#     elif test == True:
#         print 'According to the chi-squared test, the two histories do not merge.'
#         
#     print hists[0][1], hists[hist_ind][1]

for runthrough in range(10):
    num_states = len(hists)
    i = 0
    j = 1

    while i < num_states-1:
        while j < num_states:
            test = chisquared_test(hists[i], hists[j])
        
            if test == False:
                new_state = merge_states(hists[i], hists[j])
                
                # Remove the two states from hists
                
                deque((list.pop(hists, k) for k in sorted([i, j], reverse=True)), maxlen=0)
                
                # Append the merged state to hists
                
                hists.append(new_state)
            
                j = j - 1
            
                num_states = num_states - 1
            elif test == True:
                print 'Do not merge on {}, {}'.format(hists[i][1], hists[j][1])
        
            j = j + 1
        i = i + 1

print hists