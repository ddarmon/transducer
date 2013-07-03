import scipy.stats
import numpy
import copy
from collections import deque
import ipdb

def generate_hist_dict(noi_ts, sources_ts, num_symbols, L = 1):
	# generate_hist_dict takes in the adjacent
	# time series and the time series for a 
	# node of interest and outputs a dictionary
	# of the form 
	# 	{history : [# of 0s, # of 1s, ..., # of num_symbols]}
	# that is used for the state merging reconstruction of
	# the local causal states.

	n = len(noi_ts)

	# A dictionary of histories. It is of the form
	#   {history : [# of 0s, # of 1s, ..., # of num_symbols]}

	hist_dict = {}

	for start in range(0, n - L):
	    cur_string = ''

	    for source_ind in range(len(sources_ts)):
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

	return hist_dict, hists

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

def chisquared_test(state1, state2, df, alpha = 0.001):
    # Input: Two lists of the form [[hist1, hist2, ...], [# 0s, #1s, ..., #num_symbols]]
    #
    # Output: A boolean with whether we reject the null hypothesis (and do not merge),
    # or do not reject the null hypothesis (and do merge)
    
    chisquared_statistic = 0

    S1 = numpy.sum(state1[1])
    S2 = numpy.sum(state2[1])

    K1 = numpy.sqrt(S2 / float(S1))
    K2 = 1/K1

    for ind in range(df + 1):
        numerator = K1*state1[1][ind] - K2*state2[1][ind]
    
        denominator = state1[1][ind] + state2[1][ind]
    
        chisquared_statistic += numerator**2/float(denominator)

	# The quantile (inverse CDF) for a standard chi-square variable

	quantile = scipy.stats.chi.ppf
    
    test = (chisquared_statistic > quantile(1-alpha, df)) # If false, we reject the null hypothesis. Otherwise, we do not.
    
    return test

def exact_test(state1, state2, alpha = 0.001):
    # Perform Fisher's exact test.

    print 'Warning: The Exact test currently assumes a *two symbol* alphabet.'
    
    c00, c01 = state1[1]
    c10, c11 = state2[1]

    oddsratio, pvalue = scipy.stats.fisher_exact([[c00, c01], [c10, c11]])

    test = (pvalue < alpha)

    return test

def csmr(hists, num_symbols, alpha = 0.001):
	# Take in a dictionary hists of the form
	# 	{history : [# of 0s, # of 1s, ..., # of num_symbols]}
	# and output the partition of those histories into
	# inferred causal states, as well as the inferred
	# weighted predictive distribution. This comes in
	# two forms: counts and probs.

	# This algorithm for merging is due to Shalizi. The
	# initialism csmr stands for Causal State Merging 
	# Reconstruction, in homage to CSSR.

	df = num_symbols - 1

	numpy.random.shuffle(hists) # NOTE: numpy.random.shuffle shuffles in place

	# Create a list to contain the causal states.

	states = []

	nhists = len(hists) # the number of histories

	# Move the first past into its own state

	states.append(hists.pop(0))

	# Run through each of the remaining histories and see if it belongs to
	# one of the existing states.

	# print 'Using Fisher\'s exact test...'
	print 'Using chi-squared test...'

	for historyind, history in enumerate(hists):
	    nomatch = True # nomatch keeps track, for each history, whether or not we've found
	                   # a matching causal state.
	    
	    for stateind, state in enumerate(states):
	        test = chisquared_test(history, state, df = df, alpha = alpha)
	        # test = exact_test(history, state, alpha = alpha)
	        
	        if test == False: 
	            states[stateind] = merge_states(history, state)
	            
	            nomatch = False # We've found a match, so we indicate that in nomatch
	            
	            break # We've merged, so break out of the inner for loop
	        
	    if nomatch == True:
	        states.append(history)

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

	return states, states_final, hist_lookup

def filter_states(noi_ts, sources_ts, hist_lookup, L = 1):
	# Returns the local causal state sequence for the node
	# using the estimated partition returned by csmr.

	state_seq = []

	n = len(noi_ts)

	for start in range(0, n - L):
	    cur_string = ''

	    for source_ind in range(len(sources_ts)):
	        cur_string += sources_ts[source_ind][start:(start+L)]

	    cur_string += noi_ts[start:(start+L)]
	    
	    state_seq.append(hist_lookup[cur_string])
	return state_seq

def predict(noi_ts, sources_ts, hist_lookup, states_probs, L = 1):
	# Predict the next symbol in the sequence using the 
	# predictive distributions inferred from csmr.

	n = len(noi_ts)

	prediction = ''

	for start in range(0, n - L):
	    cur_string = ''

	    for source_ind in range(len(sources_ts)):
	        cur_string += sources_ts[source_ind][start:(start+L)]

	    cur_string += noi_ts[start:(start+L)]
	    
	    prediction += str(numpy.argmax(states_probs[hist_lookup[cur_string]][1]))

	return prediction