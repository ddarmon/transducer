# NOTE: stf stands for 'spatio-temporal filter'

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

	# The length of the timeseries.

	n = noi_ts.shape[1]

	# A dictionary of histories. It is of the form
	#   {history : [# of 0s, # of 1s, ..., # of num_symbols]}

	hist_dict = {}

	# Convention: we'll record a history as a tuple of tuples, i.e.
	#	((noi_ts[t-1], noi_ts[t - 2], ..., noi_ts[t - L]), (sources_ts[0][t-1], ..., sources_ts[0][t-L]), ...)
	# This is the simplest way I can think of for doing this.

	for start in range(0, n - L):
		future_symbol = noi_ts[0, start+L]

		cur_hist = []

		cur_hist.append(tuple(noi_ts[0, start:(start+L)]))

		for source_ind in range(sources_ts.shape[2]):
			cur_hist.append(tuple(sources_ts[0, start:(start+L), source_ind]))
		
		# Turn cur_hist from a list to a tuple.

		cur_hist = tuple(cur_hist)

		if cur_hist in hist_dict: # We have found an observed sequence already in our collection
			hist_dict[cur_hist][future_symbol] += 1
		else: # We have not found the observed sequence in our collection yet
			hist_dict[cur_hist] = [0 for ind in range(num_symbols)]
			hist_dict[cur_hist][future_symbol] += 1

	hists = []

	for histories in hist_dict:
		hists.append([[histories], hist_dict[histories]])

	return hist_dict, hists

def generate_hist_dict_ensemble(noi_ts, sources_ts, num_symbols):
	# generate_hist_dict_ensemble takes in an *ensemble* of time
	# series and computes the history counts across
	# ensembles instead of across time. The output is a 
	# dictionary of the form:
	# 	{history : [# of 0s, # of 1s, ..., # of num_symbols]}
	# that is used for the state merging reconstruction of
	# the local causal states.

	# Recall that the full sources_ts is
	# 	N x T x Np
	# where N is the number of trials, T is the number of timesteps
	# and Np is the number of parents of the noi.

	# We'll assume that for a fixed time t, a chunk of this data
	# structure has been passed as
	#
	# 	sources_ts[:, (t-L):t, :]
	#
	# That is, we've passed from time t - L to time t - 1,
	# or L timesteps overall.

	# Similarly, recall that noi_ts is
	# 	N x T x 1

	# We'll assume that for a fixed time t, a chunk of this data
	# structure has been passed as
	#
	# 	noi_ts[:, (t-L):(t+1), :]
	#
	# That is, we pass L timesteps back, along with the future time
	# step to predict.

	# A dictionary of histories. It is of the form
	#   {history : [# of 0s, # of 1s, ..., # of num_symbols]}

	hist_dict = {}

	# Convention: we'll record a history as a tuple of tuples, i.e.
	#	((noi_ts[t-1], noi_ts[t - 2], ..., noi_ts[t - L]), (sources_ts[0][t-1], ..., sources_ts[0][t-L]), ...)
	# This is the simplest way I can think of for doing this.

	N = noi_ts.shape[0] # The number of trials

	for n in range(N):
		future_symbol = noi_ts[n, -1]

		cur_hist = []

		cur_hist.append(tuple(noi_ts[n, 0:-1]))

		for source_ind in range(sources_ts.shape[2]):
			cur_hist.append(tuple(sources_ts[n, :, source_ind]))
		
		# Turn cur_hist from a list to a tuple.

		cur_hist = tuple(cur_hist)

		if cur_hist in hist_dict: # We have found an observed sequence already in our collection
			hist_dict[cur_hist][future_symbol] += 1
		else: # We have not found the observed sequence in our collection yet
			hist_dict[cur_hist] = [0 for ind in range(num_symbols)]
			hist_dict[cur_hist][future_symbol] += 1

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

def csmr(hists, num_symbols, alpha = 0.001, H_test = 'chisquared'):
	# Take in a dictionary hists of the form
	# 	{history : [# of 0s, # of 1s, ..., # of num_symbols]}
	# and output the partition of those histories into
	# inferred causal states, as well as the inferred
	# weighted predictive distribution. This comes in
	# two forms: counts and probs.

	# This algorithm for merging is due to Shalizi. The
	# initialism csmr stands for Causal State Merging 
	# Reconstruction, in homage to CSSR.

	# First, make a copy of hists, since it is passed by
	# reference.

	d_hists = copy.deepcopy(hists)

	df = num_symbols - 1

	numpy.random.shuffle(d_hists) # NOTE: numpy.random.shuffle shuffles in place

	# Create a list to contain the causal states.

	states = []

	nhists = len(d_hists) # the number of histories

	# Move the first past into its own state

	states.append(d_hists.pop(0))

	# Run through each of the remaining histories and see if it belongs to
	# one of the existing states.

	if H_test == 'exact':
		print 'Using Fisher\'s exact test...'
	elif H_test == 'chisquared':
		print 'Using chi-squared test...'

	for historyind, history in enumerate(d_hists):
		nomatch = True # nomatch keeps track, for each history, whether or not we've found
					   # a matching causal state.
		
		for stateind, state in enumerate(states):
			if H_test == 'exact':
				test = exact_test(history, state, alpha = alpha)
			elif H_test == 'chisquared':
				test = chisquared_test(history, state, df = df, alpha = alpha)
			
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

	n = noi_ts.shape[1]

	for start in range(0, n - L):
		cur_hist = []

		cur_hist.append(tuple(noi_ts[0, start:(start+L)]))

		for source_ind in range(sources_ts.shape[2]):
			cur_hist.append(tuple(sources_ts[0, start:(start+L), source_ind]))
		
		# Turn cur_hist from a list to a tuple.

		cur_hist = tuple(cur_hist)
		
		state_seq.append(hist_lookup[cur_hist])
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