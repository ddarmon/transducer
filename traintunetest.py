####################################################################################
#
#	A file that takes in a .dat file (a time series that has been discretized)
#	and outputs 

import sys
import pylab
import numpy
import ipdb

def create_traintunetest_single(fname, ratios = (0.8, 0.1, 0.1), toprint = False, shuffle = False):
	# This function takes in the file name for a
	# data file and outputs three files, one for
	# training, one for tuning, and one for
	# testing. This assumes that the time series
	# file contains a *single* line and is not
	# multi-line.

	ofile = open('{0}.dat'.format(fname))

	# Pull out the time series.

	timeseries = ofile.readline().rstrip('\n')

	ofile.close()

	N = len(timeseries) # The number of time points in the time series.

	assert numpy.sum(ratios) == 1, "Warning: Your train / tune / test ratios should sum to 1."

	# Get the number of training/tuning/testing
	# time points.

	ntrain = int(numpy.ceil(N*ratios[0]))
	ntune = int(numpy.ceil(N*ratios[1]))

	ntest = N - (ntrain + ntune)

	# If we should shuffle the train/tune set, 
	# create a shuffled set of indices.

	if toprint:
		print 'ntrain: {0}\nntune: {1}\nntest: {2}\n'.format(ntrain, ntune, ntest)

	# Generate the train/tune/test datasets.

	overallfile = open('{0}-train+tune.dat'.format(fname), 'w')

	trainfile = open('{0}-train.dat'.format(fname), 'w')

	trainfile.write('{0}'.format(timeseries[:ntrain]))
	overallfile.write('{0}'.format(timeseries[:ntrain]))

	trainfile.close()

	tunefile = open('{0}-tune.dat'.format(fname), 'w')

	tunefile.write('{0}'.format(timeseries[ntrain:(ntrain + ntune)]))
	overallfile.write('{0}'.format(timeseries[ntrain:(ntrain + ntune)]))

	tunefile.close()

	overallfile.close()

	testfile = open('{0}-test.dat'.format(fname), 'w')

	testfile.write('{0}'.format(timeseries[(ntrain + ntune):]))

	testfile.close()