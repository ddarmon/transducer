# This file takes in a multivariate timeseries and
# generates a new timeseries in the
#
# 	input -> transducer -> output
#
# format.
#
# For example, for network data, the the timeseries
# might be all of the behavior inputs to some node,
# and this is `filtered' down to the *number* of
# nodes which are active.

# 	DMD, 041013-13-22

import numpy
import copy
from collections import deque
import ipdb
import sys

from stf_methods import *

print '\n\n\n\n\n\nWARNING!!! This does not currently handle input symbols with more than 0..9 ... \n\n\n\n\n'

# numpy.random.seed(1) # Fix the random number generator so we get
					   # reproducible results.

adj_file = 'edge_list_3K_user_connected_directed.txt'
# adj_file = 'adj_mat_toy.txt'

# datatype = 'timeseries_synthetic/toy_transducer'
datatype = 'timeseries_synthetic/twitter_p1_i2'

# for noi_ind in range(5):
for noi_ind in ['17']:
	# Node of interest.

	noi = str(noi_ind)

	print 'On node {}...'.format(noi)

	with open(adj_file) as ofile:
		weighted = False

		ofile.readline()

		line = ofile.readline()

		sources = []

		if weighted:
			while line != '':
				if '\t' in line:
					source, dest, weight = line.strip().split('\t')
				elif ',' in line:
					source, dest, weight = line.strip().split(',')
				else:
					source, dest, weight = line.strip().split(' ')

				if dest == noi:
					sources.append(source)

				line = ofile.readline()
		else:
			while line != '':
				if '\t' in line:
					source, dest = line.strip().split('\t')
				elif ',' in line:
					source, dest = line.strip().split(',')
				else:
					source, dest = line.strip().split(' ')

				if dest == noi:
					sources.append(source)

				line = ofile.readline()

	# sources_ts contains all of the time series
	# for the inputs *into* a particular node (in
	# the case of a transducer) or *adjacent to*
	# a particular node (in the case of a spatio-
	# temporal random field).

	sources_ts = []

	print sources

	for source in sources:
		ofile = open('{}/sample{}.dat'.format(datatype, source))

		sources_ts.append(ofile.readline().rstrip('\n'))

		ofile.close()

	if len(sources_ts) == 0:
		pass
	else:
		T = len(sources_ts[0])

		# A particular compression, namely the
		# number of neighbors tweeting at the previous
		# timestep.

		with open('{}/input{}.dat'.format(datatype, noi), 'w') as wfile:
			for t in xrange(T):
				cur_tot = 0
				for source in sources_ts:
					cur_tot += int(source[t])

				if cur_tot > 9:
					print 'Sorry, the code doesn\'t work with that kind of alphabet yet.'

					sys.exit(1)

				wfile.write('{}'.format(cur_tot))