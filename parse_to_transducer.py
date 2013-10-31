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

# numpy.random.seed(1) # Fix the random number generator so we get
					   # reproducible results.

adj_file = 'edge_list_3K_user_connected_directed.txt'
# adj_file = 'adj_mat_toy.txt'

# datatype = 'timeseries_synthetic/toy_transducer'
datatype = 'timeseries_synthetic/twitter_p1_i2_ensemble'

# for noi_ind in range(5):
for noi_ind in ['0']:
	# Node of interest.

	noi = str(noi_ind)

	print 'On node {}...'.format(noi)

	# Read in the parents of the node of interest,
	# which act as *sources/inputs* to the node of
	# interest.

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

	if len(sources) == 0:
		pass
	else:
		# sources_ts contains all of the time series
		# for the inputs *into* a particular node (in
		# the case of a transducer) or *adjacent to*
		# a particular node (in the case of a spatio-
		# temporal random field).

		sample_count = 0

		with open('{}/sample{}.dat'.format(datatype, sources[0])) as ofile:
			for line in ofile:
				sample_count += 1

			if ';' in line: # If we have an alphabet of size greater than 10.
				symbol_count = len(line.strip()[:-1].split(';'))
			else: # If we have an alphabet of size less than or equal to 10.
				symbol_count = len(line.strip())

		sources_ts = numpy.zeros((sample_count, symbol_count, len(sources)), dtype = 'int8')

		print sources

		for source_index, source in enumerate(sources):
			with open('{}/sample{}.dat'.format(datatype, source)) as ofile:
				for sample_index, line in enumerate(ofile):
					if ';' in line:
						ts = line.strip()[:-1].split(';')

						sources_ts[sample_index, :, source_index] = cur_ts[sample_index, :] = map(int, ts)
					else:
						ts = line.strip()

						sources_ts[sample_index, :, source_index] = numpy.fromstring(ts, dtype = 'int8') - 48

		# A particular compression, namely the
		# number of neighbors tweeting at the previous
		# timestep.

		with open('{}/input{}.dat'.format(datatype, noi), 'w') as wfile:
			for sample_ind in range(sample_count):
				for t in range(symbol_count):
					cur_tot = numpy.sum(sources_ts[sample_ind, t, :]) # Sum up the number of inputs active in the current sample / timestep

					wfile.write('{};'.format(cur_tot))

				wfile.write('\n')