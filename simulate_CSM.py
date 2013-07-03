import cssr_interface
import pylab
import numpy
import collections
import ipdb

class state:
	def __init__(self, p_emit0 = None, s_emit0 = None, s_emit1 = None):
		self.p_emit0 = p_emit0
		self.s_emit0 = s_emit0
		self.s_emit1 = s_emit1
	
	
	def setEmit0State(self, state, prob):
		self.s_emit0 = state
		self.p_emit0 = prob
	
	def setEmit1State(self, state):
		self.s_emit1 = state

# The users with the highest statistical complexity

# suffixes = ['14664756', '17789501', '18773467', '19713521', '18049840'] 

# The users with the greatest CSM Rate - ESN Rate

# CSM_file = '/Users/daviddarmon/Dropbox/papers/socialcom-2013/paper/CSM-Simulations/best-csm.txt'

# suffixes = [line.rstrip('\n') for line in open(CSM_file)]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# User options:
num_its = 20000

to_fix_initial_state = False

num_sims = 1

tosave = True # Whether or not to save the timeseries to a file.

toplot = False # Whether or not to display a raster plot of the sequence

#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname = 'complex-csm.dot'

ofile = open(fname)

line = ofile.readline()

while line[0] != '0':
	line = ofile.readline()

CSM = collections.defaultdict(state)

while line != '}':
	lsplit = line.split()

	from_state = lsplit[0]
	to_state = lsplit[2]
	esymbol = int(lsplit[5][1])
	eprob = float(lsplit[6])
	
	if esymbol == 0:
		CSM[from_state].setEmit0State(to_state, eprob)
	elif esymbol == 1:
		CSM[from_state].setEmit1State(to_state)
	
	line = ofile.readline()

ofile.close()

n_states = len(CSM)

urand = numpy.random.rand(num_its, num_sims)

symbol_seq = numpy.empty((num_its, num_sims), dtype = 'int32')

for cur_sim in xrange(num_sims):
	if to_fix_initial_state:
		cur_state = '0' # Fix the initial state across trials
	else:
		cur_state = str(numpy.random.randint(0, n_states)) # Set the initial state at random


	for ind in range(num_its):
		# The following probability isn't foolproof, since
		# we could run into the problem of pemit0 = None

		if urand[ind, cur_sim] < CSM[cur_state].p_emit0:
			symbol_seq[ind, cur_sim] = 0
			cur_state = CSM[cur_state].s_emit0
		else:
			symbol_seq[ind, cur_sim] = 1
			cur_state = CSM[cur_state].s_emit1

if toplot:
	f, axarr = pylab.subplots(num_sims, sharex=True)

	for axind in range(num_sims):
		axarr[axind].vlines(numpy.arange(symbol_seq.shape[0])[symbol_seq[:, axind] == 1], -0.5, 0.5)
		axarr[axind].yaxis.set_visible(False)

		axarr[axind].tick_params(
	        axis='x',          # changes apply to the x-axis
	        which='both',      # both major and minor ticks are affected
	        bottom='off',      # ticks along the bottom edge are off
	        top='off',         # ticks along the top edge are off
	        color = 'white')

	pylab.show()

if tosave:
	ofile = open('{}.dat'.format(fname.split('.')[0]), 'w')

	for trial in xrange(symbol_seq.shape[1]):
		for symbol in symbol_seq[:, trial]:
			ofile.write(str(symbol))

		ofile.write('\n')

	ofile.close()