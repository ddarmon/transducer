import numpy
import sys

def categorical(weights, objects):
    """Return a random item from objects, with the weighting defined by weights 
    (which must sum to 1)."""
    cs = numpy.cumsum(weights) #An array of the weights, cumulatively summed.
    idx = sum(cs < numpy.random.rand()) #Find the index of the first weight over a random value.
    return objects[idx]

num_symbols = 100000

syms = [0, 1]

state_sequence = numpy.empty(num_symbols)
symbol_sequence = numpy.empty(num_symbols-1)

state_sequence[0] = 1

sym_probs = {0 : [0.5, 0.5], 1 : [0, 1]}

trans = {0 : {0 : 0, 1 : 1}, 1 : {1 : 0}}

for ind in xrange(num_symbols-1):
    symbol_sequence[ind] = categorical(sym_probs[state_sequence[ind]], syms)
    
    state_sequence[ind+1] = trans[state_sequence[ind]][symbol_sequence[ind]]

ofile = open('even.dat', 'w')

for sym in symbol_sequence:
    # sys.stdout.write('{0} '.format(int(sym)))
    ofile.write('{0}'.format(int(sym)))

ofile.close()