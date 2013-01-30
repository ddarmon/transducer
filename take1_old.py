import scipy.stats
import numpy

execfile("generate_examples.py")

ofile = open('seq.dat')

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
    hists.append(histories)

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

# def merge_histories():
#     # Input:
#     #
#     # Output:
#     #
    

for hist_ind in range(1, len(hists)):
    chisquared_statistic = 0

    S1 = numpy.sum(hist_dict[hists[0]])
    S2 = numpy.sum(hist_dict[hists[hist_ind]])

    K1 = numpy.sqrt(S2 / float(S1))
    K2 = 1/K1

    for ind in range(num_symbols):
        numerator = K1*hist_dict[hists[0]][ind] - K2*hist_dict[hists[hist_ind]][ind]
    
        denominator = hist_dict[hists[0]][ind] + hist_dict[hists[hist_ind]][ind]
    
        chisquared_statistic += numerator**2/float(denominator)
    
    test = (chisquared_statistic > quantile(1-alpha, df)) # If false, we reject the null hypothesis. Otherwise, we do not.
    
    if test == False:        
        print 'According to the chi-squared test, the two histories merge.'
    elif test == True:
        print 'According to the chi-squared test, the two histories do not merge.'
        
    print hist_dict[hists[0]], hist_dict[hists[hist_ind]]