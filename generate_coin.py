import numpy

seq = numpy.random.rand(100000)

seq = numpy.ceil(seq - 0.5)

ofile = open('coin.dat', 'w')

for val in seq:
    ofile.write(str(int(val)))

ofile.close()