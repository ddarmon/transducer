import numpy

seq = numpy.random.rand(10000)

seq = numpy.ceil(seq - 0.5)

ofile = open('seq.dat', 'w')

for val in seq:
    ofile.write(str(int(val)))

ofile.close()