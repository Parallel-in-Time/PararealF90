import numpy
import matplotlib as plt

Nprocs = numpy.array([2, 4, 8])
print Nprocs
timers  = numpy.zeros([3,Nprocs.size])
speedup = numpy.zeros([3,Nprocs.size])
# Load serial runtime first
filename = "timings_serial_fine.dat"
f = open(filename,'r')
time_serial_f = float(f.readline())
f.close

types  = [ 'mpi', 'openmp', 'openmp_pipe' ]
for tt in range(0,3):
  type = types.pop(0)
  for ii in range(0,Nprocs.size):
    np     = Nprocs[ii]
    nps    = '%0.2i' % np
    nps_m1 = '%0.2i' % (np-1)
    filename = "timings_"+type+nps_m1+"_"+nps+".dat"
    f = open(filename,'r')
    timers[tt,ii] = float(f.readline())
    f.close()
    speedup[tt,ii] = time_serial_f/timers[tt,ii]


