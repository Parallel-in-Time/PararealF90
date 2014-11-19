import numpy
Nprocs = 16
total  = numpy.zeros((3,Nprocs))
fine   = numpy.zeros((3,Nprocs))
coarse = numpy.zeros((3,Nprocs))
comm   = numpy.zeros((3,Nprocs))
types  = [ 'mpi', 'openmp', 'openmp_pipe' ]
#
for tt in range(0,3):
    type = types.pop(0)
    for nn in range(0,Nprocs):
        nn_str = '%0.2i' % nn
        Np_str = '%0.2i' % Nprocs
        filename = "timings_"+type+nn_str+"_"+Np_str+".dat"
        #
        f = open(filename, 'r')
        total[tt,nn]  = f.readline()
        fine[tt,nn]   = f.readline()
        coarse[tt,nn] = f.readline()
        comm[tt,nn]   = f.readline()
time_fine=numpy.zeros(3)
time_coarse = numpy.zeros(3)
time_comm = numpy.zeros(3)
time_sum = numpy.zeros(3)
time_total = numpy.zeros(3)
coverage   = numpy.zeros(3)

for tt in range(0,3):
  time_fine[tt]   = fine[tt,Nprocs-1]
  time_coarse[tt] = coarse[tt,Nprocs-1]
  time_comm[tt]   = comm[tt,Nprocs-1]
  time_sum[tt]    = time_fine[tt] + time_coarse[tt] + time_comm[tt]
  time_total[tt]  = total[tt,Nprocs-1]
  coverage[tt]    = time_sum[tt]/time_total[tt]
  #    time_fine[tt]   = numpy.amax(fine[tt,:])
  #  time_coarse[tt] = numpy.amax(coarse[tt,:])
  #  time_comm[tt]   = numpy.amax(comm[tt,:])
  #  sum            = fine[tt,:] + coarse[tt,:] + comm[tt,:]
  #  temp           = total[tt,:]
  #  time_total[tt] = numpy.amax(temp)
  #  time_sum[tt]   = numpy.amax(sum)
  #  coverage[tt]   = time_sum[tt]/time_total[tt]

print ("Timings:         MPI  --  OpenMp -- OpenMp_pipe")
print ("Total time:   %8.5f - %8.5f - %8.5f" % (time_total[0], time_total[1], time_total[2]) )
print ("Fine time:    %8.5f - %8.5f - %8.5f" % (time_fine[0], time_fine[1], time_fine[2]) )
print ("Coarse time:  %8.5f - %8.5f - %8.5f" % (time_coarse[0], time_coarse[1], time_coarse[2]) )
print ("Comm time:    %8.5f - %8.5f - %8.5f" % (time_comm[0], time_comm[1], time_comm[2]) )
print ("Summed time:  %8.5f - %8.5f - %8.5f" % (time_sum[0], time_sum[1], time_sum[2]) )
print ("Coverage:     %8.5f - %8.5f - %8.5f" % (coverage[0], coverage[1], coverage[2]) )
