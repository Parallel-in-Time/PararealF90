import numpy
Nprocs = 8
total  = numpy.zeros((3,Nprocs))
fine   = numpy.zeros((3,Nprocs))
coarse = numpy.zeros((3,Nprocs))
comm   = numpy.zeros((3,Nprocs))
types  = [ 'mpi', 'openmp', 'openmp_pipe' ]
#
for tt in range(0,3):
    type = types.pop(0)
    for nn in range(0,Nprocs):
        if nn<10:
            filename = "timings_"+type+"_0"+str(nn)+".dat"
        else:
            filename = "timings_"+type+"_0"+str(nn)+".dat"
        #        
        f = open(filename, 'r')
        total[tt,nn]  = f.readline()
        fine[tt,nn]   = f.readline()
        coarse[tt,nn] = f.readline()
        comm[tt,nn]   = f.readline()
time_sum = numpy.zeros(3)
time_total = numpy.zeros(3)
quotient   = numpy.zeros(3)
for tt in range(0,3):
    sum            = fine[tt,:] + coarse[tt,:] + comm[tt,:]
    temp           = total[tt,:]
    time_total[tt] = numpy.amax(temp)
    time_sum[tt]   = numpy.amax(sum)
    quotient[tt]   = time_sum[tt]/time_total[tt]
print time_total
print time_sum
print quotient
