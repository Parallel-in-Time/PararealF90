import numpy
from matplotlib import pyplot as plt

fs = 18

Nprocs = numpy.array([2, 4, 6, 8])
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
    total_time  = float(f.readline())
    fine_time   = float(f.readline())
    coarse_time = float(f.readline())
    comm_time   = float(f.readline())
    timers[tt,ii] = fine_time + coarse_time + comm_time
    f.close()
    speedup[tt,ii] = time_serial_f/timers[tt,ii]

fig = plt.figure(figsize=(8,8))

plt.plot(Nprocs, speedup[0,:], linewidth=0, marker='^', markersize=fs, color='b', label='MPI')
plt.plot(Nprocs, speedup[1,:], linewidth=0, marker='<', markersize=fs, color='g', label='OpenMP')
plt.plot(Nprocs, speedup[2,:], linewidth=0, marker='>', markersize=fs, color='r', label='OpenMP(pipe)')

nodes = list(Nprocs)
ymin = 0
ymax = max(map(max,speedup))+1.0

plt.xlabel('Number of cores', fontsize=fs)
plt.ylabel('Speedup', fontsize=fs)
plt.xticks(Nprocs, fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)
plt.legend(loc='upper left')
plt.ylim([ymin, ymax])

# Saveing figure
fig.savefig('Speedup.pdf', dpi=300)

plt.show()