import numpy
from matplotlib import pyplot as plt

fs = 18

Nprocs = numpy.array([2, 4, 6, 8, 12, 24])
Niter  = 2
timers  = numpy.zeros([3,Nprocs.size])
speedup = numpy.zeros([3,Nprocs.size])
bound   = numpy.zeros([1,Nprocs.size])

# Load serial runtime first
filename = "timings_serial_fine.dat"
f = open(filename,'r')
time_serial_f = float(f.readline())
f.close

# Load coarse runtime
filename = "timings_serial_coarse.dat"
f = open(filename, 'r')
time_serial_g = float(f.readline())
f.close
g_to_f = time_serial_g/time_serial_f

niter_s = '%0.2i' % Niter

types  = [ 'mpi', 'openmp', 'openmp_pipe' ]
for tt in range(0,3):
  type = types.pop(0)
  for ii in range(0,Nprocs.size):
    np     = Nprocs[ii]
    nps    = '%0.2i' % np
    nps_m1 = '%0.2i' % (np-1)
    filename = "timings_"+type+niter_s+"_"+nps_m1+"_"+nps+".dat"
    f = open(filename,'r')
    total_time  = float(f.readline())
    fine_time   = float(f.readline())
    coarse_time = float(f.readline())
    comm_time   = float(f.readline())
    timers[tt,ii] = fine_time + coarse_time + comm_time
    f.close()
    speedup[tt,ii] = time_serial_f/timers[tt,ii]
    bound[0,ii]    = 1.0/( (1.0 + float(Niter)/float(np))*g_to_f + float(Niter)/float(np))

fig = plt.figure(figsize=(8,8))

plt.plot(Nprocs, speedup[0,:], linewidth=0, marker='^', markersize=fs, color='b', label='MPI')
plt.plot(Nprocs, speedup[1,:], linewidth=0, marker='<', markersize=fs, color='g', label='OpenMP')
plt.plot(Nprocs, speedup[2,:], linewidth=0, marker='>', markersize=fs, color='r', label='OpenMP(pipe)')
plt.plot(Nprocs, bound[0,:], linewidth=1.0, color='k', label='Bound')
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

fig = plt.figure(figsize=(8,8))

plt.plot(Nprocs, timers[0,:], linewidth=1.0, marker='^', markersize=fs, color='b', label='MPI')
plt.plot(Nprocs, timers[1,:], linewidth=1.0, marker='<', markersize=fs, color='g', label='OpenMP')
plt.plot(Nprocs, timers[2,:], linewidth=.0, marker='>', markersize=fs, color='r', label='OpenMP(pipe)')
plt.plot(Nprocs, time_serial_f + 0.0*timers[0,:], linewidth=2.0, color='k')
nodes = list(Nprocs)
#ymin = 0
#ymax = max(map(max,speedup))+1.0

plt.xlabel('Number of cores', fontsize=fs)
plt.ylabel('Runtime [sec.]', fontsize=fs)
plt.xticks(Nprocs, fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)
plt.legend(loc='upper right')
#plt.ylim([ymin, ymax])

# Saveing figure
fig.savefig('Runtime.pdf', dpi=300)

plt.show()