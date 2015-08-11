import numpy
from matplotlib import pyplot as plt
from pylab import rcParams

fs = 8

Nsamples   = 4

machine = "cub"

if machine=="dora":
  Nprocs     = numpy.array([2, 4, 6, 8, 10, 12, 24])
if machine=="cub":
  Nprocs     = numpy.array([2, 4, 6, 8])
Niter      = 4
timers     = numpy.zeros([3, Nprocs.size, Nsamples])
timers_avg = numpy.zeros([3, Nprocs.size])
speedup    = numpy.zeros([3, Nprocs.size])
bound      = numpy.zeros([1, Nprocs.size])
stand_dev  = numpy.zeros([3, Nprocs.size])
confidence = numpy.zeros([3,numpy.size(Nprocs)])

# Load serial runtime first
if Nsamples==1:
  filename = "timings_serial_fine.dat"
else:
  filename = "0_timings_serial_fine.dat"

f = open(filename,'r')
time_serial_f = float(f.readline())
f.close

# Load coarse runtime
if Nsamples==1:
  filename = "timings_serial_coarse.dat"
else:
  filename = "0_timings_serial_coarse.dat" 
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
    for jj in range(0,Nsamples):
      if Nsamples==1:
        filename = "timings_"+type+niter_s+"_"+nps_m1+"_"+nps+".dat"
      else:
        filename = str(jj)+"_timings_"+type+niter_s+"_"+nps_m1+"_"+nps+".dat"
      f = open(filename,'r')
      total_time  = float(f.readline())
      fine_time   = float(f.readline())
      coarse_time = float(f.readline())
      comm_time   = float(f.readline())
      timers[tt,ii,jj] = fine_time + coarse_time + comm_time
      f.close()
    # Compute averages
    for jj in range(0,Nsamples):
      timers_avg[tt,ii] = timers_avg[tt,ii] + timers[tt,ii,jj]/float(Nsamples)
    # Compute standard deviation
    stand_dev[tt,ii] = 0.0
    for jj in range(0,Nsamples):
        stand_dev[tt,ii] = stand_dev[tt,ii] + (timers[tt,ii,jj] - timers_avg[tt,ii])**2/float(Nsamples)
    #    
    stand_dev[tt,ii] = numpy.sqrt(stand_dev[tt,ii]) 
    print ("relative standard deviation: %9.3f" % (stand_dev[tt,ii]/timers_avg[tt,ii]) )
    #  
    confidence[tt,ii] = 1.96*stand_dev[tt,ii]/numpy.sqrt(float(Nsamples))
    
    speedup[tt,ii] = time_serial_f/timers_avg[tt,ii]
    bound[0,ii]    = 1.0/( (1.0 + float(Niter)/float(np))*g_to_f + float(Niter)/float(np))

rcParams['figure.figsize'] = 2.5, 2.5    

fig = plt.figure()
plt.plot(Nprocs, speedup[0,:], linewidth=1.0, marker='^', markersize=fs, color='b', label='MPI')
plt.plot(Nprocs, speedup[1,:], linewidth=1.0, marker='<', markersize=fs, color='g', label='OpenMP')
plt.plot(Nprocs, speedup[2,:], linewidth=1.0, marker='>', markersize=fs, color='r', label='OpenMP(pipe)')
plt.plot(Nprocs,   bound[0,:], linewidth=1.0, color='k', label='Bound')
nodes = list(Nprocs)
ymin  = 0
ymax  = max(map(max,speedup))+1.0

plt.xlabel('Number of cores', fontsize=fs)
plt.ylabel('Speedup', fontsize=fs, labelpad=2)
if machine=="dora":
  plt.xticks([2,6,10,14,18,22], fontsize=fs)
if machine=="cub":
  plt.xticks(Nprocs, fontsize=fs)

plt.yticks(fontsize=fs)
plt.grid(True)
plt.legend(loc='upper left', fontsize=fs, prop={'size':fs})
plt.ylim([ymin, ymax])

# Saveing figure
fig.savefig('Speedup.pdf',bbox_inches='tight')

plt.show()

fig = plt.figure()
plt.plot(Nprocs, timers_avg[0,:], linewidth=1.0, marker='^', markersize=fs, color='b', label='MPI')
plt.plot(Nprocs, timers_avg[1,:], linewidth=1.0, marker='<', markersize=fs, color='g', label='OpenMP')
plt.plot(Nprocs, timers_avg[2,:], linewidth=1.0, marker='>', markersize=fs, color='r', label='OpenMP(pipe)')
plt.plot(Nprocs, time_serial_f + 0.0*timers_avg[0,:], linewidth=1.0, color='k')
nodes = list(Nprocs)
ymax = max(map(max,timers_avg))+1.0
 
NN = Nprocs[ numpy.size(Nprocs) - 2] - 0.75
  
plt.gca().annotate('Serial runtime', xy=( NN, 1.075*time_serial_f), xytext=( NN, 1.075*time_serial_f ), fontsize=fs-1)
plt.gca().set_yscale('log')
#plt.gca().set_xscale('log')
plt.xlabel('Number of cores', fontsize=fs)
plt.ylabel('Runtime [sec.] (log-scaled)', fontsize=fs, labelpad=2)

plt.tick_params(axis='both', which='major', labelsize=fs)

plt.gca().set_ylim([2.0, 50.0])
if machine=="dora":
  plt.gca().set_xticks([2,6,10,14,18,22])
if machine=="cub":
  plt.gca().set_xticks(Nprocs)
plt.gca().set_yticks([5, 10, 20, 50])
plt.gca().set_yticklabels(["5", "10", "20", "50"])
plt.gca().get_yaxis().get_major_formatter().labelOnlyBase = False
plt.grid(True)
if machine=="dora":
  plt.legend(loc='upper right', prop={'size':fs})
if machine=="cub":
  plt.legend(loc='lower left', prop={'size':fs})

# Saveing figure
fig.savefig('Runtime.pdf',bbox_inches='tight')

plt.show()

print "For MPI ..."
for i in range(0,numpy.size(Nprocs)):
  print "Average runtime for #proc = %2i : %5.3f +/- %5.3f" % (Nprocs[i], timers_avg[0,i], confidence[0,i])
print ""
print "For OpenMP ..."
for i in range(0,numpy.size(Nprocs)):
  print "Average runtime for #proc = %2i : %5.3f +/- %5.3f" % (Nprocs[i], timers_avg[1,i], confidence[1,i])
print ""
print "For OpenMP-pipe ..."
for i in range(0,numpy.size(Nprocs)):
  print "Average runtime for #proc = %2i : %5.3f +/- %5.3f" % (Nprocs[i], timers_avg[2,i], confidence[2,i])
