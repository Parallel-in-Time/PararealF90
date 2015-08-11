#!/usr/bin/python
import numpy
from matplotlib import pyplot as plt
from pylab import rcParams

def extract_energy(line):
  str_to_find = "energy_used', "
  ind1        = line.find(str_to_find)
  ind1        = ind1 + len(str_to_find)
  ind2        = line.find("]")
  energy      = int(line[ind1:ind2])
  return energy
#
# Main part of script
#

fs = 8
Nsamples = 50

#Nprocs  = numpy.array([2, 4, 6, 8, 12, 24])
Nprocs     = numpy.array([24])
energy     = numpy.zeros([3, Nprocs.size, Nsamples])
energy_avg = numpy.zeros([3, Nprocs.size])
stand_dev  = numpy.zeros([3, Nprocs.size])
confidence = numpy.zeros([3, Nprocs.size])

# Load serial energy
#filename = "9_serial_f_Np1.rur"
#f = open(filename,'r')
#line1 = f.readline()
# Extract energy entry in RUR file
#line2 = f.readline()
#energy_fine = extract_energy(line2)
#f.close

# Load coarse runtime
#filename = "9_serial_g_Np1.rur"
#f = open(filename, 'r')
#line1 = f.readline()
#line2 = f.readline()
#energy_coarse = extract_energy(line2)
#f.close

types  = [ 'mpi', 'openmp', 'openmp_pipe' ]
for tt in range(0,3):
  type = types.pop(0)
  for ii in range(0,Nprocs.size):
    np     = Nprocs[ii]
    for jj in range(0,Nsamples):
      filename = str(jj)+"_"+type+"_Np"+str(np)+".rur"
      f = open(filename,'r')
      line1 = f.readline()
      line2 = f.readline()
      energy[tt,ii,jj] = extract_energy(line2)
      f.close()
    # Compute averages
    for jj in range(0,Nsamples):
      energy_avg[tt,ii] = energy_avg[tt,ii] + energy[tt,ii,jj]/float(Nsamples)
    # Compute standard deviation
    for jj in range(0,Nsamples):
      stand_dev[tt,ii] = stand_dev[tt,ii] + (energy[tt,ii,jj] - energy_avg[tt,ii])**2/float(Nsamples)
    stand_dev[tt,ii] = numpy.sqrt(stand_dev[tt,ii])
    print ("relative standard deviation: %9.3f" % (stand_dev[tt,ii]/energy_avg[tt,ii]) )

    confidence[tt,ii] = 1.96*stand_dev[tt,ii]/numpy.sqrt(float(Nsamples))
    #print ("95 percent confidence: +/- %9.3f" % (confidence[tt,ii]) )
  
print ('Average energy-to-solution for MPI: %7.2f +/- %5.2f joule' % (energy_avg[0,0], confidence[0,0]))  
print ('Average energy-to-solution for OpenMP: %7.2f +/- %5.2f joule' % (energy_avg[1,0], confidence[1,0]))  
print ('Average energy-to-solution for OpenMP-pipe: %7.2f +/- %5.2f joule' % (energy_avg[2,0], confidence[2,0]))  
rcParams['figure.figsize'] = 6, 2.5   
fig, ax = plt.subplots()
ind     = numpy.arange(3)
width   = 0.5
rects_1 = ax.bar(0.5-0.5*width, energy_avg[0,:]/1000, width, color='b', yerr = confidence[0,0]/1000, error_kw=dict(ecolor='k', lw=1, capsize=8, capthick=1.0))
rects_2 = ax.bar(1.5-0.5*width, energy_avg[1,:]/1000, width, color='g', yerr = confidence[1,0]/1000, error_kw=dict(ecolor='k', lw=1, capsize=8, capthick=1.0))
rects_3 = ax.bar(2.5-0.5*width, energy_avg[2,:]/1000, width, color='r', yerr = confidence[2,0]/1000, error_kw=dict(ecolor='k', lw=1, capsize=8, capthick=1.0))
ax.set_xlim([0, 3])
ax.grid()
#ax.set_yscale('log')
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.set_xticks([0.5, 1.5, 2.5])
ax.set_xticklabels( ('MPI', 'OpenMP', 'OpenMP(pipe)') , fontsize=fs)
ax.set_ylabel('Energy-to-solution (kilojoule)', fontsize=fs, labelpad=15)
fig.savefig('Energy.pdf', bbox_inches='tight')
plt.show()