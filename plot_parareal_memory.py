#!/usr/bin/python
import numpy
from matplotlib import pyplot as plt
from pylab import rcParams
from matplotlib.patches import Ellipse, Polygon
from subprocess import call

fs = 8
compiler = 'cray'

# Two functions to be used below
def extract_memory(line):
  str_to_find = "max_rss', "
  ind1 = line.find(str_to_find)
  ind1 = ind1 + len(str_to_find)
  ind2 = line.find("rchar")
  # To ignore the , ' in the RUR file before rchar, ignore the last three entries
  memory = int(line[ind1:ind2-3])
  return memory
  
Nprocs  = numpy.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24])

memory  = numpy.zeros([3,Nprocs.size])

# Load serial energy
filename = "serial_f_Np1.rur"
f = open(filename,'r')

# Extract memory entry in RUR file
line1 = f.readline()
memory_fine = extract_memory(line1)/1024
f.close

# Load coarse runtime
filename = "serial_g_Np1.rur"
f = open(filename, 'r')
line1 = f.readline()
memory_coarse = extract_memory(line1)/1024
f.close

types  = [ 'mpi', 'openmp', 'openmp_pipe' ]
for tt in range(0,3):
  type = types.pop(0)
  for ii in range(0,Nprocs.size):
    np     = Nprocs[ii]
    filename = type+"_Np"+str(np)+".rur"
    f = open(filename,'r')
    line1 = f.readline()
    if type=='mpi':
      memory[tt,ii] = np*extract_memory(line1)/1024
    else:
      memory[tt,ii] = extract_memory(line1)/1024
    line2 = f.readline()
    f.close()
    
types  = [ 'mpi', 'openmp', 'openmp_pipe' ]
print ("Memory fine serial propagator (MB): %4i" % memory_fine)
print ("Memory fine coarse propagator (MB): %4i" % memory_coarse)
for tt in range(0,3):
  type = types.pop(0)
  for ii in range(0,Nprocs.size):
    print ("Memory for Parareal "+type+" on "+str(ii)+" cores (MB): %4i" % memory[tt,ii])
    
rcParams['figure.figsize'] = 2.5, 2.5   
fig, ax = plt.subplots()
ind = numpy.arange(Nprocs.size)
width = 0.4
rects1 = ax.bar( ind,         memory[0,:], width, color='b', hatch='x')
rects2 = ax.bar( ind+width,   memory[2,:], width, color='r', hatch='\\')
#rects3 = ax.bar( ind+2*width, memory[2,:], width, color='r', hatch='-')

ax.plot( ind+1.0*width, memory_fine*Nprocs, 'k-', linewidth=1.0)
ax.legend( (rects1[0], rects2[0]), ('MPI','OpenMP','OpenMP(pipe)'), loc=2, fontsize=fs)
ax.set_xticks(ind+1.0*width)
ax.set_xticklabels( ('2', '4', '6', '8', '10', '12', '14', '16', '18', '20', '22', '24'))
for label in ax.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=fs)
ax.set_xlabel(r'Number of cores $P$', fontsize=fs)
ax.set_ylabel(r'Memory in MByte $m(P)$', fontsize=fs, labelpad=5)
ax.set_ylim([0, 350])
#plt.show()
filename = 'memory_dora_'+compiler+'.pdf'
fig.savefig(filename, bbox_inches='tight')
call(["pdfcrop", filename, filename])
