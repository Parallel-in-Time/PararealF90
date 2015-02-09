#!/usr/bin/python
import numpy
from matplotlib import pyplot as plt

# Two functions to be used below
def extract_memory(line):
  str_to_find = "max_rss', "
  ind1 = line.find(str_to_find)
  ind1 = ind1 + len(str_to_find)
  ind2 = line.find("rchar")
  # To ignore the , ' in the RUR file before rchar, ignore the last three entries
  memory = int(line[ind1:ind2-3])
  return memory

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

fs = 18

Nprocs  = numpy.array([2, 4, 6, 8, 12, 24])
energy  = numpy.zeros([3,Nprocs.size])
memory  = numpy.zeros([3,Nprocs.size])

# Load serial energy
filename = "serial_f_Np1.rur"
f = open(filename,'r')

# Extract memory entry in RUR file
line1 = f.readline()
memory_fine = extract_memory(line1)

# Extract energy entry in RUR file
line2 = f.readline()
energy_fine = extract_energy(line2)
f.close

# Load coarse runtime
filename = "serial_g_Np1.rur"
f = open(filename, 'r')
line1 = f.readline()
memory_coarse = extract_memory(line1)
line2 = f.readline()
energy_coarse = extract_energy(line2)
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
      memory[tt,ii] = np*extract_memory(line1)
    else:
      memory[tt,ii] = extract_memory(line1)
    line2 = f.readline()
    energy[tt,ii] = extract_energy(line2)
    f.close()

print memory_fine
print memory_coarse
print memory

print ""
print energy_fine
print energy_coarse
print energy