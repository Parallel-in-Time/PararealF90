import os
for np in range(1,9):
    cc = 'mpirun -n ' + str(np)  + ' scaling/bin/advection_weakscaling_mpi.out'
    print cc
    os.system(cc)