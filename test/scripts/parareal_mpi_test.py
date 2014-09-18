import os, sys
sys.path.append('../../scripts')
from build_namelist import build_namelist
from generate_q0 import generate_q0
nu = 0.005
Nx = 32
Ny = 32
Nz = 32
N_fine = 10 
N_coarse = 1
Niter = 0
Tend  = 1.0
do_io = True
build_namelist(nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io)
generate_q0(Nx, Ny, Nz)
Nproc = 8
os.system('cp ../../bin/parareal_mpi.out .')
os.system("mpirun -n "+str(Nproc)+" ./parareal_mpi.out")