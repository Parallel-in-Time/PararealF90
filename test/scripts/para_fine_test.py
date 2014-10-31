##
# Verifies that all three version of Parareal reproduce the
# fine integrator if the number of iterations is equal to the
# number of time slices.
#
import os, sys, numpy
from build_namelist import build_namelist
from generate_q0 import generate_q0

def para_fine_test():
  # define parameter
  nu         = 0.005
  Nx         = 16
  Ny         = 22
  Nz         = 19
  N_fine     = 200
  N_coarse   = 100
  Tend       = 0.2
  do_io      = True
  be_verbose = False
  Np         = 2
  Niter      = 2
  run_cmd    = 'mpirun'
  Np_s       = '%0.2i' % (Np-1)
  Np_s_p1    = '%0.2i' % (Np)

  # Prepare
  os.system('rm -f *.dat')
  generate_q0(Nx, Ny, Nz)
  build_namelist(nu, Nx, Ny, Nz, Np*N_fine, Np*N_coarse, Niter, Tend, do_io, be_verbose)

  # Compute fine reference
  os.system('./bin/run_timestepper.out F')

  # Build namelist for Parareal
  build_namelist(nu, Nx, Ny, Nz,    N_fine,    N_coarse,    Np, Tend, do_io, be_verbose)

  # Parareal-MPI
  os.system(run_cmd+' -n '+str(Np)+' ./bin/run_parareal_mpi.out')
  fser = open('qend.dat','r')
  fpar = open('q_final_'+Np_s+'_'+Np_s_p1+'_mpi.dat')
  max_err = 0.0
  for i in range(0,Nx):
    for j in range(0,Ny):
      for k in range(0,Nz):
        Qser = float(fser.readline())
        Qpar = float(fpar.readline())
        max_err = max( abs(Qser-Qpar), max_err )
  if max_err>1e-14:
    print max_err
    sys.exit("ERROR: Parareal-MPI with Nit=Nproc and fine integrator did not yield identical results.")
  elif numpy.isnan(max_err):
    sys.exit("ERROR: Parareal-MPI with Nit=Nproc and fine integrator produced NaN error.")

  # Parareal-OpenMP
  os.system('OMP_NUM_THREADS='+str(Np)+' '+run_cmd+' -n '+str(Np)+' ./bin/parareal_openmp.out')
  fser = open('qend.dat','r')
  fpar = open('q_final_'+Np_s+'_'+Np_s_p1+'_openmp.dat')
  max_err = 0.0
  for i in range(0,Nx):
    for j in range(0,Ny):
      for k in range(0,Nz):
        Qser = float(fser.readline())
        Qpar = float(fpar.readline())
        max_err = max( abs(Qser-Qpar), max_err )
  if max_err>1e-14:
    print max_err
    sys.exit("ERROR: Parareal-OpenMP with Nit=Nproc and fine integrator did not yield identical results.")
  elif numpy.isnan(max_err):
    sys.exit("ERROR: Parareal-OpenMP with Nit=Nproc and fine integrator produced NaN error.")

  # Parareal-OpenMP-pipe
  os.system('OMP_NUM_THREADS='+str(Np)+' '+run_cmd+' -n '+str(Np)+' ./bin/parareal_openmp_pipe.out')
  fser = open('qend.dat','r')
  fpar = open('q_final_'+Np_s+'_'+Np_s_p1+'_openmp_pipe.dat')
  max_err = 0.0
  for i in range(0,Nx):
    for j in range(0,Ny):
      for k in range(0,Nz):
        Qser = float(fser.readline())
        Qpar = float(fpar.readline())
        max_err = max( abs(Qser-Qpar), max_err )
  if max_err>1e-14:
    print max_err
    sys.exit("ERROR: Parareal-OpenMP-pipe with Nit=Nproc and fine integrator did not yield identical results.")
  elif numpy.isnan(max_err):
    sys.exit("ERROR: Parareal-OpenMP-pipe with Nit=Nproc and fine integrator produced NaN error.")

  print " [0] -- Successful: All three versions of Parareal (MPI, OpenMP, OpenMP-pipe) reproduce the fine integrator for Nit = Np"