##
# Verifies that all three variants of Parareal (MPI, OpenMP and
# OpenMP-pipelined) yield identical results. In order to rule out
# possible race conditions in the OpenMP versions, the same setup
# is run and compared multiple times (default 5).
#
import os, sys, numpy
from build_namelist import build_namelist
from generate_q0 import generate_q0
import multiprocessing
from termcolor import colored

def para_compare_test(run_cmd):
  #
  nu = 0.005
  Nx = 13
  Ny = 15
  Nz = 14
  N_fine   = 40
  N_coarse = 20
  Niter = 1
  Tend  = 0.2
  do_io = True
  be_verbose = False
  Np       = multiprocessing.cpu_count()
  Ntests   = 10
  Np_s     = '%0.2i' % (Np-1)
  Np_s_p1  = '%0.2i' % Np

  for nn in range(0,Ntests):
    sys.stdout.write('Running comparison test %2i of %2i. \r' % (nn, Ntests) )
    sys.stdout.flush()
    #
    #
    os.system('rm -f *.dat')
    generate_q0(Nx, Ny, Nz)
    build_namelist(nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose)
    os.system('OMP_NUM_THREADS=1 '+run_cmd+' -n '+str(Np)+' ./bin/run_parareal_mpi.out')
    os.system('OMP_NUM_THREADS='+str(Np)+' '+run_cmd+' -n 1 ./bin/parareal_openmp_pipe.out')
    os.system('OMP_NUM_THREADS='+str(Np)+' '+run_cmd+' -n 1 ./bin/parareal_openmp.out')
    
    # Compare MPI to OpenMP-pipe
    Np_s = '%0.2i' % (Np)
    for nt in range(0,Np):
      if nt<10:
        nt_s = '0'+str(nt)
      else:
        nt_s = str(nt)
      fmpi    = open('q_final_'+nt_s+'_'+Np_s+'_mpi.dat')
      fopenmp = open('q_final_'+nt_s+'_'+Np_s+'_openmp_pipe.dat')
      max_err = 0.0
      for i in range(0,Nx):
        for j in range(0,Ny):
          for z in range(0,Nz):
            Qmpi    = float(fmpi.readline())
            Qopenmp = float(fopenmp.readline())
            max_err = max( abs(Qmpi - Qopenmp), max_err )
      if max_err>1e-14:
        print 'Timeslice: '+str(nt)
        print max_err
        sys.exit(colored('ERROR: Parareal-MPI and Parareal-OpenMP-pipe did not yield identical results.','red'))
      elif numpy.isnan(max_err):
        sys.exit(colored('ERROR: Parareal-MPI and Parareal-OpenMP-pipe produced NaN error','red'))

    # Compare MPI to OpenMP
    Np_s = '%0.2i' % (Np)
    for nt in range(0,Np):
      if nt<10:
        nt_s = '0'+str(nt)
      else:
        nt_s = str(nt)
        fmpi    = open('q_final_'+nt_s+'_'+Np_s+'_mpi.dat')
        fopenmp = open('q_final_'+nt_s+'_'+Np_s+'_openmp.dat')
        max_err = 0.0
        for i in range(0,Nx):
          for j in range(0,Ny):
            for z in range(0,Nz):
              Qmpi    = float(fmpi.readline())
              Qopenmp = float(fopenmp.readline())
              max_err = max( abs(Qmpi - Qopenmp), max_err )
      if max_err>1e-14:
        print 'Timeslice: '+str(nt)
        print max_err
        sys.exit(colored('ERROR: Parareal-MPI and Parareal-OpenMP did not yield identical results.','red'))
      elif numpy.isnan(max_err):
        sys.exit(colored('ERROR: Parareal-MPI and Parareal-OpenMP produced NaN error','red'))

  print colored(" [0] -- Successful: All three versions of Parareal (MPI, OpenMP, OpenMP-pipe) give identical results.",'green')
