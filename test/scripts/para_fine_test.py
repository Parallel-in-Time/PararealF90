##
# Verifies that all three version of Parareal reproduce the
# fine integrator if the number of iterations is equal to the
# number of time slices.
#
import os, sys, numpy
from build_namelist import build_namelist
from generate_q0 import generate_q0
from get_run_cmd import get_run_cmd
import random as rnd
import multiprocessing
from termcolor import colored

def para_fine_test(system, run_cmd, max_cpu, be_verbose, Ntests):
  for nn in range(0,Ntests):
    sys.stdout.write('Running fine test %2i of %2i. \r' % (nn, Ntests) )
    sys.stdout.flush()
    
    # define parameter
    nu         = rnd.uniform(0.001, 0.009)
    Nx         = rnd.randint(16, 25)
    Ny         = rnd.randint(16, 25)
    Nz         = rnd.randint(16, 25)
    N_fine     = rnd.randint(40, 60)
    N_coarse   = rnd.randint(25, 30)
    
    Tend       = 0.2
    do_io      = True
    if (max_cpu==1):
      Np = 2
    else:
      Np         = rnd.randint(2,max_cpu)
    
    Niter      = Np
    Np_s       = '%0.2i' % (Np-1)
    Np_s_p1    = '%0.2i' % (Np)

    # Prepare
    generate_q0(Nx, Ny, Nz)
    build_namelist(nu, Nx, Ny, Nz, Np*N_fine, Np*N_coarse, Niter, Tend, do_io, be_verbose)

    # Compute fine reference
    run_cmd_full = get_run_cmd(system, run_cmd, True, 1)
    if be_verbose:
      print run_cmd_full
    os.system(run_cmd_full+' ./bin/run_timestepper.out F')

    # Build namelist for Parareal
    build_namelist(nu, Nx, Ny, Nz,    N_fine,    N_coarse,    Np, Tend, do_io, be_verbose)

    # Parareal-MPI
    run_cmd_full = get_run_cmd(system, run_cmd, True, Np)
    if be_verbose:
      print run_cmd_full
    os.system(run_cmd_full+' ./bin/run_parareal_mpi.out')
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
      sys.exit(colored("ERROR: Parareal-MPI with Nit=Nproc and fine integrator did not yield identical results.",'red'))
    elif numpy.isnan(max_err):
      sys.exit(colored("ERROR: Parareal-MPI with Nit=Nproc and fine integrator produced NaN error.",'red'))

    # Parareal-OpenMP
    run_cmd_full = get_run_cmd(system, run_cmd, False, Np)
    if be_verbose:
      print run_cmd_full
    os.system(run_cmd_full+' ./bin/run_parareal_openmp.out')
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
      sys.exit(colored("ERROR: Parareal-OpenMP with Nit=Nproc and fine integrator did not yield identical results.",'red'))
    elif numpy.isnan(max_err):
      sys.exit(colored("ERROR: Parareal-OpenMP with Nit=Nproc and fine integrator produced NaN error.",'red'))

    # Parareal-OpenMP-pipe
    run_cmd_full = get_run_cmd(system, run_cmd, False, Np)
    if be_verbose:
      print run_cmd_full
    os.system(run_cmd_full+' ./bin/run_parareal_openmp_pipe.out')
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
      sys.exit(colored("ERROR: Parareal-OpenMP-pipe with Nit=Nproc and fine integrator did not yield identical results.",'red'))
    elif numpy.isnan(max_err):
      sys.exit(colored("ERROR: Parareal-OpenMP-pipe with Nit=Nproc and fine integrator produced NaN error.",'red'))

  print colored(" [0] -- Successful: All three versions of Parareal (MPI, OpenMP, OpenMP-pipe) reproduce the fine integrator for Nit = Np", 'green')