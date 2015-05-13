import os, sys
sys.path.append('./scripts')
from build_namelist import build_namelist
from generate_q0 import generate_q0
from build_runscript import build_runscript
from generate_timemesh import generate_timemesh
from get_parameter import get_parameter

nu, Nx, Ny, Nz, dt_fine, dt_coarse, Niter, Tend, do_io, be_verbose = get_parameter()
#
generate_q0(Nx, Ny, Nz)
#Nproc = [2, 4, 8]
Nproc = [8]

# read the run command to use plus possible options
with open("system.defs", "r") as rfile:
    system = rfile.readline()
    system = system.rstrip()
    runcmd = rfile.readline()
    runcmd = runcmd.rstrip()
    rfile.close()
for np in Nproc:
  #types = [ 'mpi', 'openmp', 'openmp_pipe' ]
  types = [ 'mpi' ]
  timemesh = generate_timemesh(0.0, Tend, dt_fine, dt_coarse, np)
  Nfine = timemesh.get('Nfine')
  Ncoarse = timemesh.get('Ncoarse')
  for ii in range(0,len(types)):
      type=types.pop(0)
      param_file = "param_para_"+type+"_Np"+str(np)+".in"
      build_namelist(nu, Nx, Ny, Nz, Nfine, Ncoarse, Niter, Tend, do_io, be_verbose, param_file)
      if system=="mac":
          if type=="mpi":
              os.system("time "+runcmd+" -n "+str(np)+" bin/run_parareal_"+type+".out "+param_file)
          elif type=="openmp":
              os.system("time OMP_NUM_THREADS="+str(np)+" mpirun -n 1 bin/run_parareal_"+type+".out "+param_file)
          elif type=="openmp_pipe":
              os.system("time OMP_NUM_THREADS="+str(np)+" mpirun -n 1 bin/run_parareal_"+type+".out "+param_file)
      else:    
          jobname="parareal_"+type+"_Np"+str(np)
          build_runscript(np, jobname, type, system, param_file)
          os.system("sbatch submit_"+type+"_Np"+str(np)+".sh")
  
