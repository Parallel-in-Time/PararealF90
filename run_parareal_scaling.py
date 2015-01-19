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
Nproc = [2, 4, 6, 8]

# read the run command to use plus possible options
with open("system.defs", "r") as rfile:
    system = rfile.readline()
    system = system.rstrip()
    runcmd = rfile.readline()
    runcmd = runcmd.rstrip()
    rfile.close()

timemesh = generate_timemesh(0.0, Tend, dt_fine, dt_coarse, 1)
Nfine = timemesh.get('Nfine')
Ncoarse = timemesh.get('Ncoarse')
param_file = "param_fine_ref.in"
build_namelist(nu, Nx, Ny, Nz, Nfine, Ncoarse, Niter, Tend, do_io, be_verbose, param_file)

# Run serial reference
if system=="mac":
  os.system("time bin/run_timestepper.out "+param_file+" F")
else:
  jobname="fine_serial"
  build_runscript(1, jobname, "serial_f", system, param_file)
  os.system("sbatch submit_serial_f_Np1.sh")

# To compute coarse-to-fine ratio, run also coarse propagator
if system=="mac":
  os.system("time bin/run_timestepper.out "+param_file+" C")
else:
  jobname="coarse_serial"
  build_runscript(1, jobname, "coarse_f", system, param_file)
  os.system("sbatch submit_seria_g_Np1.sh")

for np in Nproc:
  types = [ 'mpi', 'openmp', 'openmp_pipe' ]
  #types = ['mpi']
  timemesh = generate_timemesh(0.0, Tend, dt_fine, dt_coarse, np)
  Nfine = timemesh.get('Nfine')
  Ncoarse = timemesh.get('Ncoarse')
  for ii in range(0,len(types)):
      type=types.pop(0)
      param_file = "param_para_"+type+"_Np"+str(np)+".in"
      build_namelist(nu, Nx, Ny, Nz, Nfine, Ncoarse, Niter, Tend, do_io, be_verbose, param_file)
      if system=="mac":
          if type=="mpi":
              os.system("time mpirun -n "+str(np)+" bin/run_parareal_"+type+".out "+param_file)
          elif type=="openmp":
              os.system("time OMP_NUM_THREADS="+str(np)+" mpirun -n 1 bin/run_parareal_"+type+".out "+param_file)
          elif type=="openmp_pipe":
              os.system("time OMP_NUM_THREADS="+str(np)+" mpirun -n 1 bin/run_parareal_"+type+".out "+param_file)
      else:    
          jobname="parareal_"+type+"_Np"+str(np)
          build_runscript(np, jobname, type, system, param_file)
          os.system("sbatch submit_"+type+"_Np"+str(np)+".sh")
#  os.system('mv timings*.dat data/')
