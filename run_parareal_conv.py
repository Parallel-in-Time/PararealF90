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
Nproc = 8

# read the run command to use plus possible options
with open("system.defs", "r") as rfile:
    system = rfile.readline()
    system = system.rstrip()
    runcmd = rfile.readline()
    runcmd = runcmd.rstrip()
    rfile.close()

# Run serial reference
timemesh = generate_timemesh(0.0, Tend, dt_fine, dt_coarse, 1)
Nfine = timemesh.get('Nfine')
Ncoarse = timemesh.get('Ncoarse')
build_namelist(nu, Nx, Ny, Nz, Nfine, Ncoarse, Niter, Tend, do_io, be_verbose)
if system=="mac":
  os.system("time bin/run_timestepper.out F")
else:
  jobname="fine_serial"
  build_runscript(1, jobname, "serial_f", system)
  os.system("sbatch submit_serial_f_Np1.sh")

#
for Niter in range(1,Nproc):
  #types = [ 'mpi', 'openmp', 'openmp_pipe' ]
  types = ['mpi']
  timemesh = generate_timemesh(0.0, Tend, dt_fine, dt_coarse, Nproc)
  Nfine = timemesh.get('Nfine')
  Ncoarse = timemesh.get('Ncoarse')
  build_namelist(nu, Nx, Ny, Nz, Nfine, Ncoarse, Niter, Tend, do_io, be_verbose)
  for ii in range(0,len(types)):
      type=types.pop(0)
      if system=="mac":
          if type=="mpi":
              os.system("time mpirun -n "+str(Nproc)+" bin/run_parareal_"+type+".out")
          elif type=="openmp":
              os.system("time OMP_NUM_THREADS="+str(Nproc)+" mpirun -n 1 bin/run_parareal_"+type+".out")
          elif type=="openmp_pipe":
              os.system("time OMP_NUM_THREADS="+str(Nproc)+" mpirun -n 1 bin/run_parareal_"+type+".out")
      else:    
          jobname="parareal_"+type+"_Np"+str(Nproc)
          build_runscript(Nproc, jobname, type, system)
          os.system("sbatch submit_"+type+"_Np"+str(Nproc)+".sh")
#  os.system('mv timings*.dat data/')
