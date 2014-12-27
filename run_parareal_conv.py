import os, sys
sys.path.append('./scripts')
from build_namelist import build_namelist
from generate_q0 import generate_q0
from build_runscript import build_runscript
from generate_timemesh import generate_timemesh
nu = 0.005
Nx = 32
Ny = 33
Nz = 34
dt_fine   = 1.0/200.0
dt_coarse = 1.0/40.0
Niter = 1
Tend  = 2.0
do_io = True
be_verbose = False
#
generate_q0(Nx, Ny, Nz)
Nproc = 2

# read the run command to use plus possible options
with open("system.defs", "r") as rfile:
    system = rfile.readline()
    system = system.rstrip()
    runcmd = rfile.readline()
    runcmd = runcmd.rstrip()
    rfile.close()

# Run serial reference
timemesh = generate_timemesh(0.0, Tend, dt_fine, dt_coarse, Nproc)
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
  types = [ 'mpi', 'openmp', 'openmp_pipe' ]
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
          build_runscript(np, jobname, type, system)
          os.system("sbatch submit_"+type+"_Np"+str(Nproc)+".sh")
#  os.system('mv timings*.dat data/')
