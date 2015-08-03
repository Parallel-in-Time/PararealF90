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
Nproc = 24
Nsamples = 40

# Prepare generation of a bash script to rename RUR output...
# These log files are identified using the JobID. Running rename.sh after all jobs are finished will
# rename these from something like rur.123456 into e.g. serial_f_Np1.rur
os.system("rm -f rename.sh")
os.system("touch rename.sh")
os.system("chmod u=rwx rename.sh")

# read the run command to use plus possible options
with open("system.defs", "r") as rfile:
    system = rfile.readline()
    system = system.rstrip()
    runcmd = rfile.readline()
    runcmd = runcmd.rstrip()
    rfile.close()

timemesh = generate_timemesh(0.0, Tend, dt_fine, dt_coarse, Nproc)
Nfine = timemesh.get('Nfine')
Ncoarse = timemesh.get('Ncoarse')
  
for kk in range(14,15):
  types = [ 'mpi', 'openmp', 'openmp_pipe' ]
  for ii in range(0,len(types)):
      
      type=types.pop(0)
      param_file = "param_para_"+type+"_Np"+str(Nproc)+".in"
      build_namelist(nu, Nx, Ny, Nz, Nfine, Ncoarse, Niter, Tend, do_io, be_verbose, param_file)
      
      if system=="mac":
          print "Can run energy measurements on Cray only"
      else:    
          jobname=str(kk)+"_"+"parareal_"+type+"_Np"+str(Nproc)
          build_runscript(Nproc, jobname, type, system, param_file, kk)
          os.system("sbatch submit_"+str(kk)+"_"+type+"_Np"+str(Nproc)+".sh")
print ""
print "=========================================================================================="
print "REMEMBER: Execute script rename.sh after all jobs completed to rename RUR output files ..."
print "=========================================================================================="
