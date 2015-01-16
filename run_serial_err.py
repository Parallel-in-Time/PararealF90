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
dx    = 1.0/(max(Nx, Ny, Nz))
print ('Approximate coarse advective CFL number: %6.2f' % (dt_coarse/dx) )
print ('Approximate coarse diffusive CFL number: %6.2f' % (nu*dt_coarse/(dx*dx)))
print ('Approximate fine advective CFL number:   %6.2f' % (dt_fine/dx) )
print ('Approximate fine diffusive CFL number:   %6.2f' % (nu*dt_fine/(dx*dx)))

# read the run command to use plus possible options
with open("system.defs", "r") as rfile:
    system = rfile.readline()
    system = system.rstrip()
    runcmd = rfile.readline()
    runcmd = runcmd.rstrip()
    rfile.close()

# Run serial reference with extra fine time step
timemesh = generate_timemesh(0.0, Tend, dt_fine, dt_coarse, 1)
Nfine    = timemesh.get('Nfine')
Nfine    = 10*Nfine
Ncoarse  = timemesh.get('Ncoarse')
build_namelist(nu, Nx, Ny, Nz, Nfine, Ncoarse, Niter, Tend, do_io, be_verbose)
if system=="mac":
  os.system("time bin/run_timestepper.out F")
  # rename output file
  os.system('mv q_final_fine.dat q_final_fine_ref.dat')
else:
  jobname="fine_serial"
  build_runscript(1, jobname, "serial_f_ref", system)
  # append line to rename output file
  with open('submit_serial_f_ref_Np1.sh','a') as fileobj:
    fileobj.write("mv q_final_fine.dat q_final_fine_ref.dat \n")
    fileobj.close()
  os.system("sbatch submit_serial_f_ref_Np1.sh")

# Run serial fine scheme with normal time step
timemesh = generate_timemesh(0.0, Tend, dt_fine, dt_coarse, 1)
Nfine    = timemesh.get('Nfine')
Ncoarse  = timemesh.get('Ncoarse')
build_namelist(nu, Nx, Ny, Nz, Nfine, Ncoarse, Niter, Tend, do_io, be_verbose)
if system=="mac":
  os.system("time bin/run_timestepper.out F")
else:
  jobname="fine_serial"
  build_runscript(1, jobname, "serial_f", system)
  os.system("sbatch submit_serial_f_Np1.sh")

# Run serial coarse scheme
timemesh = generate_timemesh(0.0, Tend, dt_fine, dt_coarse, 1)
Nfine    = timemesh.get('Nfine')
Ncoarse  = timemesh.get('Ncoarse')
build_namelist(nu, Nx, Ny, Nz, Nfine, Ncoarse, Niter, Tend, do_io, be_verbose)
if system=="mac":
  os.system("time bin/run_timestepper.out C")
else:
  jobname="coarse_serial"
  build_runscript(1, jobname, "serial_g", system)
  os.system("sbatch submit_serial_g_Np1.sh")