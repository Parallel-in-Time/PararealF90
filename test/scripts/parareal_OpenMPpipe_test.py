import os, sys, numpy
sys.path.append('../../scripts')
from build_namelist import build_namelist
from generate_q0 import generate_q0
# Use to switch on or off some tests

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
Np    = 8
Ntests_f = 5
Ntests_c = 5
Ntests   = 5
run_cmd  = 'mpirun'
Np_s = '%0.2i' % (Np-1)
Np_s_p1 = '%0.2i' % Np

os.system('cp ../../bin/run_timestepper.out .')
os.system('cp ../../bin/run_parareal_mpi.out .')
os.system('cp ../../bin/parareal_openmp.out .')
os.system('cp ../../bin/parareal_openmp_pipe.out .')

for nn in range(0,Ntests_c):
    sys.stdout.write('Running test %2i of %2i. \r' % (nn, Ntests_c) )
    sys.stdout.flush()
    #
    os.system('rm -f *.dat')
    # MPI version
    generate_q0(Nx, Ny, Nz)
    build_namelist(nu, Nx, Ny, Nz, Np*N_fine, Np*N_coarse, Niter, Tend, do_io, be_verbose)
    os.system('./run_timestepper.out C')
    # Run Parareal with Niter=0
    build_namelist(nu, Nx, Ny, Nz,    N_fine,    N_coarse,     0, Tend, do_io, be_verbose)
    os.system('OMP_NUM_THREADS='+str(Np)+' ./parareal_openmp_pipe.out')
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
        sys.exit("ERROR: Parareal-OpenMP-Pipe with Nit=0 and coarse integrator did not yield identical results.")
    elif numpy.isnan(max_err):
        sys.exit("ERROR: Parareal-OpenMP-Pipe with Nit=0 and coarse integrator produced NaN error.")

for nn in range(0,Ntests_f):
    sys.stdout.write('Running test %2i of %2i. \r' % (nn, Ntests_f) )
    sys.stdout.flush()
    #
    #
    os.system('rm -f *.dat')
    # OpenMP-pipe version
    generate_q0(Nx, Ny, Nz)
    build_namelist(nu, Nx, Ny, Nz, Np*N_fine, Np*N_coarse, Niter, Tend, do_io, be_verbose)
    os.system('./run_timestepper.out F')    
    # run parareal with Niter=Np
    build_namelist(nu, Nx, Ny, Nz,    N_fine,    N_coarse,    Np, Tend, do_io, be_verbose)
    os.system('OMP_NUM_THREADS='+str(Np)+' ./parareal_openmp_pipe.out')
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
        print "\n"
        print max_err
        sys.exit("ERROR: Parareal-OpenMP-pipe with Nit=Nproc and fine integrator did not yield identical results.")
    elif numpy.isnan(max_err):
        sys.exit("ERROR: Parareal-OpenMP-pipe with Nit=Nproc and fine integrator produced NaN error.")

for nn in range(0,Ntests):
    sys.stdout.write('Running test %2i of %2i. \r' % (nn, Ntests) )
    sys.stdout.flush()
    #
    #
    os.system('rm -f *.dat')
    generate_q0(Nx, Ny, Nz)
    build_namelist(nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose)
    os.system(run_cmd+' -n '+str(Np)+' ./run_parareal_mpi.out')
    os.system('OMP_NUM_THREADS='+str(Np)+' ./parareal_openmp_pipe.out')
    # add OpenMP-pipe
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
            sys.exit('ERROR: Parareal-MPI and Parareal-OpenMP-pipe did not yield identical results.')
        elif numpy.isnan(max_err):
            sys.ext('ERROR: Parareal-MPI and Parareal-OpenMP-pipe produced NaN error')
   
print "[0] -- Test of Parareal-OpenMP-pipe successful"
