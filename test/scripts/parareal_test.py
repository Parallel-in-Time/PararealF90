import os, sys, numpy
sys.path.append('../../scripts')
from build_namelist import build_namelist
from generate_q0 import generate_q0
# Use to switch on or off some tests
T = True
F = False
do_test = [ T, T, T, T, T, T, T, T ]
#
nu = 0.005
Nx = 16
Ny = 22
Nz = 19
N_fine   = 200
N_coarse = 100
Niter = 3
Tend  = 0.2
do_io = T
be_verbose = F
Np    = 2
run_cmd  = 'mpirun'
Np_s = '%0.2i' % (Np-1)
Np_s_p1 = '%0.2i' % (Np)

os.system('cp ../../bin/run_timestepper.out .')
os.system('cp ../../bin/run_parareal_mpi.out .')
os.system('cp ../../bin/parareal_openmp.out .')
os.system('cp ../../bin/parareal_openmp_pipe.out .')

#    
# (A) - Test that Parareal with Nit=0 reproduces the coarse integrator
#

# Test 0
test = 0
if do_test[test]:
    os.system('rm -f *.dat')
    # MPI version
    generate_q0(Nx, Ny, Nz)
    build_namelist(nu, Nx, Ny, Nz, Np*N_fine, Np*N_coarse, Niter, Tend, do_io, be_verbose)
    os.system('./run_timestepper.out C')
    # Run Parareal with Niter=0
    build_namelist(nu, Nx, Ny, Nz,    N_fine,    N_coarse,     0, Tend, do_io, be_verbose)
    os.system(run_cmd+' -n '+str(Np)+' ./run_parareal_mpi.out')
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
        sys.exit("ERROR: Parareal-MPI with Nit=0 and coarse integrator did not yield identical results.")
    elif numpy.isnan(max_err):
        sys.exit("ERROR: Parareal-MPI with Nit=0 and coarse integrator produced NaN error.")
else:
    print 'WARNING: Test '+str(test)+' switched off (Parareal-MPI=coarse for Nit=0)'
    
# Test 1    
test += 1    
if do_test[test]:
    os.system('rm -f *.dat')
    # OpenMP version
    generate_q0(Nx, Ny, Nz)
    build_namelist(nu, Nx, Ny, Nz, Np*N_fine, Np*N_coarse, Niter, Tend, do_io, be_verbose)
    os.system('./run_timestepper.out C')
    # Run Parareal with Niter=0
    build_namelist(nu, Nx, Ny, Nz,    N_fine,    N_coarse,     0, Tend, do_io, be_verbose)
    os.system('OMP_NUM_THREADS='+str(Np)+' ./parareal_openmp.out')
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
        sys.exit("ERROR: Parareal-OpenMP with Nit=0 and coarse integrator did not yield identical results.")
    elif numpy.isnan(max_err):
        sys.exit("ERROR: Parareal-OpenMP with Nit=0 and coarse integrator produced NaN error.")
else:
    print 'WARNING: Test '+str(test)+' switched off (Parareal-OpenMP = coarse for Nit=0)'

# Test 2
test += 1        
if do_test[test]:
    os.system('rm -f *.dat')
    # OpenMP-pipe version
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
        sys.exit("ERROR: Parareal-OpenMP-pipe with Nit=0 and coarse integrator did not yield identical results.")
    elif numpy.isnan(max_err):
        sys.exit("ERROR: Parareal-OpenMP-pipe with Nit=0 and coarse integrator produced NaN error.")
else:
    print 'WARNING: Test '+str(test)+' switched off (Parareal-OpenMP-pipe=coarse for Nit=0)'
        
#
# (B) Now test that Parareal with Nit=Nproc replicates fine integrator
#

# Test 3
test += 1
if do_test[test]:
# MPI version
    os.system('rm -f *.dat')
    generate_q0(Nx, Ny, Nz)
    build_namelist(nu, Nx, Ny, Nz, Np*N_fine, Np*N_coarse, Niter, Tend, do_io, be_verbose)
    os.system('./run_timestepper.out F')    
    # run parareal with Niter=Np
    build_namelist(nu, Nx, Ny, Nz,    N_fine,    N_coarse,    Np, Tend, do_io, be_verbose)
    os.system(run_cmd+' -n '+str(Np)+' ./run_parareal_mpi.out')
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
else:
    print 'WARNING: Test '+str(test)+' switched off (Parareal-MPI=Fine for Nit=Np)'
  
# Test 4    
test += 1    
if do_test[test]:
# OpenMP version
    os.system('rm -f *.dat')
    generate_q0(Nx, Ny, Nz)
    build_namelist(nu, Nx, Ny, Nz, Np*N_fine, Np*N_coarse, Niter, Tend, do_io, be_verbose)
    os.system('./run_timestepper.out F')    
    # run parareal with Niter=Np
    build_namelist(nu, Nx, Ny, Nz,    N_fine,    N_coarse,    Np, Tend, do_io, be_verbose)
    os.system('OMP_NUM_THREADS='+str(Np)+' ./parareal_openmp.out')
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
else:
    print 'WARNING: Test '+str(test)+' switched off (Parareal-OpenMP=Fine for Nit=Np)'
    
# Test 5
test += 1
if do_test[test]:
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
        print max_err
        sys.exit("ERROR: Parareal-OpenMP-pipe with Nit=Nproc and fine integrator did not yield identical results.")
    elif numpy.isnan(max_err):
        sys.exit("ERROR: Parareal-OpenMP-pipe with Nit=Nproc and fine integrator produced NaN error.")
else:
    print 'WARNING: Test '+str(test)+' switched off (Parareal-OpenMP-pipe=Fine for Nit=Np'
     
#
# (C) Make sure that OpenMP, OpenMP-pipe and MPI version all give identical results
#

# Test 6
test += 1
if do_test[test]:
    os.system('rm -f *.dat')
    generate_q0(Nx, Ny, Nz)
    build_namelist(nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose)
    os.system(run_cmd+' -n '+str(Np)+' ./run_parareal_mpi.out')
    os.system('OMP_NUM_THREADS='+str(Np)+' ./parareal_openmp.out')
    # add OpenMP-pipe
    Np_s = '%0.2i' % (Np)
    for nt in range(0,Np):
        nt_s = '%0.2i' % (nt)
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
            print max_err
            sys.exit('ERROR: Parareal-MPI and Parareal-OpenMP did not yield identical results.')
        elif numpy.isnan(max_err):
            sys.ext('ERROR: Parareal-MPI and Parareal-OpenMP produced NaN error')
else:   
    print 'WARNING: Test '+str(test)+' switched off (Parareal-MPI=Parareal-OpenMP)'

# Test 7
test += 1
if do_test[test]:
    os.system('rm -f *.dat')
    generate_q0(Nx, Ny, Nz)
    build_namelist(nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose)
    os.system(run_cmd+' -n '+str(Np)+' ./run_parareal_mpi.out')
    os.system('OMP_NUM_THREADS='+str(Np)+' ./parareal_openmp_pipe.out')
    Np_s = '%0.2i' % (Np)
    # add OpenMP-pipe
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
else:   
    print 'WARNING: Test '+str(test)+' switched off (Parareal-MPI=Parareal-OpenMP-pipe'
    
print "[0] -- Test of Parareal-MPI and Parareal-OpenMP successful"