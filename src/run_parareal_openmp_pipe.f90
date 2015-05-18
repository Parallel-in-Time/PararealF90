PROGRAM run_parareal_openmp

USE parareal_openmp_pipe, only: InitializePararealOpenMP_Pipe, FinalizePararealOpenMP_Pipe, PararealOpenMP_Pipe
USE params, only : Nx, Ny, Nz, dx, dy, dz, nu, N_coarse, N_fine, Niter, Tend, do_io, be_verbose, ReadParameter

IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q

! -- CODE: --

CALL ReadParameter()

! Initialize
CALL InitializePararealOpenMP_Pipe(nu, Nx, Ny, Nz)

! Load initial data
ALLOCATE(Q(-2:Nx+3,-2:Ny+3,-2:Nz+3))
OPEN(unit=20, FILE='q0.dat', ACTION='read', STATUS='old')
READ(20, '(F35.25)') Q
CLOSE(20)

CALL PararealOpenMP_Pipe(Q, Tend, N_fine, N_coarse, Niter, dx, dy, dz, do_io, be_verbose)

! Finalize
CALL FinalizePararealOpenMP_Pipe;


END PROGRAM run_parareal_openmp