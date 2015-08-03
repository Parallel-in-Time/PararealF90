PROGRAM run_parareal_mpi

USE parareal_mpi, only: InitializePararealMPI, FinalizePararealMPI, PararealMPI
USE params, only : Nx, Ny, Nz, dx, dy, dz, nu, N_coarse, N_fine, Niter, Tend, do_io, be_verbose, ReadParameter

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER :: ierr

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q


!> @todo docu

CALL ReadParameter()

! Initialize
CALL MPI_INIT(ierr)
CALL InitializePararealMPI(nu, Nx, Ny, Nz)

! Load initial data
ALLOCATE(Q(-2:Nx+3,-2:Ny+3,-2:Nz+3))
OPEN(unit=20, FILE='q0.dat', ACTION='read', STATUS='old')
READ(20, '(F35.25)') Q
CLOSE(20)

CALL PararealMPI(Q, Tend, N_fine, N_coarse, Niter, dx, dy, dz, do_io, be_verbose)

! Finalize
CALL FinalizePararealMPI;
CALL MPI_FINALIZE(ierr)


END PROGRAM run_parareal_mpi