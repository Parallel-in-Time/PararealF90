PROGRAM run_parareal_mpi

USE parareal_mpi, only: InitializePararealMPI, FinalizePararealMPI, PararealMPI;

IMPLICIT NONE

INCLUDE 'mpif.h'

!> @todo docu
LOGICAL :: do_io, be_verbose

DOUBLE PRECISION :: nu, dx, dy, dz, Tend

INTEGER :: Nx, Ny, Nz, N_fine, N_coarse, Niter, ierr, mpi_thread_provided

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q

!> @todo docu
NAMELIST /param/ nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose

! CODE:

! Read parameters
OPEN(unit=20, FILE='parameter.in', ACTION='read', STATUS='old')
READ(20,NML=param)
CLOSE(20)

! Computational domain is unit cube
dy = DBLE(1.0)/DBLE(Ny)
dz = DBLE(1.0)/DBLE(Nz)
dx = DBLE(1.0)/DBLE(Nx)

! Initialize
CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)
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