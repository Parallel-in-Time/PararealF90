PROGRAM run_parareal_openmp

USE parareal_openmp_pipe, only: InitializePararealOpenMP_Pipe, FinalizePararealOpenMP_Pipe, PararealOpenMP_Pipe

IMPLICIT NONE

!> @todo docu
LOGICAL :: do_io, be_verbose

DOUBLE PRECISION :: nu, dx, dy, dz, Tend

INTEGER :: Nx, Ny, Nz, N_fine, N_coarse, Niter

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q

!> @todo docu
NAMELIST /param/ nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose

! -- CODE: --

! Read parameters
OPEN(unit=20, FILE='parameter.in', ACTION='read', STATUS='old')
READ(20,NML=param)
CLOSE(20)

! Computational domain is unit cube
dy = DBLE(1.0)/DBLE(Ny)
dz = DBLE(1.0)/DBLE(Nz)
dx = DBLE(1.0)/DBLE(Nx)

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