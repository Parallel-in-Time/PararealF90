PROGRAM run_parareal_mpi

USE parareal_mpi, only: InitializePararealMPI, FinalizePararealMPI, PararealMPI;

!> @todo docu
LOGICAL :: do_io, be_verbose

DOUBLE PRECISION :: nu, dx, dy, dz, Tend

INTEGER :: Nx, Ny, Nz, N_fine, N_coarse, Niter



!> @todo docu
NAMELIST /param/ nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose

! CODE:

OPEN(unit=20, FILE='parameter.in', ACTION='read', STATUS='old')
READ(20,NML=param)
CLOSE(20)



CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)


IF ((myrank==0) .AND. (be_verbose)) THEN
WRITE(*,'(A, I2)') '--- Running MPI parareal, no. of processes: ', Nproc
END IF


dim = SIZE(Q,1)*SIZE(Q,2)*SIZE(Q,3)

dy = DBLE(1.0)/DBLE(Ny)
dz = DBLE(1.0)/DBLE(Nz)
dx = DBLE(1.0)/DBLE(Nx)

OPEN(unit=20, FILE='q0.dat', ACTION='read', STATUS='old')
READ(20, '(F35.25)') Q
CLOSE(20)


CALL MPI_FINALIZE(ierr)


END PROGRAM run_parareal_mpi