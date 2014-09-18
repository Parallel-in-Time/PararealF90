PROGRAM WENO5_Weakscaling

USE Advection, only : GetRHSAdvection, InitializeAdvection, FinalizeAdvection
USE OMP_LIB

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER :: Nthreads = 8, Nx_factor = 10, Nx_base = 20
            
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, RQ 
DOUBLE PRECISION                                  :: dx, dy, dz, T0, T1, runtime, runtime_max
INTEGER                                           :: r, ierr, mpi_thread_provided

INTEGER :: Nx, Ny, Nz, i_end, j_end, i_start, j_start, k_start, k_end, order, ntasks, myrank

CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, ntasks, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

i_start = 1
j_start = 1
k_start = 1

order=5

DO r=Nx_factor,Nx_factor

    Nx = r*Nx_base
    Ny = r*Nx_base
    Nz = r*Nx_base
    dx = 1.0/Nx
    dy = 1.0/Ny
    dz = 1.0/Nz

    i_end = i_start+Nx-1
    j_end = j_start+Ny-1
    k_end = k_start+Nz-1

    ALLOCATE(Q( j_start-3:j_end+3,   i_start-3:i_end+3,  k_start-3:k_end+3, 0:0))
    ALLOCATE(RQ(j_start-3:j_end+3,   i_start-3:i_end+3,  k_start-3:k_end+3, 0:0))

    CALL InitializeAdvection(i_start-3, i_end+3, j_start-3, j_end+3, k_start-3, k_end+3, 1, .true.)

    T0 = MPI_WTIME()
    CALL GetRHSAdvection(Q(:,:,:,0), RQ(:,:,:,0), dx, dy, dz, i_start, i_end, j_start, j_end, k_start, k_end, order)
    T1 = MPI_WTIME()

    runtime = T1-T0

    DEALLOCATE(Q)
    DEALLOCATE(RQ)
    CALL FinalizeAdvection()

END DO

CALL MPI_REDUCE(runtime, runtime_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD)

IF (myrank==0) THEN
    WRITE(*,'(A, F5.3)') 'Runtime in seconds: ', runtime_max
END IF
CALL MPI_FINALIZE(ierr)

END PROGRAM WENO5_Weakscaling