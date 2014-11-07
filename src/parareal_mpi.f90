!>
!! @todo docu
!!
!! \\( y^{k+1}_{n+1} = G(y^{k+1}_{n}) + F(y^k_n) - G(y^k_n) \\)
MODULE parareal_mpi

USE omp_lib
USE timestepper, only : Euler, Rk3Ssp, InitializeTimestepper, FinalizeTimestepper

IMPLICIT NONE

INCLUDE 'mpif.h'

!> @todo docu
INTEGER, PARAMETER :: Nthreads = 1, & ! In MPI version, don't need multiple threads
    order_adv_c = 1, order_diff_c = 2, order_adv_f = 5, order_diff_f = 4

!> @todo docu
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q, GQ, Qend

!> @todo docu
DOUBLE PRECISION :: dt_coarse

!> @todo docu
DOUBLE PRECISION :: dt_fine

!> @todo docu
DOUBLE PRECISION :: dt_slice

!> @todo docu
DOUBLE PRECISION :: tstart_myslice

!> @todo docu
DOUBLE PRECISION :: tend_myslice

! Internal variables to be used for timers.
DOUBLE PRECISION :: timer_coarse, timer_fine, timer_comm, timer_all, T0, T1

!> @todu docu
INTEGER :: mpi_thread_provided, dim, ierr, k, Nproc, myrank, recv_status(MPI_STATUS_SIZE), send_status(MPI_STATUS_SIZE)

!> @todo docu
CHARACTER(len=64) :: filename


TYPE parareal_parameter
    INTEGER Nx, Ny, Nz
END TYPE

TYPE(parareal_parameter) :: param

CONTAINS

  !> @todo docu
  SUBROUTINE InitializePararealMPI(nu, Nx, Ny, Nz)
    DOUBLE PRECISION, INTENT(IN) :: nu
    INTEGER, INTENT(IN) :: Nx, Ny, Nz

    CALL InitializeTimestepper(nu, Nx, Ny, Nz, Nthreads)
    ALLOCATE(Q(   -2:Nx+3,-2:Ny+3,-2:Nz+3))
    ALLOCATE(Qend(-2:Nx+3,-2:Ny+3,-2:Nz+3))
    ALLOCATE(GQ(  -2:Nx+3,-2:Ny+3,-2:Nz+3))

    dim = SIZE(Q,1)*SIZE(Q,2)*SIZE(Q,3)

    param%Nx = Nx
    param%Ny = Ny
    param%Nz = Nz

  END SUBROUTINE InitializePararealMPI

  !> @todo docu
  SUBROUTINE FinalizePararealMPI()
    DEALLOCATE(Q)
    DEALLOCATE(GQ)
    DEALLOCATE(Qend)
    CALL FinalizeTimestepper()
  END SUBROUTINE FinalizePararealMPI

  !> @todo docu
  SUBROUTINE PararealMPI(Q_initial, Tend, N_fine, N_coarse, Niter, dx, dy, dz, do_io, be_verbose)

    DOUBLE PRECISION, DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: Q_initial
    DOUBLE PRECISION, INTENT(IN) :: Tend, dx, dy, dz
    INTEGER,                                  INTENT(IN)    :: N_fine, N_coarse, Niter
    LOGICAL,                                  INTENT(IN)    :: do_io, be_verbose

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nproc, ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    ! Output
    IF ((myrank==0) .AND. (be_verbose)) THEN
      WRITE(*,'(A, I2)') '--- Running MPI parareal, no. of processes: ', Nproc
    END IF

    timer_all    = MPI_WTIME()
    timer_fine   = 0.0
    timer_coarse = 0.0
    timer_comm   = 0.0

    ! Divide time interval [0,T] into Nproc many timeslices
    dt_slice  = Tend/DBLE(Nproc)
    dt_fine   = dt_slice/DBLE(N_fine)
    dt_coarse = dt_slice/DBLE(N_coarse)

    tstart_myslice = DBLE(myrank)*dt_slice
    tend_myslice   = DBLE(myrank+1)*dt_slice

    IF (be_verbose) THEN
        WRITE(*,'(A, I2, A, F5.3, A, F5.3)') 'Process ', myrank, ' from ', tstart_myslice, ' to ', tend_myslice
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
        IF (myrank==0) THEN
            WRITE(*,'(A, F9.5)') 'Time slice length:  ', dt_slice
            WRITE(*,'(A, F9.5)') 'Fine step length:   ', dt_fine
            WRITE(*,'(A, F9.5)') 'Coarse step length: ', dt_coarse
        END IF
    END IF

    ! --- START PARAREAL --- !

    T0 = MPI_WTIME()

    ! Q <- y^0_0
    Q = Q_initial

    IF (myrank>0) THEN
        ! coarse predictor to get initial guess at beginning of slice
        ! Q <- y^0_n = G(y^0_(n-1)) = G(G(y^0_(n-2))) = ...
        CALL Euler(Q, DBLE(0.0), tstart_myslice, myrank*N_coarse, dx, dy, dz, order_adv_c, order_diff_c)
    END IF

    ! one more coarse step to get coarse value at end of slice
    ! GQ <- G(y^0_n)
    GQ = Q

    CALL Euler(GQ, tstart_myslice, tend_myslice, N_coarse, dx, dy, dz, order_adv_c, order_diff_c)

    ! Qend <- G(y^0_n) = y^0_(n+1)
    Qend = GQ

    T1 = MPI_WTIME()
    timer_coarse = timer_coarse + (T1-T0)

    DO k=1,Niter   ! Parareal iteration
         
      ! Initial state: 
      ! Q     <- y^(k-1)_n 
      ! GQ    <- G(y^(k-1)_n) 
      ! Qend  <- y(k-1)_(n+1)
      !
      ! End state:
      ! Q     <- y^k_n
      ! GQ    <- G(y^k_n)
      ! Qend  <- y^k_(n+1) 
           

      T0 = MPI_WTIME()

      ! Run fine integrator:
      ! Q <- F(y^(k-1)_n)
      CALL Rk3SSp(Q, tstart_myslice, tend_myslice, N_fine, dx, dy, dz, order_adv_f, order_diff_f)
      
      ! Compute difference fine minus coarse
      ! Qend <- F(y^(k-1)_n) - G(y^(k-1)_n)
      Qend = Q - GQ
      
       T1 = MPI_WTIME()
       timer_fine = timer_fine + (T1-T0)

      ! Fetch updated value from previous process
      IF(myrank>0) THEN
        
        ! Q <- y^k_n 
        T0 = MPI_WTIME()
        CALL MPI_RECV(Q, dim, MPI_DOUBLE_PRECISION, myrank-1, 0, MPI_COMM_WORLD, recv_status, ierr)
        T1 = MPI_WTIME()
        timer_comm = timer_comm + (T1-T0)
           
      ELSE IF (myrank==0) THEN

        ! Fetch initial value again
        ! Q <- y^k_n with n=0
        T0 = MPI_WTIME()
        Q = Q_initial
        T1 = MPI_WTIME()
        timer_comm = timer_comm + (T1-T0)

      ELSE
        WRITE(*,*) 'Found negative value for myrank, now exiting.'
        STOP  
      END IF
      
      ! GQ <- y^k_n
      GQ = Q
      
      ! Run correction
      ! GQ <- G(y^k_n)
      T0 = MPI_WTIME()
      CALL Euler(GQ, tstart_myslice, tend_myslice, N_coarse, dx, dy, dz, order_adv_c, order_diff_c)
      T1 = MPI_WTIME()
       
      ! Correct
      ! Qend <- G(y^k_n) + F(y^(k-1))_n - G(y^(k-1)_n) = y^k_(n+1)
      Qend = GQ + Qend

      timer_coarse = timer_coarse + (T1-T0)
        
      ! Send forward updated value
      IF(myrank<Nproc-1) THEN
        T0 = MPI_WTIME()
        CALL MPI_SEND(Qend, dim, MPI_DOUBLE_PRECISION, myrank+1, 0, MPI_COMM_WORLD, send_status, ierr)
        T1 = MPI_WTIME()
        timer_comm = timer_comm + (T1-T0)
      END IF 
      
      ! Final state:
      ! Q    <- y^k_n
      ! GQ   <- G(y^k_n)
      ! Qend <- y^k_(n+1)
      
    END DO ! End of parareal loop

    ! --- END PARAREAL --- !

    ! Return final value in Q_initial
    Q_initial = Qend

    IF(do_io) THEN
        WRITE(filename, '(A,I0.2,A,I0.2,A)') 'q_final_', myrank, '_', Nproc, '_mpi.dat'
        OPEN(UNIT=myrank, FILE=filename, ACTION='write', STATUS='replace')
        WRITE(myrank, '(F35.25)') Qend(1:param%Nx, 1:param%Ny, 1:param%Nz)
        CLOSE(myrank)
    END IF

    timer_all = MPI_WTIME() - timer_all

    IF(do_io) THEN
        WRITE(filename, '(A,I0.2,A,I0.2,A)') 'timings_mpi', myrank, '_', Nproc, '.dat'
        OPEN(UNIT=myrank, FILE=filename, ACTION='write', STATUS='replace')
        WRITE(myrank, '(F8.2)') timer_all
        WRITE(myrank, '(F8.2)') timer_fine
        WRITE(myrank, '(F8.2)') timer_coarse
        WRITE(myrank, '(F8.2)') timer_comm
        CLOSE(myrank)
    END IF

  END SUBROUTINE PararealMPI

END MODULE parareal_mpi