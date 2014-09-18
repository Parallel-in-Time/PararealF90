PROGRAM parareal_mpi

USE omp_lib
USE timestepper, only : Euler, Rk3Ssp, InitializeTimestepper, FinalizeTimestepper

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER :: Nthreads = 1, & ! In MPI version, don't need multiple threads
    order_adv_c = 1, order_diff_c = 2, order_adv_f = 5, order_diff_f = 4

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q, GQ, Qend

DOUBLE PRECISION :: nu, dx, dy, dz, Tend, dt_coarse, dt_fine, dt_slice, &
    tstart_myslice, tend_myslice, timer_coarse, timer_fine, timer_comm, timer_all, T0, T1
    
INTEGER :: mpi_thread_provided, dim, ierr, Nx, Ny, Nz, k, N_fine,&
     N_coarse, Nproc, myrank, Niter, recv_status(MPI_STATUS_SIZE), send_status(MPI_STATUS_SIZE)
     
LOGICAL :: do_io, be_verbose

CHARACTER(len=64) :: filename

NAMELIST /param/ nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose

! -------------- here comes the code....

OPEN(unit=20, FILE='parameter.in', ACTION='read', STATUS='old')
READ(20,NML=param)
CLOSE(20)

CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, Nproc, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

timer_all    = MPI_WTIME()
timer_fine   = 0.0
timer_coarse = 0.0
timer_comm   = 0.0

IF ((myrank==0) .AND. (be_verbose)) THEN
    WRITE(*,'(A, I2)') '--- Running MPI parareal, no. of processes: ', Nproc
END IF

CALL InitializeTimestepper(nu, Nx, Ny, Nz, Nthreads)
ALLOCATE(Q(    -2:Nx+3,-2:Ny+3,-2:Nz+3))
ALLOCATE(Qend(-2:Nx+3,-2:Ny+3,-2:Nz+3))
ALLOCATE(GQ(   -2:Nx+3,-2:Ny+3,-2:Nz+3))

dim = SIZE(Q,1)*SIZE(Q,2)*SIZE(Q,3)

dy = DBLE(1.0)/DBLE(Ny)
dz = DBLE(1.0)/DBLE(Nz)
dx = DBLE(1.0)/DBLE(Nx)
    
OPEN(unit=20, FILE='q0.dat', ACTION='read', STATUS='old')
READ(20, '(F35.25)') Q
CLOSE(20)

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

! --- START ACTUAL PARAREAL ---

T0 = MPI_WTIME()

IF (myrank>0) THEN
    ! coarse predictor to get initial guess at beginning of slice
    ! Q <- y^0_n
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

! now iterate
DO k=1,Niter
     
  ! Initial state: 
  ! Q     <- y^(k-1)_n 
  ! GQ    <- G(y^(k-1)_n) 
  ! Qend  <- y(k-1)_(n+1)
  !
  ! End state:
  ! Q     <- y^k_n
  ! GQ    <- G(y^k_n)
  ! Qend  <- y^k_(n+1) 
       
  ! Run fine integrator:
  ! Q <- F(y^(k-1)_n)
  T0 = MPI_WTIME()
  
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
    OPEN(unit=20, FILE='q0.dat', ACTION='read', STATUS='old')
    READ(20, '(F35.25)') Q
    CLOSE(20)    
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
  
  END DO

IF(do_io) THEN
    IF (myrank<10) THEN
      WRITE(filename, '(A,I1,A,I2,A)') 'q_final_0', myrank, '_', Nproc,'_mpi.dat'
    ELSE
      WRITE(filename, '(A,I2,A,I2,A)') 'q_final_', myrank, '_', Nproc, '_mpi.dat'
    END IF 
    OPEN(unit=myrank, FILE=filename)
    WRITE(myrank, '(F35.25)') Qend(1:Nx, 1:Ny, 1:Nz)
    CLOSE(myrank)
END IF

! END ACTUAL PARAREAL

CALL FinalizeTimestepper

timer_all = MPI_WTIME() - timer_all

IF(do_io) THEN
    IF (myrank<10) THEN
      WRITE(filename, '(A,I1,A,I2,A)') 'timings_mpi_0', myrank, '_', Nproc, '.dat'
    ELSE
      WRITE(filename, '(A,I2,A,I2,A)') 'timings_mpi', myrank, '_', Nproc, '.dat'
    END IF 
    OPEN(unit=myrank, FILE=filename)
    WRITE(myrank, '(F8.2)') timer_all
    WRITE(myrank, '(F8.2)') timer_fine
    WRITE(myrank, '(F8.2)') timer_coarse
    WRITE(myrank, '(F8.2)') timer_comm
    CLOSE(myrank)
END IF

CALL MPI_FINALIZE(ierr)

END PROGRAM parareal_mpi