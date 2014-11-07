! 
!
!
PROGRAM parareal_openmp

USE omp_lib
USE timestepper, only : Euler, Rk3Ssp, InitializeTimestepper, FinalizeTimestepper

IMPLICIT NONE

INTEGER, PARAMETER :: order_adv_c = 1, order_diff_c = 2, order_adv_f = 5, order_diff_f = 4

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q0
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, GQ, Qend

DOUBLE PRECISION :: nu, dx, dy, dz, Tend, dt_coarse, dt_fine, dt_slice, timer_all, &
        T0, T1, timer_coarse, timer_comm

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tstart_myslice, tend_myslice, &
    timer_fine
    
INTEGER :: Nx, Ny, Nz, k, N_fine, N_coarse, Niter, Nthreads, nt
     
LOGICAL :: do_io, be_verbose

CHARACTER(len=64) :: filename

NAMELIST /param/ nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose

! -------------- here comes the code....

OPEN(unit=20, FILE='parameter.in', ACTION='read', STATUS='old')
READ(20,NML=param)
CLOSE(20)

timer_all = OMP_GET_WTIME()
timer_coarse = 0.0
timer_comm   = 0.0

Nthreads = OMP_GET_MAX_THREADS()

ALLOCATE(timer_fine(0:Nthreads-1))
ALLOCATE(tstart_myslice(0:Nthreads-1))
ALLOCATE(tend_myslice(0:Nthreads-1))

IF (be_verbose) THEN
    WRITE(*,'(A, I2)') '--- Running OpenMP parareal, no. of threads: ', Nthreads
END IF

CALL InitializeTimestepper(nu, Nx, Ny, Nz, Nthreads)

dy = DBLE(1.0)/DBLE(Ny)
dz = DBLE(1.0)/DBLE(Nz)
dx = DBLE(1.0)/DBLE(Nx)

! Divide time interval [0,T] into Nproc many timeslices
dt_slice  = Tend/DBLE(Nthreads)
dt_fine   = dt_slice/DBLE(N_fine)
dt_coarse = dt_slice/DBLE(N_coarse)

ALLOCATE(Q(    -2:Nx+3,-2:Ny+3,-2:Nz+3, 0:Nthreads-1))
ALLOCATE(Qend( -2:Nx+3,-2:Ny+3,-2:Nz+3, 0:Nthreads-1))
ALLOCATE(GQ(   -2:Nx+3,-2:Ny+3,-2:Nz+3, 0:Nthreads-1))
ALLOCATE(Q0(   -2:Nx+3,-2:Ny+3,-2:Nz+3))
Q0(:,:,:) = 0.0
!$OMP PARALLEL DO schedule(static)
DO nt=0,Nthreads-1
    Q(:,:,:,nt)    = 0.0
    Qend(:,:,:,nt) = 0.0
    GQ(:,:,:,nt)   = 0.0
    timer_fine(nt)   = 0.0
    tstart_myslice(nt) = DBLE(nt)*dt_slice
    tend_myslice(nt)   = DBLE(nt+1)*dt_slice

    IF (be_verbose) THEN
        WRITE(*,'(A, I2, A, F5.3, A, F5.3)') 'Thread ', nt, ' from ', tstart_myslice(nt), ' to ', tend_myslice(nt)
    END IF

END DO
!$OMP END PARALLEL DO

dt_slice  = Tend/DBLE(Nthreads)
dt_fine   = dt_slice/DBLE(N_fine)
dt_coarse = dt_slice/DBLE(N_coarse)
IF(be_verbose) THEN
    WRITE(*,'(A, F9.5)') 'Time slice length:  ', dt_slice
    WRITE(*,'(A, F9.5)') 'Fine step length:   ', dt_fine
    WRITE(*,'(A, F9.5)') 'Coarse step length: ', dt_coarse   
END IF

! Set initial value:
! Q(0) <- y^0_0
OPEN(unit=20, FILE='q0.dat', ACTION='read', STATUS='old')
READ(20, '(F35.25)') Q0
CLOSE(20)
Q(:,:,:,0) = Q0

T0 = OMP_GET_WTIME()

DO nt=0,Nthreads-1

    ! GQ(nt) <- Q(nt) = y^0_nt
    GQ(:,:,:,nt) = Q(:,:,:,nt)
    
    ! GQ(nt) <- G(y^0_nt)
    CALL Euler(GQ(:,:,:,nt), tstart_myslice(nt), tend_myslice(nt), N_coarse, dx, dy, dz, order_adv_c, order_diff_c)
    
    ! Q(nt+1) <- y^0_(nt+1) = G(y^0_nt), if not last slice
    IF (nt<Nthreads-1) THEN
        Q(:,:,:,nt+1) = GQ(:,:,:,nt)
    END IF

    
    ! Q(nt+1) <- y^0_nt+1 = G(y^0_nt)
    ! GQ(nt)  <- G(y^0_nt)
END DO

! timer_coarse(0) = timer_coarse(0) + OMP_GET_WTIME() - T0 
! For some reason, can't write timer_coarse(0), but need 0:0

! Right now, the coarse guess is the closest thing to a solution
! Qend(nt) <- G(y^0_nt)
Qend = GQ

timer_coarse = timer_coarse + OMP_GET_WTIME() - T0

DO k=1,Niter

    ! Initial state
    ! Q(nt)    <- y^(k-1)_nt
    ! GQ(nt)   <- G(y^(k-1)_nt)
    ! Qend(nt) <- y^(k-1)_(nt+1)
    
    ! Run fine integrator
    ! Q(nt) <- F(y^(k-1)_nt)
    
    !$OMP PARALLEL DO schedule(static) private(T0, T1)
    DO nt=0,Nthreads-1
        T0 = OMP_GET_WTIME()
        CALL Rk3Ssp(Q(:,:,:,nt), tstart_myslice(nt), tend_myslice(nt), N_fine, dx, dy, dz, order_adv_f, order_diff_f)
        T1 = OMP_GET_WTIME()
        
        ! Compute difference fine minus coarse
        ! Qend(nt) <- F(y^(k-1)_nt) - G(y^(k-1)_nt)
        Qend(:,:,:,nt) = Q(:,:,:,nt) - GQ(:,:,:,nt)
                    
        timer_fine(nt) = timer_fine(nt) + T1 - T0
    END DO
    !$OMP END PARALLEL DO
    
    ! Now read initial value into Q(0) again:
    ! Q(0) <- y^k_0
    T0 = OMP_GET_WTIME()
    Q(:,:,:,0) = Q0
    timer_comm = timer_comm + OMP_GET_WTIME() - T0
    
    DO nt=0,Nthreads-1

        T0 = OMP_GET_WTIME()        
        ! GQ(nt) <- y^k_nt
        GQ(:,:,:,nt) = Q(:,:,:,nt)
        
        ! GQ(nt) <- G(y^k_nt)
        CALL Euler(GQ(:,:,:,nt), tstart_myslice(nt), tend_myslice(nt), N_coarse, dx, dy, dz, order_adv_c, order_diff_c)
        timer_coarse = timer_coarse + OMP_GET_WTIME() - T0
        
        ! Qend(nt) <- y^k_(nt+1) = G(y^k_nt) + F(y^(k-1)_nt) - G(y^(k-1)_nt)
        Qend(:,:,:,nt) = GQ(:,:,:,nt) + Qend(:,:,:,nt)
        
        ! Now update initial value for next slice
        ! Q(nt+1) = Qend(nt) <- y^k_(n+1)
        IF (nt<Nthreads-1) THEN
            T0 = OMP_GET_WTIME()
            Q(:,:,:,nt+1) = Qend(:,:,:,nt)
            timer_comm = timer_comm + OMP_GET_WTIME() - T0
        END IF
    END DO
END DO

IF(do_io) THEN
    DO nt=0,Nthreads-1
        WRITE(filename, '(A,I0.2,A,I0.2,A)') 'q_final_', nt, '_', Nthreads, '_openmp.dat'
        OPEN(unit=nt, FILE=filename, ACTION='write', STATUS='replace')
        WRITE(nt, '(F35.25)') Qend(1:Nx, 1:Ny, 1:Nz, nt)
        CLOSE(nt)
    END DO
END IF
CALL FinalizeTimestepper()

timer_all = OMP_GET_WTIME() - timer_all

IF(do_io) THEN
    DO nt=0,Nthreads-1
        WRITE(filename, '(A,I0.2,A,I0.2,A)') 'timings_openmp', nt, '_', Nthreads, '.dat'
        OPEN(unit=nt, FILE=filename, ACTION='write', STATUS='replace')
        WRITE(nt, '(F8.2)') timer_all
        WRITE(nt, '(F8.2)') timer_fine(nt)
        WRITE(nt, '(F8.2)') timer_coarse
        WRITE(nt, '(F8.2)') timer_comm
        CLOSE(nt)
    END DO
END IF

END PROGRAM parareal_openmp