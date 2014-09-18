PROGRAM parareal_openmp

USE omp_lib
USE timestepper, only : Euler, Rk3Ssp, InitializeTimestepper, FinalizeTimestepper

IMPLICIT NONE

INTEGER, PARAMETER :: order_adv_c = 1, order_diff_c = 2, order_adv_f = 5, order_diff_f = 4

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q0
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, GQ, Qend

DOUBLE PRECISION :: nu, dx, dy, dz, Tend, dt_coarse, dt_fine, dt_slice, timer_all, &
        T0, T1 

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tstart_myslice, tend_myslice, &
    timer_fine, timer_coarse, timer_comm
    
INTEGER :: Nx, Ny, Nz, k, N_fine, N_coarse, Niter, Nthreads, nt, thread_nr
INTEGER(kind=OMP_LOCK_KIND), DIMENSION(:), ALLOCATABLE :: nlocks

LOGICAL :: do_io, be_verbose

CHARACTER(len=64) :: filename

NAMELIST /param/ nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose

! -------------- here comes the code....

OPEN(unit=20, FILE='parameter.in', ACTION='read', STATUS='old')
READ(20,NML=param)
CLOSE(20)

timer_all = OMP_GET_WTIME()

Nthreads = OMP_GET_MAX_THREADS()
!call omp_set_num_threads(1)

ALLOCATE(timer_fine(0:Nthreads-1))
ALLOCATE(timer_coarse(0:Nthreads-1))
ALLOCATE(timer_comm(0:Nthreads-1))
ALLOCATE(tstart_myslice(0:Nthreads-1))
ALLOCATE(tend_myslice(0:Nthreads-1))
ALLOCATE(nlocks(0:Nthreads-1))

IF (be_verbose) THEN
    WRITE(*,'(A, I2)') '--- Running OpenMP-pipe parareal, no. of threads: ', Nthreads
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
    timer_coarse(nt) = 0.0
    timer_comm(nt)   = 0.0
    tstart_myslice(nt) = DBLE(nt)*dt_slice
    tend_myslice(nt)   = DBLE(nt+1)*dt_slice

    IF (be_verbose) THEN
        WRITE(*,'(A, I2, A, F5.3, A, F5.3)') 'Thread ', nt, ' from ', tstart_myslice(nt), ' to ', tend_myslice(nt)
    END IF

    ! Initialize locks, one for each timeslice
    CALL OMP_INIT_LOCK(nlocks(nt))
    
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

!$OMP PARALLEL private(thread_nr, k, T0, T1)

thread_nr = omp_get_thread_num()

!$OMP DO schedule(static)
DO nt=0,Nthreads-1

    T0 = OMP_GET_WTIME()

    Q(:,:,:,nt) = Q(:,:,:,0)

    IF (nt>0) THEN

        ! Q(nt)  <- y^0_nt = G(y^0_(nt-1))
        CALL Euler(Q(:,:,:,nt), DBLE(0.0), tstart_myslice(nt), nt*N_coarse, dx, dy, dz, order_adv_c, order_diff_c)

    END IF

    ! GQ(nt) <- G(y^0_nt)
    GQ(:,:,:,nt) = Q(:,:,:,nt)

    CALL Euler(GQ(:,:,:,nt), tstart_myslice(nt), tend_myslice(nt), N_coarse, dx, dy, dz, order_adv_c, order_diff_c)

    ! Qend(nt) <- G(y^0_nt) = y^0_(nt+1)    
    Qend(:,:,:,nt) = GQ(:,:,:,nt)

    timer_coarse(thread_nr) = timer_coarse(thread_nr) + OMP_GET_WTIME() - T0    

END DO
!$OMP END DO NOWAIT

DO k=1,Niter

    ! Initial state
    ! Q(nt)    <- y^(k-1)_nt
    ! GQ(nt)   <- G(y^(k-1)_nt)
    ! Qend(nt) <- y^(k-1)_(nt+1)

    ! Q(nt) <- F(y^(k-1)_nt)

    !$OMP DO schedule(static) ORDERED
    DO nt=0,Nthreads-1

        T0 = OMP_GET_WTIME()

        ! Lock part of Q corresponding to time-slice nt while writing
        CALL OMP_SET_LOCK(nlocks(nt))

        CALL RK3Ssp(Q(:,:,:,nt), tstart_myslice(nt), tend_myslice(nt), N_fine, dx, dy, dz, order_adv_f, order_diff_f)

        ! Qend(nt) <- F(y^(k-1)_nt) - G(y^(k-1)_nt)
        Qend(:,:,:,nt) = Q(:,:,:,nt) - GQ(:,:,:,nt)

        CALL OMP_UNSET_LOCK(nlocks(nt))

        timer_fine(thread_nr) = timer_fine(thread_nr) + OMP_GET_WTIME() - T0

        !$OMP ORDERED

        IF (nt==0) THEN

            ! Set Q(0) to initial value again
            T0 = OMP_GET_WTIME()

            CALL OMP_SET_LOCK(nlocks(0))
            Q(:,:,:,0) = Q0
            CALL OMP_UNSET_LOCK(nlocks(0))

            timer_comm(0) = timer_comm(0) + OMP_GET_WTIME() - T0
            
        END IF

        T0 = OMP_GET_WTIME()

        ! GQ(nt) <- y^k_nt
        GQ(:,:,:,nt) = Q(:,:,:,nt)

        ! GQ(nt) <- G(y^k_nt)
        CALL Euler( GQ(:,:,:,nt), tstart_myslice(nt), tend_myslice(nt), N_coarse, dx, dy, dz, order_adv_c, order_diff_c)

        ! Qend(nt) <- y^k_(nt+1) G(y^k_nt) + F(y^(k-1)_nt) - G(y^(k-1)_nt)
        Qend(:,:,:,nt) = GQ(:,:,:,nt) + Qend(:,:,:,nt)

        timer_coarse(thread_nr) = timer_coarse(thread_nr) + OMP_GET_WTIME() - T0

        IF (nt<Nthreads-1) THEN
            ! Q(nt+1) = Qend(nt) = y^k_(nt+1)
            T0 = OMP_GET_WTIME()

            ! Update value in Q(nt+1), set lock to avoid starting updating it while
            ! thread nt+1 is still busy writing the update fine value in it (in call to F above)
            CALL OMP_SET_LOCK(nlocks(nt+1))
            Q(:,:,:,nt+1) = Qend(:,:,:,nt)
            CALL OMP_UNSET_LOCK(nlocks(nt+1))
        
            timer_comm(thread_nr) = timer_comm(thread_nr) + OMP_GET_WTIME() - T0
        END IF

        !$OMP END ORDERED

    END DO
    !$OMP END DO NOWAIT

END DO
!$OMP END PARALLEL

IF(do_io) THEN
    DO nt=0,Nthreads-1
        CALL OMP_DESTROY_LOCK(nlocks(nt))
        IF (nt<10) THEN
          WRITE(filename, '(A,I1,A,I2,A)') 'q_final_0', nt, '_', Nthreads, '_openmp_pipe.dat'
        ELSE
          WRITE(filename, '(A,I2,A,I2,A)') 'q_final_', nt, '_', Nthreads, '_openmp_pipe.dat'
        END IF 
        OPEN(unit=nt, FILE=filename)
        WRITE(nt, '(F35.25)') Qend(1:Nx, 1:Ny, 1:Nz, nt)
        CLOSE(nt)
    END DO
END IF
CALL FinalizeTimestepper()

timer_all = OMP_GET_WTIME() - timer_all

IF(do_io) THEN
    DO nt=0,Nthreads-1
        IF (nt<10) THEN
          WRITE(filename, '(A,I1,A,I2,A)') 'timings_openmp_pipe_0', nt, '_', Nthreads, '.dat'
        ELSE
          WRITE(filename, '(A,I2,A,I2,A)') 'timings_openmp_pipe', nt, '_', Nthreads, '.dat'
        END IF 
        OPEN(unit=nt, FILE=filename)
        WRITE(nt, '(F8.2)') timer_all
        WRITE(nt, '(F8.2)') timer_fine(nt)
        WRITE(nt, '(F8.2)') timer_coarse(nt)
        WRITE(nt, '(F8.2)') timer_comm(nt)
        CLOSE(nt)
    END DO
END IF

END PROGRAM parareal_openmp