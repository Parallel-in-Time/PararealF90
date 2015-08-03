!>
!! A OpenMP based implementation of the parallel-in-time Parareal method with pipelining Here, solutions on all time slices are hold in shared memory so that a thread computing a time slice can 
!! fetch the update computed by its predecessor directly from the shared memory without the need for message passing. While in the non-pipelined version only the fine propagator is evaluated in
!! parallel, here the different threads work more like processes in the MPI implementation. OpenMP locks are set and released manually to ensure the necessary synchronization.
!!
MODULE parareal_openmp_pipe

USE omp_lib
USE timestepper, only : Euler, Rk3Ssp, InitializeTimestepper, FinalizeTimestepper

IMPLICIT NONE

PRIVATE
PUBLIC :: InitializePararealOpenMP_Pipe, FinalizePararealOpenMP_Pipe, PararealOpenMP_Pipe

!> Fixed orders of spatial discretization
INTEGER, PARAMETER :: order_adv_c = 1, order_diff_c = 2, order_adv_f = 5, order_diff_f = 4

!> Three solution buffers used in Parareal
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, GQ, Qend

!> Coarse time step
DOUBLE PRECISION :: dt_coarse

!> Fine time step
DOUBLE PRECISION :: dt_fine

!> Lenth of time slice
DOUBLE PRECISION :: dt_slice

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tstart_myslice, tend_myslice, &
    timer_fine, timer_coarse, timer_comm

! Internal timer variables
DOUBLE PRECISION :: timer_all, T0, T1

INTEGER :: k, Nthreads, nt, thread_nr

!> For the OpenMP version with pipelining, the communication/update steps need to be manually coordinated with OpenMP locks
INTEGER(kind=OMP_LOCK_KIND), DIMENSION(:), ALLOCATABLE :: nlocks

CHARACTER(len=64) :: filename

TYPE parareal_parameter
  INTEGER Nx, Ny, Nz
END TYPE

TYPE(parareal_parameter) :: param

CONTAINS

  !> Initialize and allocate needed buffers
  SUBROUTINE InitializePararealOpenMP_Pipe(nu, Nx, Ny, Nz)
    DOUBLE PRECISION, INTENT(IN) :: nu
    INTEGER, INTENT(IN) :: Nx, Ny, Nz

    Nthreads = OMP_GET_MAX_THREADS()

    param%Nx = Nx
    param%Ny = Ny
    param%Nz = Nz

    ALLOCATE(timer_fine(0:Nthreads-1))
    ALLOCATE(timer_coarse(0:Nthreads-1))
    ALLOCATE(timer_comm(0:Nthreads-1))
    ALLOCATE(tstart_myslice(0:Nthreads-1))
    ALLOCATE(tend_myslice(0:Nthreads-1))
    ALLOCATE(nlocks(0:Nthreads-1))

    CALL InitializeTimestepper(nu, Nx, Ny, Nz, Nthreads)

    ALLOCATE(Q(    -2:Nx+3,-2:Ny+3,-2:Nz+3, 0:Nthreads-1))
    ALLOCATE(Qend( -2:Nx+3,-2:Ny+3,-2:Nz+3, 0:Nthreads-1))
    ALLOCATE(GQ(   -2:Nx+3,-2:Ny+3,-2:Nz+3, 0:Nthreads-1))

    ! First-touch allocation of auxiliary buffers

    !$OMP PARALLEL DO schedule(static)
    DO nt=0,Nthreads-1
        Q(:,:,:,nt)    = 0.0
        Qend(:,:,:,nt) = 0.0
        GQ(:,:,:,nt)   = 0.0
        timer_fine(nt)   = 0.0
        timer_coarse(nt) = 0.0
        timer_comm(nt)   = 0.0

        ! Initialize locks, one for each timeslice
        CALL OMP_INIT_LOCK(nlocks(nt))
        
    END DO
    !$OMP END PARALLEL DO


  END SUBROUTINE InitializePararealOpenMP_Pipe

  !> Key routine to run Parareal with OpenMP and pipelining.
  !> @param[in] Q_initial Initial value
  !> @param[in] Tend Final simulation time
  !> @param[in] N_fine Number of fine steps *per time slice*
  !> @param[in] N_coarse Number of coarse steps *per time slice*
  !> @param[in] Niter Number of Parareal iterations
  !> @param[in] do_io Whether to perform IO or not
  !> @param[in] be_verbose If true, Parareal gives several status messages that can aid in debugging
  SUBROUTINE PararealOpenMP_Pipe(Q_initial, Tend, N_fine, N_coarse, Niter, dx, dy, dz, do_io, be_verbose)
    DOUBLE PRECISION, DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: Q_initial
    DOUBLE PRECISION,                         INTENT(IN)    :: Tend, dx, dy, dz
    INTEGER,                                  INTENT(IN)    :: N_fine, N_coarse, Niter
    LOGICAL,                                  INTENT(IN)    :: do_io, be_verbose

    ! Divide time interval [0,T] into Nproc many timeslices
    dt_slice  = Tend/DBLE(Nthreads)
    dt_fine   = dt_slice/DBLE(N_fine)
    dt_coarse = dt_slice/DBLE(N_coarse)

    ! Compute timeslices
    DO nt=0,Nthreads-1
        tstart_myslice(nt) = DBLE(nt)*dt_slice
        tend_myslice(nt)   = DBLE(nt+1)*dt_slice
    END DO

    IF (be_verbose) THEN
        WRITE(*,'(A, I2)') '--- Running OpenMP-pipe parareal, no. of threads: ', Nthreads

        WRITE(*,'(A, F9.5)') 'Time slice length:  ', dt_slice
        WRITE(*,'(A, F9.5)') 'Fine step length:   ', dt_fine
        WRITE(*,'(A, F9.5)') 'Coarse step length: ', dt_coarse

        DO nt=0,Nthreads-1
            WRITE(*,'(A, I2, A, F5.3, A, F5.3)') 'Thread ', nt, ' from ', tstart_myslice(nt), ' to ', tend_myslice(nt)
        END DO

    END IF

    ! --- START PARAREAL --- !

    timer_all = OMP_GET_WTIME()

    ! Set initial value:
    ! Q(0) <- y^0_0
    Q(:,:,:,0) = Q_initial

    !$OMP PARALLEL private(thread_nr, k, T0, T1)

    thread_nr = omp_get_thread_num()

    T0 = OMP_GET_WTIME()

    !$OMP DO schedule(static)
    DO nt=0,Nthreads-1

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

    END DO
    !$OMP END DO NOWAIT

    timer_coarse(thread_nr) = timer_coarse(thread_nr) + OMP_GET_WTIME() - T0

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
                Q(:,:,:,0) = Q_initial
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

            T0 = OMP_GET_WTIME()

            IF (nt<Nthreads-1) THEN
                ! Q(nt+1) = Qend(nt) = y^k_(nt+1)

                ! Update value in Q(nt+1), set lock to avoid starting updating it while
                ! thread nt+1 is still busy writing the update fine value in it (in call to F above)
                CALL OMP_SET_LOCK(nlocks(nt+1))
                Q(:,:,:,nt+1) = Qend(:,:,:,nt)
                CALL OMP_UNSET_LOCK(nlocks(nt+1))
            
            END IF

            timer_comm(thread_nr) = timer_comm(thread_nr) + OMP_GET_WTIME() - T0

            !$OMP END ORDERED

        END DO
        !$OMP END DO NOWAIT

    END DO
    !$OMP END PARALLEL

    ! Return final value in Q_initial
    Q_initial = Qend(:,:,:,Nthreads-1)

    timer_all = OMP_GET_WTIME() - timer_all

    IF(do_io) THEN
        DO nt=0,Nthreads-1
            WRITE(filename, '(A,I0.2,A,I0.2,A,I0.2,A)') 'q_final_', Niter, '_', nt, '_', Nthreads, '_openmp_pipe.dat'
            OPEN(unit=nt, FILE=filename, ACTION='write', STATUS='replace')
            WRITE(nt, '(F35.25)') Qend(1:param%Nx, 1:param%Ny, 1:param%Nz, nt)
            CLOSE(nt)
        END DO
    END IF

    DO nt=0,Nthreads-1
        WRITE(filename, '(A,I0.2,A,I0.2,A,I0.2,A)') 'timings_openmp_pipe', Niter, '_', nt, '_', Nthreads, '.dat'
        OPEN(unit=nt, FILE=filename, ACTION='write', STATUS='replace')
        WRITE(nt, '(F8.2)') timer_all
        WRITE(nt, '(F8.2)') timer_fine(nt)
        WRITE(nt, '(F8.2)') timer_coarse(nt)
        WRITE(nt, '(F8.2)') timer_comm(nt)
        CLOSE(nt)
    END DO

  END SUBROUTINE PararealOpenMP_Pipe

  !> Finalize module and deallocate buffers
  SUBROUTINE FinalizePararealOpenMP_Pipe()
      CALL FinalizeTimestepper()
      DEALLOCATE(Q)
      DEALLOCATE(Qend)
      DEALLOCATE(GQ)
      DEALLOCATE(timer_fine)
      DEALLOCATE(Timer_coarse)
      DEALLOCATE(timer_comm)
      DEALLOCATE(tstart_myslice)
      DEALLOCATE(tend_myslice)
      DO nt=0,Nthreads-1
        CALL OMP_DESTROY_LOCK(nlocks(nt))
      END DO
      DEALLOCATE(nlocks)
  END SUBROUTINE FinalizePararealOpenMP_Pipe

END MODULE parareal_openmp_pipe