!>
!! @todo docu
!!
MODULE parareal_openmp

USE omp_lib
USE timestepper, only : Euler, Rk3Ssp, InitializeTimestepper, FinalizeTimestepper

IMPLICIT NONE

PRIVATE
PUBLIC :: InitializePararealOpenMP, FinalizePararealOpenMP, PararealOpenMP

!> @todo docu
INTEGER, PARAMETER :: order_adv_c = 1, order_diff_c = 2, order_adv_f = 5, order_diff_f = 4

!> @todo docu
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, GQ, Qend

!> @todo docu
DOUBLE PRECISION :: dt_coarse

!> @todo docu
DOUBLE PRECISION :: dt_fine

!> @todo docu
DOUBLE PRECISION :: dt_slice

!> @todo docu
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: tstart_myslice, tend_myslice, timer_fine

! Internal variables to be used for timers
DOUBLE PRECISION :: timer_all, T0, T1, timer_coarse, timer_comm

INTEGER :: k, Nthreads, nt
     
LOGICAL :: do_io, be_verbose

CHARACTER(len=64) :: filename

TYPE parareal_parameter
  INTEGER Nx, Ny, Nz
END TYPE

TYPE(parareal_parameter) :: param

CONTAINS

  !> @todo docu
  !!
  SUBROUTINE InitializePararealOpenMP(nu, Nx, Ny, Nz)
    DOUBLE PRECISION, INTENT(IN) :: nu
    INTEGER, INTENT(IN) :: Nx, Ny, Nz

    Nthreads = OMP_GET_MAX_THREADS()

    param%Nx = Nx
    param%Ny = Ny
    param%Nz = Nz

    ALLOCATE(Q(    -2:Nx+3,-2:Ny+3,-2:Nz+3, 0:Nthreads-1))
    ALLOCATE(Qend( -2:Nx+3,-2:Ny+3,-2:Nz+3, 0:Nthreads-1))
    ALLOCATE(GQ(   -2:Nx+3,-2:Ny+3,-2:Nz+3, 0:Nthreads-1))

    ALLOCATE(timer_fine(0:Nthreads-1))
    ALLOCATE(tstart_myslice(0:Nthreads-1))
    ALLOCATE(tend_myslice(0:Nthreads-1))

    CALL InitializeTimestepper(nu, Nx, Ny, Nz, Nthreads)

    ! First-touch allocation of auxiliary buffers; computing time slices

    !$OMP PARALLEL DO schedule(static)
    DO nt=0,Nthreads-1
      Q(:,:,:,nt)      = 0.0
      Qend(:,:,:,nt)   = 0.0
      GQ(:,:,:,nt)     = 0.0
      timer_fine(nt)   = 0.0
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE InitializePararealOpenMP


  !>
  !! @todo docu
  SUBROUTINE PararealOpenMP(Q_initial, Tend, N_fine, N_coarse, Niter, dx, dy, dz, do_io, be_verbose)

    DOUBLE PRECISION, DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: Q_initial
    DOUBLE PRECISION,                         INTENT(IN)    :: Tend, dx, dy, dz
    INTEGER,                                  INTENT(IN)    :: N_fine, N_coarse, Niter
    LOGICAL,                                  INTENT(IN)    :: do_io, be_verbose

    ! Divide time interval [0,Tend] into Nthreads many timeslices
    dt_slice  = Tend/DBLE(Nthreads)
    dt_fine   = dt_slice/DBLE(N_fine)
    dt_coarse = dt_slice/DBLE(N_coarse)

    DO nt=0,Nthreads-1
      tstart_myslice(nt) = DBLE(nt)*dt_slice
      tend_myslice(nt)   = DBLE(nt+1)*dt_slice
    END DO

    IF (be_verbose) THEN

        WRITE(*,'(A, I2)') '--- Running OpenMP parareal, no. of threads: ', Nthreads

        WRITE(*,'(A, F9.5)') 'Time slice length:  ', dt_slice
        WRITE(*,'(A, F9.5)') 'Fine step length:   ', dt_fine
        WRITE(*,'(A, F9.5)') 'Coarse step length: ', dt_coarse

        DO nt=0,Nthreads-1
          WRITE(*,'(A, I2, A, F5.3, A, F5.3)') 'Thread ', nt, ' from ', tstart_myslice(nt), ' to ', tend_myslice(nt)
        END DO

    END IF

    ! --- START PARAREAL --- !

    timer_all    = OMP_GET_WTIME()
    timer_coarse = 0.0
    timer_comm   = 0.0

    ! Set initial value:
    ! Q(0) <- y^0_0
    Q(:,:,:,0) = Q_initial

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

            ! Compute difference fine minus coarse
            ! Qend(nt) <- F(y^(k-1)_nt) - G(y^(k-1)_nt)
            Qend(:,:,:,nt) = Q(:,:,:,nt) - GQ(:,:,:,nt)

            T1 = OMP_GET_WTIME()
            timer_fine(nt) = timer_fine(nt) + T1 - T0

        END DO
        !$OMP END PARALLEL DO
        
        ! Now read initial value into Q(0) again:
        ! Q(0) <- y^k_0
        T0 = OMP_GET_WTIME()
        Q(:,:,:,0) = Q_initial
        timer_comm = timer_comm + OMP_GET_WTIME() - T0
        
        DO nt=0,Nthreads-1

            T0 = OMP_GET_WTIME()        
            ! GQ(nt) <- y^k_nt
            GQ(:,:,:,nt) = Q(:,:,:,nt)
            
            ! GQ(nt) <- G(y^k_nt)
            CALL Euler(GQ(:,:,:,nt), tstart_myslice(nt), tend_myslice(nt), N_coarse, dx, dy, dz, order_adv_c, order_diff_c)

            ! Qend(nt) <- y^k_(nt+1) = G(y^k_nt) + F(y^(k-1)_nt) - G(y^(k-1)_nt)
            Qend(:,:,:,nt) = GQ(:,:,:,nt) + Qend(:,:,:,nt)

            timer_coarse = timer_coarse + OMP_GET_WTIME() - T0

            ! Now update initial value for next slice
            ! Q(nt+1) = Qend(nt) <- y^k_(n+1)
            IF (nt<Nthreads-1) THEN
                T0 = OMP_GET_WTIME()
                Q(:,:,:,nt+1) = Qend(:,:,:,nt)
                timer_comm = timer_comm + OMP_GET_WTIME() - T0
            END IF
        END DO
    END DO ! Enf of parareal loop

    ! --- END PARAREAL --- !

    ! Return final value in Q_initial
    Q_initial = Qend(:,:,:,Nthreads-1)

    timer_all = OMP_GET_WTIME() - timer_all

    IF(do_io) THEN
        DO nt=0,Nthreads-1
            WRITE(filename, '(A,I0.2,A,I0.2,A)') 'q_final_', nt, '_', Nthreads, '_openmp.dat'
            OPEN(unit=nt, FILE=filename, ACTION='write', STATUS='replace')
            WRITE(nt, '(F35.25)') Qend(1:param%Nx, 1:param%Ny, 1:param%Nz, nt)
            CLOSE(nt)
        END DO
    END IF


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

  END SUBROUTINE PararealOpenMP

  !> @todo docu
  !!
  SUBROUTINE FinalizePararealOpenMP()
      CALL FinalizeTimestepper()
      DEALLOCATE(Q)
      DEALLOCATE(Qend)
      DEALLOCATE(GQ)
      DEALLOCATE(timer_fine)
      DEALLOCATE(tstart_myslice)
      DEALLOCATE(tend_myslice)
  END SUBROUTINE FinalizePararealOpenMP

END MODULE parareal_openmp