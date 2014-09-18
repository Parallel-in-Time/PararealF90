PROGRAM timestepper_test

USE timestepper, only : Euler, Rk3Ssp, InitializeTimestepper, FinalizeTimestepper

IMPLICIT NONE

INCLUDE 'mpif.h'

DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932_8 ! underscore indicates rounding to real(8) precision

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, Qref
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: err_v, dt_v, convrate

DOUBLE PRECISION :: dt, dx, dy, dz, x, y, z, T0, T1, Tend
INTEGER :: Nx, Ny, Nz, Nthreads, i, j, k, Nsteps, order_adv, nt, &
    mpi_thread_provided, ierr, method, kk, oodt, order_diff
    
INTEGER, DIMENSION(7) :: N_v = (/ 260, 270, 280, 290, 300, 310, 620 /)

CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)

Nthreads = 1
Tend     = 0.2

DO method = 1,3,2

    IF (method==1) THEN
        order_adv  = 1
        order_diff = 2
    ELSE
        order_adv  = 5
        order_diff = 4
    END IF    

       
    ALLOCATE(err_v(1:SIZE(N_v)-1))
    ALLOCATE(dt_v(1:SIZE(N_v)-1))
    ALLOCATE(convrate(1:SIZE(N_v)-2))
              
    Nx = 50
    Ny = 51
    Nz = 52
  
    dy = 1.0/DBLE(Ny)
    dz = 1.0/DBLE(Nz)
    dx = 1.0/DBLE(Nx)

    ALLOCATE(Q(   -2:Nx+3, -2:Ny+3, -2:Nz+3, 0:Nthreads-1))
    ALLOCATE(Qref(-2:Nx+3, -2:Ny+3, -2:Nz+3, 0:Nthreads-1))

    CALL InitializeTimestepper(0.5*dz, Nx, Ny, Nz, Nthreads)
                        
    DO kk=SIZE(N_v),1,-1

        oodt = N_v(kk)

        dt     = 1.0/DBLE(oodt)
        Nsteps = NINT(Tend/dt)
        dt     = Tend/DBLE(Nsteps)

        ! Reset Q
        DO k=1,Nz
            DO j=1,Ny
                DO i=1,Nx

                    x = 0.5*dx + DBLE(i - 1)*dx
                    y = 0.5*dy + DBLE(j - 1)*dy
                    z = 0.5*dz + DBLE(k - 1)*dz

                    Q(i,j,k,:) = SIN(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z)

                END DO
            END DO
        END DO

        T0 = MPI_WTIME()
        !$OMP PARALLEL DO schedule(static)
        DO nt=0,Nthreads-1       
            IF (method==1) THEN
                CALL Euler(Q(:,:,:,nt), DBLE(0.0), Tend, Nsteps, dx, dy, dz, order_adv, order_diff)
        
            ELSE IF (method==3) THEN
                CALL Rk3Ssp(Q(:,:,:,nt), DBLE(0.0), Tend, Nsteps, dx, dy, dz, order_adv, order_diff)
            ELSE
                WRITE(*,*) 'Method can only have values 1 and 3. Now exiting...'
                STOP
            END IF
        END DO
        !$OMP END PARALLEL DO
        T1 = MPI_WTIME()

        DO nt=1,Nthreads-1
            IF (MAXVAL(ABS(Q(:,:,:,nt)-Q(:,:,:,nt-1)))>1e-14) THEN
                WRITE(*,'(A, I1, A)') 'For method=', method, ' not all threads computed the same result, this should not have happened. Now exiting...'
                STOP
            END IF
        END DO

        IF (ANY(ISNAN(Q))) THEN
            WRITE(*,'(A, I1, A)') 'For method=', method, ' found NAN in solution, run is probably unstable. Now exiting...'
            STOP
        END IF
        ! For nonlinear advection, no analytic solution is available to use solution
        ! with finest time-step as reference
        IF (kk==SIZE(N_v)) THEN
            Qref = Q
        ELSE            
            dt_v(kk)  = dt
            err_v(kk) = MAXVAL(ABS(Q(1:Nx,1:Ny,1:Nz,:) - Qref(1:Nx,1:Ny,1:Nz,:)))
            !WRITE(*,'(E9.3)') err_v(kk)
        END IF
    

    END DO ! N_v

    DO kk=1,SIZE(N_v)-2
        convrate(kk) = LOG10( err_v(kk+1)/err_v(kk) )/LOG10( dt_v(kk+1)/dt_v(kk) )
    END DO

    IF (MINVAL(convrate)<=0.95*method) THEN
        WRITE(*,'(A, I1, A, I1)') 'Failed to verify convergence order for method = ', method, ', order = ', order_adv
        DO kk=1,SIZE(N_v)-2
           WRITE(*,'(F6.3)') convrate(kk)
        END DO
        STOP
    END IF

    !DO kk=1,SIZE(N_v)-2
    !   WRITE(*,'(F6.3)') convrate(kk)
    !END DO
        
    CALL FinalizeTimestepper

    DEALLOCATE(Q)
    DEALLOCATE(Qref)    
    DEALLOCATE(err_v)
    DEALLOCATE(dt_v)
    DEALLOCATE(convrate)
                   
END DO ! method

CALL MPI_FINALIZE(ierr)

WRITE(*,*) '**** Timestepper: Test successful ****'

END PROGRAM timestepper_test