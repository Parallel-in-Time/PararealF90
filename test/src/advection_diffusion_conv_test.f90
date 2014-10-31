!>
! Runs a linear advection diffusion equation with increasing resolution and checks that the observed
! convergence orders match the order of the involved discretizations.
!
PROGRAM advection_diffusion_test

USE Timestepper, only : Euler, RK3Ssp, InitializeTimestepper, FinalizeTimestepper
IMPLICIT NONE

INCLUDE 'mpif.h'

#include <preprocessor.f90>

#if(linear==0)
    WRITE(*,*) 'WARNING: advection_diffusion_test requires linear advection to have analytical solution available... recompile with flag linear=1 to run this test. Now exiting.'
    STOP
#else

DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932_8 ! underscore indicates rounding to real(8) precision
INTEGER, PARAMETER, DIMENSION(5) :: N_v = (/ 20, 30, 40, 50, 60 /)

DOUBLE PRECISION, PARAMETER :: Tend = 0.05

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, Qref
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: error

DOUBLE PRECISION :: dx, dy, dz, x, y, z, dt, nu, convrate
INTEGER :: method, i, j, k, Nx, Ny, Nz, nn, order_adv, order_diff, Nsteps, Nthreads, mpi_thread_provided, ierr

CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)


DO method=1,3,2
    
        IF  (method==1) THEN
            order_adv  = 1
            order_diff = 2
        ELSE IF (method==3) THEN
            order_adv  = 5
            order_diff = 4
        END IF
        
        ALLOCATE(error(SIZE(N_v)))

        DO nn=1,SIZE(N_v)

            Nx = N_v(nn)-3
            Ny = N_v(nn)+1
            Nz = N_v(nn)+3
        
            dx = 1.0/DBLE(Nx)
            dy = 1.0/DBLE(Ny)
            dz = 1.0/DBLE(Nz)
                   
            !nu = 0.0025
            nu = 0.0
            !dt = 0.5*(dz*dz)/nu
            dt = 0.25*dz
            Nsteps = CEILING(Tend/dt)
            dt     = Tend/DBLE(Nsteps)

           Nthreads=1

          ALLOCATE(Q(   -2:Nx+3, -2:Ny+3, -2:Nz+3, 0:Nthreads-1))
          ALLOCATE(Qref(-2:Nx+3, -2:Ny+3, -2:Nz+3, 0:Nthreads-1))

          CALL InitializeTimestepper(nu, Nx, Ny, Nz, Nthreads)
      
          ! Compute initial and analytic solution 
          DO k=-2,Nz+3
              DO j=-2,Ny+3
                  DO i=-2,Nx+3

                  x = 0.5*dx + DBLE(i - 1)*dx
                  y = 0.5*dy + DBLE(j - 1)*dy
                  z = 0.5*dz + DBLE(k - 1)*dz
              
                  Q(i,j,k,:)    = SIN(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z)    
                  Qref(i,j,k,:) = EXP(-12*nu*pi*pi*Tend)*SIN(2.0*pi*(x-Tend))*SIN(2.0*pi*(y-Tend))*SIN(2.0*pi*(z-Tend))
                              
                  END DO
              END DO
          END DO

      
          IF (method==1) THEN
              CALL Euler( Q(:,:,:,0), DBLE(0.0), Tend, Nsteps, dx, dy, dz, order_adv, order_diff)
          ELSE IF (method==3) THEN
              CALL Rk3Ssp(Q(:,:,:,0), DBLE(0.0), Tend, Nsteps, dx, dy, dz, order_adv, order_diff)
          END IF

          error(nn) = MAXVAL(ABS(Q(1:Nx,1:Ny,1:Nz,0) - Qref(1:Nx,1:Ny,1:Nz,0)))/MAXVAL(ABS(Qref(1:Nx,1:Ny,1:Nz,0)))

          CALL FinalizeTimestepper  
                
          DEALLOCATE(Q)
          DEALLOCATE(Qref)

        END DO ! N_v
        
        DO nn=1,SIZE(N_v)-1
            convrate = LOG10( error(nn+1)/error(nn) )/LOG10( DBLE(N_v(nn))/DBLE(N_v(nn+1)) )
           IF (convrate < 0.95*DBLE(method)) THEN
              WRITE(*,'(A, I2, A)') 'ERROR: Method ', method, ' failed to reproduce expected order of convergence. Now exiting.'
              STOP
           END IF
        END DO
                
        DEALLOCATE(error)

END DO ! method

CALL MPI_FINALIZE(ierr)

WRITE(*,*) '[0] -- Successful: Both Euler and Rk3Ssp produced the expected order of convergence for the linear advection-diffusion problem.'

#endif

END PROGRAM advection_diffusion_test