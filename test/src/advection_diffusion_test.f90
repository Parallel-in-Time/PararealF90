!>
! Runs a linear advection diffusion equation with increasing resolution and checks that the observed
! convergence orders match the order of the involved discretizations.
! Additionally, the same test is run N times by N threads and it is checked that all threads compute the same results.
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
INTEGER, PARAMETER, DIMENSION(5) :: N_v = (/ 40, 60, 80, 100, 120 /)

DOUBLE PRECISION, PARAMETER :: Tend = 0.05

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, RQ, Qref
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: error, convrate

DOUBLE PRECISION :: dx, dy, dz, x, y, z, dt, T0, T1, nu
INTEGER :: method, i, j, k, Nx, Ny, Nz, nt, nn, order_adv, order_diff, Nsteps, Nthreads, mpi_thread_provided, ierr

CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)


DO method=3,3,2
    
        IF  (method==1) THEN
            order_adv  = 1
            order_diff = 2
        ELSE IF (method==3) THEN
            order_adv  = 5
            order_diff = 4
        END IF
        
        ALLOCATE(error(SIZE(N_v)))
        ALLOCATE(convrate(SIZE(N_v)-1))
                  
        DO nn=1,SIZE(N_v)
        
            Nx = N_v(nn)-3
            Ny = N_v(nn)+1
            Nz = N_v(nn)+3
        
            dx = 1.0/DBLE(Nx)
            dy = 1.0/DBLE(Ny)
            dz = 1.0/DBLE(Nz)
                   
            nu = 0.0025
            dt = 0.5*(dz*dz)/nu
            Nsteps = CEILING(Tend/dt)
            dt     = Tend/DBLE(Nsteps)

           DO Nthreads=1,6

                ALLOCATE(Q(   -2:Nx+3, -2:Ny+3, -2:Nz+3, 0:Nthreads-1))
                ALLOCATE(RQ(  -2:Nx+3, -2:Ny+3, -2:Nz+3, 0:Nthreads-1))
                ALLOCATE(Qref(-2:Nx+3, -2:Ny+3, -2:Nz+3, 0:Nthreads-1))

                CALL InitializeTimestepper(nu, Nx, Ny, Nz, Nthreads)
            
                ! Compute initial and analytic solution 
                DO k=1,Nz
                    DO j=1,Ny
                        DO i=1,Nx

                        x = 0.5*dx + DBLE(i - 1)*dx
                        y = 0.5*dy + DBLE(j - 1)*dy
                        z = 0.5*dz + DBLE(k - 1)*dz
                    
                        Q(i,j,k,:)    = SIN(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z)    
                        Qref(i,j,k,:) = EXP(-12*nu*pi*pi*Tend)*SIN(2.0*pi*(x-Tend))*SIN(2.0*pi*(y-Tend))*SIN(2.0*pi*(z-Tend))
                                    
                        END DO
                    END DO
                END DO
            
                T0 = MPI_WTIME()
                !$OMP PARALLEL DO schedule(static)
                DO nt=0,Nthreads-1
            
                    IF (method==1) THEN
                        CALL Euler( Q(:,:,:,nt), DBLE(0.0), Tend, Nsteps, dx, dy, dz, order_adv, order_diff)
                    ELSE IF (method==3) THEN
                        CALL Rk3Ssp(Q(:,:,:,nt), DBLE(0.0), Tend, Nsteps, dx, dy, dz, order_adv, order_diff)
                    END IF
            
                END DO
                !$OMP END PARALLEL DO
                T1 = MPI_WTIME()
                            
                DO nt=1,Nthreads-1
                     IF (MAXVAL(ABS(Q(:,:,:,nt)-Q(:,:,:,nt-1)))>1e-14) THEN
                        WRITE(*,'(A, I1, A)') 'ERROR: For method=', method, ' not all threads computed the same result, this should not have happened. Now exiting...'
                        STOP
                    END IF               
                END DO

                error(nn) = MAXVAL(ABS(Q(1:Nx,1:Ny,1:Nz,0) - Qref(1:Nx,1:Ny,1:Nz,0)))/MAXVAL(ABS(Qref(1:Nx,1:Ny,1:Nz,0)))
                
                WRITE(*,'(E9.3)') error(nn)
                
                CALL FinalizeTimestepper  
                      
                DEALLOCATE(Q)
                DEALLOCATE(RQ)
                DEALLOCATE(Qref)
            
            END DO ! Nthreads
        END DO ! N_v
        
        DO nn=1,SIZE(N_v)-1
            convrate(nn) = LOG10( error(nn+1)/error(nn) )/LOG10( DBLE(N_v(nn))/DBLE(N_v(nn+1)) )
            WRITE(*,'(F9.5)') convrate(nn)
        END DO
                
        DEALLOCATE(error)
        DEALLOCATE(convrate)
        
END DO ! method

CALL MPI_FINALIZE(ierr)

#endif

END PROGRAM advection_diffusion_test