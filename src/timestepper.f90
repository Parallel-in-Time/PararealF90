!>
!! @todo docu
!!
!! \\( q_t = F(q) \\)
!!
MODULE timestepper

USE omp_lib,   only : omp_get_thread_num
USE Advection, only : GetRHSAdvection, InitializeAdvection, FinalizeAdvection
USE Diffusion, only : GetRHSDiffusion, InitializeDiffusion
USE boundaries, only : periodic, InitializeBoundaries, FinalizeBoundaries

IMPLICIT NONE

PRIVATE
PUBLIC :: Euler, Rk3Ssp, InitializeTimestepper, FinalizeTimestepper


!> @todo docu
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q_aux, FQ

!> @todo docu
TYPE timestepper_parameter
    INTEGER :: i_max, j_max, k_max
END TYPE

!> @todo docu
TYPE(timestepper_parameter) :: param

CONTAINS

    !> @todo docu
    SUBROUTINE Euler(Q, T0, T1, Nsteps, dx, dy, dz, order_adv, order_diff)
    
        DOUBLE PRECISION, DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: Q
        DOUBLE PRECISION,                         INTENT(IN)    :: dx, dy, dz, T0, T1
        INTEGER,                                  INTENT(IN)    :: Nsteps, order_adv, order_diff
        
        DOUBLE PRECISION :: dt
        INTEGER :: n, thread_nr
        
        thread_nr = omp_get_thread_num();
        dt = (T1-T0)/DBLE(Nsteps)
        
        DO n=1,Nsteps

            CALL periodic(Q)           

            CALL GetRHSAdvection(Q, FQ(:,:,:,thread_nr), dx, dy, dz, 1, param%i_max, 1, param%j_max, 1, param%k_max, order_adv)
            CALL GetRHSDiffusion(Q, FQ(:,:,:,thread_nr), dx, dy, dz, 1, param%i_max, 1, param%j_max, 1, param%k_max, order_diff)
            Q = Q + dt*FQ(:,:,:,thread_nr)
            
        END DO

    END SUBROUTINE Euler
    
    !> Third order strong stability preserving Runge-Kutta scheme from Shu and Osher (1988),
    !! see e.g. Durran, "Numerical Methods for Fluid Dynamics", pp. 55f
    !! @todo add docu
    SUBROUTINE Rk3Ssp(Q, T0, T1, Nsteps, dx, dy, dz, order_adv, order_diff)
    
        DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(INOUT) :: Q
        DOUBLE PRECISION,                   INTENT(IN)    :: dx, dy, dz, T0, T1
        INTEGER,                            INTENT(IN)    :: Nsteps, order_adv, order_diff
            
        DOUBLE PRECISION :: dt    
        INTEGER :: i, thread_nr
        
        thread_nr = omp_get_thread_num();
        dt = (T1-T0)/DBLE(Nsteps)
        
        DO i=1,Nsteps

                CALL periodic(Q)
        
                ! Q_1 = Q_n + dt*f(Q_n)
                CALL GetRHSAdvection(Q, FQ(:,:,:,thread_nr), dx, dy, dz, 1, param%i_max, 1, param%j_max, 1, param%k_max, order_adv)
                CALL GetRHSDiffusion(Q, FQ(:,:,:,thread_nr), dx, dy, dz, 1, param%i_max, 1, param%j_max, 1, param%k_max, order_diff)
                Q_aux(:,:,:,thread_nr) = Q + dt*FQ(:,:,:,thread_nr)

                CALL periodic(Q_aux(:,:,:,thread_nr))
                ! Q_2 = (3/4)*Q_n + (1/4)*( Q_1 + dt*f(Q_1) )
                CALL GetRHSAdvection(Q_aux(:,:,:,thread_nr), FQ(:,:,:,thread_nr), dx, dy, dz, 1, param%i_max, 1, param%j_max, 1, param%k_max, order_adv)
                CALL GetRHSDiffusion(Q_aux(:,:,:,thread_nr), FQ(:,:,:,thread_nr), dx, dy, dz, 1, param%i_max, 1, param%j_max, 1, param%k_max, order_diff)               
                Q_aux(:,:,:,thread_nr) = (3.0/4.0)*Q + (1.0/4.0)*( Q_aux(:,:,:,thread_nr) + dt*FQ(:,:,:,thread_nr) )  
                
                CALL periodic(Q_aux(:,:,:,thread_nr))
                ! Q_n+1 = (1/3)*Q_n + (2/3)*( Q_2 + dt*f(Q_2) )
                CALL GetRHSAdvection(Q_aux(:,:,:,thread_nr), FQ(:,:,:,thread_nr), dx, dy, dz, 1, param%i_max, 1, param%j_max, 1, param%k_max, order_adv)
                CALL GetRHSDiffusion(Q_aux(:,:,:,thread_nr), FQ(:,:,:,thread_nr), dx, dy, dz, 1, param%i_max, 1, param%j_max, 1, param%k_max, order_diff)                             
                Q = (1.0/3.0)*Q + (2.0/3.0)*( Q_aux(:,:,:,thread_nr) + dt*FQ(:,:,:,thread_nr) )       
        END DO
            
    END SUBROUTINE Rk3Ssp

    !> @todo docu
    SUBROUTINE InitializeTimestepper(nu, i_max, j_max, k_max, nthreads)
    
        DOUBLE PRECISION, INTENT(IN) :: nu
        INTEGER, INTENT(IN) :: i_max, j_max, k_max, nthreads
        
        INTEGER :: i
        
        param%i_max = i_max
        param%j_max = j_max
        param%k_max = k_max
        
        ! Use indexing from 1,...,i_max and account for three halo cells in each direction
        CALL InitializeAdvection(-2, i_max+3, -2, j_max+3, -2, k_max+3, nthreads, .false.)
        CALL InitializeDiffusion(-2, i_max+3, -2, j_max+3, -2, k_max+3, nu)
        CALL InitializeBoundaries(i_max, j_max, k_max)
        
        ALLOCATE(Q_aux(-2:i_max+3, -2:j_max+3, -2:k_max+3,0:nthreads-1))
        ALLOCATE(FQ(   -2:i_max+3, -2:j_max+3, -2:k_max+3,0:nthreads-1))

        !$OMP PARALLEL DO schedule(static)
        DO i=0,nthreads-1
            Q_aux(:,:,:,i) = 0.0
            FQ(:,:,:,i)    = 0.0
        END DO
        !$OMP END PARALLEL DO
        
    END SUBROUTINE InitializeTimestepper

    !> @todo docu
    SUBROUTINE FinalizeTimestepper()
            CALL FinalizeAdvection
            CALL FinalizeBoundaries
            DEALLOCATE(Q_aux)
            DEALLOCATE(FQ)    
    END SUBROUTINE FinalizeTimestepper
    
END MODULE timestepper