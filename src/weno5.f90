MODULE weno5
! This module provides the implementation of a two-dimensional, finite difference WENO-5 scheme.
!
! Daniel Ruprecht, 8.2.2012
! ICS Lugano

USE omp_lib, only : omp_get_thread_num
USE fluxes,  only : FluxI, FluxJ, FluxK, param, FluxI_Cell, FluxJ_Cell, FluxK_Cell

IMPLICIT NONE

#include <preprocessor.f90>

PRIVATE
PUBLIC :: WenoFluxes

! Define fixed parameters used by the WENO-5 method (see Shu .... e.g.)
DOUBLE PRECISION, PARAMETER, DIMENSION(3) :: weights_plus = (/ 0.3_8, 0.6_8, 0.1_8 /)
DOUBLE PRECISION, PARAMETER, DIMENSION(5) :: stencil_weights = (/ 2.0_8, 5.0_8, -1.0_8, -7.0_8, 11.0_8 /)*(1.0_8/6.0_8)
DOUBLE PRECISION, PARAMETER               :: coeff_1 = 13.0_8/12.0_8, coeff_2 = 1.0_8/4.0_8
DOUBLE PRECISION, PARAMETER               :: weno_tol = 1.0e-6_8
DOUBLE PRECISION, PARAMETER               :: weno_n   = 2.0_8

CONTAINS

    SUBROUTINE WenoFluxes(Q)
    
        DOUBLE PRECISION, DIMENSION(param%i_min:, param%j_min:, param%k_min:), INTENT(IN)  :: Q

        INTEGER :: thread_nr

        thread_nr = omp_get_thread_num()

#if(linear==1)
        ! Evaluate flux function (here f(q) = q) to get cell values of flux
        FluxI_Cell(:,:,:,thread_nr) = Q
        FluxJ_Cell(:,:,:,thread_nr) = Q
        FluxK_Cell(:,:,:,thread_nr) = Q
#elif(linear==0)     
        ! Nonlinear advection: f(q)  = 0.5*q**2
        FluxI_Cell(:,:,:,thread_nr) = 0.5*Q*Q
        FluxJ_Cell(:,:,:,thread_nr) = 0.5*Q*Q
        FluxK_Cell(:,:,:,thread_nr) = 0.5*Q*Q
#else
        WRITE(*,*) 'Found value other than 0 or 1 for preprocessor flag <linear>. Now exiting...'
        STOP        
#endif            
        ! Now update interface values of horizontal and vertical flux
        CALL UpdateFluxI(Q, 1.0_8 )
        CALL UpdateFluxJ(Q, 1.0_8 )
        CALL UpdateFluxK(Q, 1.0_8 )
            
    END SUBROUTINE WenoFluxes
        
    ! Updates the horizontal fluxes using WENO5
    SUBROUTINE UpdateFluxI(Qcell, max_vel)
        
        DOUBLE PRECISION, DIMENSION(param%i_min:,param%j_min:,param%k_min:), INTENT(IN) :: Qcell
        DOUBLE PRECISION,                                                                      INTENT(IN) :: max_vel
        
        DOUBLE PRECISION, DIMENSION(6) :: Qcell_local, Fcell_local
        INTEGER :: i, j, k, thread_nr
        
        thread_nr = omp_get_thread_num()
        
        ! Out of the global fields Qcell and FluxQcell, updated interface
        ! values of the flux are computed

        DO k=param%k_start,param%k_end
            DO j=param%j_start,param%j_end
                
                ! Left boundary
                i=param%i_start-1
                
                ! Initialize "moving window" arrays
                Qcell_local(1) = Qcell(i-2,j,k)
                Qcell_local(2) = Qcell(i-1,j,k)
                Qcell_local(3) = Qcell(i+0,j,k)
                Qcell_local(4) = Qcell(i+1,j,k)
                Qcell_local(5) = Qcell(i+2,j,k)
                Qcell_local(6) = Qcell(i+3,j,k)

                Fcell_local(1) = FluxI_Cell(i-2,j,k,thread_nr)
                Fcell_local(2) = FluxI_Cell(i-1,j,k,thread_nr)
                Fcell_local(3) = FluxI_Cell(i+0,j,k,thread_nr)
                Fcell_local(4) = FluxI_Cell(i+1,j,k,thread_nr)
                Fcell_local(5) = FluxI_Cell(i+2,j,k,thread_nr)
                Fcell_local(6) = FluxI_Cell(i+3,j,k,thread_nr)
                
                ! Reconstruct interface value from local stencil
                FluxI(i+1,j,k,thread_nr) = ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel)

                ! Inner cells
                DO i=param%i_start,param%i_end
                        ! Update local "moving window" arrays
                        CALL ShiftLocalValues(Qcell_local, Fcell_local, Qcell(i+3,j,k), FluxI_Cell(i+3,j,k,thread_nr))
                        ! Reconstruct flux value at interface from current local window
                        FluxI(i+1,j,k,thread_nr) =  ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel)
                END DO
            
            END DO
        END DO

    END SUBROUTINE UpdateFluxI
    
    ! Update vertical fluxes using WENO5
    SUBROUTINE UpdateFluxJ(Qcell, max_vel)
    
        DOUBLE PRECISION, DIMENSION(param%i_min:,param%j_min:,param%k_min:), INTENT(IN) :: Qcell
        DOUBLE PRECISION,                                                                      INTENT(IN) :: max_vel
        
        DOUBLE PRECISION, DIMENSION(6) :: Qcell_local, Fcell_local
        INTEGER :: i, j, k, thread_nr
        
        thread_nr = omp_get_thread_num()
        
        ! Out of the global fields Qcell and FluxQcell, updated interface
        ! values of the flux are computed

        DO k=param%k_start, param%k_end
            DO i=param%i_start,param%i_end
            
                j=param%j_start-1
                Qcell_local(1) = Qcell(i,j-2,k)
                Qcell_local(2) = Qcell(i,j-1,k)
                Qcell_local(3) = Qcell(i,j+0,k)
                Qcell_local(4) = Qcell(i,j+1,k)
                Qcell_local(5) = Qcell(i,j+2,k)
                Qcell_local(6) = Qcell(i,j+3,k)

                Fcell_local(1) = FluxJ_Cell(i,j-2,k,thread_nr)
                Fcell_local(2) = FluxJ_Cell(i,j-1,k,thread_nr)
                Fcell_local(3) = FluxJ_Cell(i,j+0,k,thread_nr)
                Fcell_local(4) = FluxJ_Cell(i,j+1,k,thread_nr)
                Fcell_local(5) = FluxJ_Cell(i,j+2,k,thread_nr)
                Fcell_local(6) = FluxJ_Cell(i,j+3,k,thread_nr)
                
                FluxJ(i,j+1,k,thread_nr) = ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel)

                DO j=param%j_start,param%j_end
                    CALL ShiftLocalValues(Qcell_local, Fcell_local, Qcell(i,j+3,k), FluxJ_Cell(i,j+3,k,thread_nr))
                    FluxJ(i,j+1,k,thread_nr) = ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel)
                END DO
                               
            END DO
        END DO
    
    END SUBROUTINE UpdateFluxJ

    SUBROUTINE UpdateFluxK(Qcell, max_vel)

        DOUBLE PRECISION, DIMENSION(param%i_min:,param%j_min:,param%k_min:), INTENT(IN) :: Qcell
        DOUBLE PRECISION,                                                              INTENT(IN) :: max_vel
        
        DOUBLE PRECISION, DIMENSION(6) :: Qcell_local, Fcell_local
        INTEGER :: i, j, k, thread_nr
        
        thread_nr = omp_get_thread_num()
        
        ! Out of the global fields Qcell and FluxQcell, updated interface
        ! values of the flux are computed


        DO i=param%i_start,param%i_end
            DO j=param%j_start,param%j_end

                k=param%k_start-1
                Qcell_local(1) = Qcell(i,j,k-2)
                Qcell_local(2) = Qcell(i,j,k-1)
                Qcell_local(3) = Qcell(i,j,k+0)
                Qcell_local(4) = Qcell(i,j,k+1)
                Qcell_local(5) = Qcell(i,j,k+2)
                Qcell_local(6) = Qcell(i,j,k+3)

                Fcell_local(1) = FluxK_Cell(i,j,k-2,thread_nr)
                Fcell_local(2) = FluxK_Cell(i,j,k-1,thread_nr)
                Fcell_local(3) = FluxK_Cell(i,j,k+0,thread_nr)
                Fcell_local(4) = FluxK_Cell(i,j,k+1,thread_nr)
                Fcell_local(5) = FluxK_Cell(i,j,k+2,thread_nr)
                Fcell_local(6) = FluxK_Cell(i,j,k+3,thread_nr)
                
                FluxK(i,j,k+1,thread_nr) = ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel)

                DO k=param%k_start,param%k_end
                    CALL ShiftLocalValues(Qcell_local, Fcell_local, Qcell(i,j,k+3), FluxK_Cell(i,j,k+3,thread_nr))
                    FluxK(i,j,k+1,thread_nr) = ReconstructInterfaceValue(Qcell_local, Fcell_local, max_vel)
                END DO

            END DO
        END DO

    END SUBROUTINE UpdateFluxK

! 
! Below are three auxiliary functions used by the WENO method that all operate on the
! local solution values: At each point, WENO5 uses information from a six-point wide 
! stencil. The buffer Qcell_local contains the cell-values of the solution Q at these six
! points, Fcell_local contains the cell-values of the flux at these six points. Using
! the WENO procedure, from these 2x6 values a new value for Flux(I,J,K)_Int, that is
! a flux value at the interface is reconstructed.
!
    
    ! Moves all values in Qcell_local and Fcell_local one index down and inserts new
    ! values at index 6.
    PURE SUBROUTINE ShiftLocalValues(Qcell_local, Fcell_local, Qcell_new, Fcell_new)

    DOUBLE PRECISION, DIMENSION(6), INTENT(INOUT) :: Qcell_local, Fcell_local
    DOUBLE PRECISION,               INTENT(IN)    :: Qcell_new, Fcell_new
    
    INTEGER :: i
                    
    ! Shift all values one to the left              
    DO i=1,5
        Qcell_local(i) = Qcell_local(i+1)
        Fcell_local(i) = Fcell_local(i+1)
    END DO
    
    ! ..and fill in the new values to the right
    Qcell_local(6) = Qcell_new
    Fcell_local(6) = Fcell_new              
                
    END SUBROUTINE ShiftLocalValues
    
    ! This routine performs the actual reconstruction of the interface values    
    FUNCTION ReconstructInterfaceValue(Qcell_local, Fcell_local, local_vel) RESULT(Fint)
    
        DOUBLE PRECISION, DIMENSION(6) :: Qcell_local, Fcell_local
        DOUBLE PRECISION               :: local_vel, Fint
        
        DOUBLE PRECISION, DIMENSION(5) :: Fcell_plus, Fcell_minus
        DOUBLE PRECISION               :: Fint_plus,  Fint_minus
        
        ! From the local values of Q and flux of Q, retrieve the positive and negative
        ! components of the flux function Fcell_plus and Fcell_minus
        CALL GetLocalFluxSplit(Qcell_local, Fcell_local, local_vel, Fcell_plus, Fcell_minus)
    
        ! Perform a WENO reconstruction of the interface value for the positive and negative
        ! component
        CALL GetWenoInterfaceValue(Fcell_plus, Fcell_minus, Fint_plus, Fint_minus)
        
        ! The final interface value is the sum of the values reconstructed from the positive
        ! and negative flux component
        Fint = Fint_plus + Fint_minus
        
        CONTAINS
        
            ! Local Lax-Friedrichs flux-splitting
            PURE SUBROUTINE GetLocalFluxSplit(Qcell_local, Fcell_local, local_vel, Fcell_plus, Fcell_minus)
    
                DOUBLE PRECISION, DIMENSION(6), INTENT(IN)  :: Qcell_local, Fcell_local
                DOUBLE PRECISION,               INTENT(IN)  :: local_vel
                DOUBLE PRECISION, DIMENSION(5), INTENT(OUT) :: Fcell_plus, Fcell_minus
            
                ! Lax-Friedrichs flux-splitting
                Fcell_plus( 1:5) = 0.5_8*( Fcell_local(1:5) + local_vel*Qcell_local(1:5) )
                Fcell_minus(1:5) = 0.5_8*( Fcell_local(2:6) - local_vel*Qcell_local(2:6) )
                                        
            END SUBROUTINE GetLocalFluxSplit
    
            ! In this routine, the WENO magic happens: The three candidate stencils are
            ! calculated, the smoothness measures beta and then the final weights
            PURE SUBROUTINE GetWenoInterfaceValue(Fcell_plus, Fcell_minus, Fint_plus, Fint_minus)
        
                DOUBLE PRECISION, DIMENSION(5), INTENT(IN)  :: Fcell_plus, Fcell_minus
                DOUBLE PRECISION,               INTENT(OUT) :: Fint_plus, Fint_minus
                
                DOUBLE PRECISION, DIMENSION(3) :: pol_values, beta, alpha
                DOUBLE PRECISION               :: alpha_sum_inv
                
                ! Compute three possible stencils with indices 1,2,3 // 2,3,4 // 3,4,5
                pol_values(1) = stencil_weights(1)*Fcell_plus(3) + stencil_weights(2)*Fcell_plus(4) + stencil_weights(3)*Fcell_plus(5)
                pol_values(2) = stencil_weights(3)*Fcell_plus(2) + stencil_weights(2)*Fcell_plus(3) + stencil_weights(1)*Fcell_plus(4)
                pol_values(3) = stencil_weights(1)*Fcell_plus(1) + stencil_weights(4)*Fcell_plus(2) + stencil_weights(5)*Fcell_plus(3) 
                
                ! Second, compute smoothness measures (magic formula from WENO paper)
                beta(1) = coeff_1*( Fcell_plus(3) - 2.0_8*Fcell_plus(4) + Fcell_plus(5) )**2.0 + coeff_2*( 3.0_8*Fcell_plus(3) - 4.0_8*Fcell_plus(4) +     Fcell_plus(5) )**2.0
                beta(2) = coeff_1*( Fcell_plus(2) - 2.0_8*Fcell_plus(3) + Fcell_plus(4) )**2.0 + coeff_2*(       Fcell_plus(2)                       -     Fcell_plus(4) )**2.0
                beta(3) = coeff_1*( Fcell_plus(1) - 2.0_8*Fcell_plus(2) + Fcell_plus(3) )**2.0 + coeff_2*(       Fcell_plus(1) - 4.0_8*Fcell_plus(2) + 3.0_8*Fcell_plus(3) )**2.0
                                
                ! Third, compute weights out of the smoothness measures
                alpha(1) = weights_plus(1)/( beta(1) + weno_tol )**weno_n   
                alpha(2) = weights_plus(2)/( beta(2) + weno_tol )**weno_n
                alpha(3) = weights_plus(3)/( beta(3) + weno_tol )**weno_n
                                    
                ! Fourth, normalize the weights.                    
                alpha_sum_inv = 1.0_8/SUM(alpha)
                alpha         = alpha_sum_inv*alpha

                ! Finally, compute the superposition of the three candidate-stencil values using the computed weights
                Fint_plus = alpha(1)*pol_values(1) + alpha(2)*pol_values(2) + alpha(3)*pol_values(3)
                
                ! *** Now perform corresponding computation for Fcell_minus ****
                pol_values(1) = stencil_weights(5)*Fcell_minus(3) + stencil_weights(4)*Fcell_minus(4) + stencil_weights(1)*Fcell_minus(5) 
                pol_values(2) = stencil_weights(1)*Fcell_minus(2) + stencil_weights(2)*Fcell_minus(3) + stencil_weights(3)*Fcell_minus(4)
                pol_values(3) = stencil_weights(3)*Fcell_minus(1) + stencil_weights(2)*Fcell_minus(2) + stencil_weights(1)*Fcell_minus(3)
                
                ! Smoothness measures
                beta(1) = coeff_1*(       Fcell_minus(3) - 2.0*Fcell_minus(4) +     Fcell_minus(5) )**2.0 + coeff_2*(   3.0_8*Fcell_minus(3) - 4.0_8*Fcell_minus(4) +       Fcell_minus(5) )**2.0
                beta(2) = coeff_1*(       Fcell_minus(2) - 2.0*Fcell_minus(3) +     Fcell_minus(4) )**2.0 + coeff_2*(         Fcell_minus(2)                        -       Fcell_minus(4) )**2.0
                beta(3) = coeff_1*(       Fcell_minus(1) - 2.0*Fcell_minus(2) +     Fcell_minus(3) )**2.0 + coeff_2*(         Fcell_minus(1) - 4.0_8*Fcell_minus(2) + 3.0_8*Fcell_minus(3) )**2.0
                    
                ! Compute weights   
                alpha(1) = weights_plus(3)/( beta(1) + weno_tol )**weno_n
                alpha(2) = weights_plus(2)/( beta(2) + weno_tol )**weno_n
                alpha(3) = weights_plus(1)/( beta(3) + weno_tol )**weno_n

                ! Normalize weights
                alpha_sum_inv = 1.0_8/SUM(alpha)
                alpha = alpha_sum_inv*alpha                                 
                
                Fint_minus = alpha(1)*pol_values(1) + alpha(2)*pol_values(2) + alpha(3)*pol_values(3)
            
            END SUBROUTINE GetWenoInterfaceValue

    END FUNCTION ReconstructInterfaceValue
 
END MODULE WENO5
