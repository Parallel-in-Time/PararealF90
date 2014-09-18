MODULE upwind
! This module provides the implementation of a two-dimensional, finite difference upwind scheme.
!
! Daniel Ruprecht, 17.4.2014
! ICS Lugano

USE omp_lib,       only : omp_get_thread_num
USE fluxes,        only : param, FluxI, FluxJ, FluxK

IMPLICIT NONE

#include <preprocessor.f90>

PRIVATE
PUBLIC :: UpwindFluxes

CONTAINS

    SUBROUTINE UpwindFluxes(Q)
    
        DOUBLE PRECISION, DIMENSION(param%i_min:, param%j_min:, param%k_min:), &
                                            INTENT(IN)  :: Q

        ! Now update interface values of horizontal and vertical flux
        CALL UpdateFluxI(Q)
        CALL UpdateFluxJ(Q)
        CALL UpdateFluxK(Q)

    END SUBROUTINE UpwindFluxes

    ! Updates the horizontal fluxes using first order upwinding
    SUBROUTINE UpdateFluxI(Q)
        
        DOUBLE PRECISION, DIMENSION(param%i_min:, param%j_min:, param%k_min:), &
                        INTENT(IN) :: Q
        
        DOUBLE PRECISION :: Uadv
        INTEGER          :: i, j, k, thread_nr
#if(linear==0)
        DOUBLE PRECISION :: switch
#endif     
   
        thread_nr = omp_get_thread_num()

        DO k=param%k_start, param%k_end       
            DO j=param%j_start,param%j_end
                DO i=param%i_start,param%i_end+1
                  
#if(linear==1)
                    ! Linear advection with uadv = 1
                    Uadv = 1.0_8
                    FluxI(i,j,k,thread_nr) =  0.5_8*Uadv*( Q(i-1,j,k) + Q(i,j,k) ) - 0.5_8*ABS(Uadv)*( Q(i,j,k) - Q(i-1,j,k) )                    
                        
#elif(linear==0)                        
                    ! Nonlinear advection
                    Uadv = 0.5*(Q(i-1,j,k)+Q(i,j,k))
                        
                    ! It is switch=0 if Uadv<0 and switch=1 if Uadv>0
                    switch = 0.5*(SIGN(1.0_8,Uadv)+1.0)

                    FluxI(i,j,k,thread_nr) = (1.0-switch)*0.5*Q(i,j,k)*Q(i,j,k) + switch*0.5*Q(i-1,j,k)*Q(i-1,j,k)
#else
                    WRITE(*,*) 'Found value other than 0 or 1 for preprocessor flag <linear>. Now exiting...'
                    STOP
#endif
  
               END DO
            END DO
        END DO
                
    END SUBROUTINE UpdateFluxI
    
    ! Update vertical fluxes using first order upwinding
    SUBROUTINE UpdateFluxJ(Q)
    
        DOUBLE PRECISION, DIMENSION(param%i_min:, param%j_min:, param%k_min:), &
                        INTENT(IN) :: Q
        
        DOUBLE PRECISION :: Uadv
        INTEGER :: i, j, k, thread_nr  
#if(linear==0)
        DOUBLE PRECISION :: switch
#endif     
        
        thread_nr = omp_get_thread_num()
        
        DO k=param%k_start, param%k_end
            DO j=param%j_start, param%j_end+1
                DO i=param%i_start, param%i_end
#if(linear==1)
                    Uadv = 1.0_8
                    FluxJ(i,j,k,thread_nr) = 0.5_8*Uadv*( Q(i,j-1,k) + Q(i,j,k) ) - 0.5_8*ABS(Uadv)*( Q(i,j,k) - Q(i,j-1,k) )
#elif(linear==0)                                    
                    Uadv = 0.5*(Q(i,j-1,k)+Q(i,j,k))
                    switch = 0.5*(SIGN(1.0_8,Uadv)+1.0)
                    FluxJ(i,j,k,thread_nr) = (1.0-switch)*0.5*Q(i,j,k)**2 + switch*0.5*Q(i,j-1,k)**2
#else
                    WRITE(*,*) 'Found value other than 0 or 1 for preprocessor flag <linear>. Now exiting...'
                    STOP
#endif                    

                END DO
            END DO
        END DO
    
    END SUBROUTINE UpdateFluxJ
    
    SUBROUTINE UpdateFluxK(Q)
    
        DOUBLE PRECISION, DIMENSION(param%i_min:, param%j_min:, param%k_min:), &
                        INTENT(IN) :: Q
        
        DOUBLE PRECISION :: Uadv
        INTEGER :: i, j, k, thread_nr  
#if(linear==0)
        DOUBLE PRECISION :: switch
#endif     
        
        thread_nr = omp_get_thread_num()
        
        DO k=param%k_start, param%k_end+1
            DO j=param%j_start, param%j_end
                DO i=param%i_start, param%i_end
                
#if(linear==1)                    
                    Uadv = 1.0_8
                    FluxK(i,j,k,thread_nr) = 0.5_8*Uadv*( Q(i,j,k-1) + Q(i,j,k) ) - 0.5_8*ABS(Uadv)*( Q(i,j,k) - Q(i,j,k-1) )
#elif(linear==0)                                    
                    Uadv = 0.5*(Q(i,j,k-1)+Q(i,j,k))
                    switch = 0.5*(SIGN(1.0_8,Uadv)+1.0)
                    FluxK(i,j,k,thread_nr) = (1.0-switch)*0.5*Q(i,j,k)**2 + switch*0.5*Q(i,j,k-1)**2
#else
                    WRITE(*,*) 'Found value other than 0 or 1 for preprocessor flag <linear>. Now exiting...'
                    STOP
#endif                                    
                END DO
            END DO
        END DO
            
    END SUBROUTINE UpdateFluxK

END MODULE upwind