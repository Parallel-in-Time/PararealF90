!>
!! This module handles the diffusion terms in the equation. It provides a second and fourth order centred discretization of
!!
!! \\( \\Delta q = q_{xx} + q_{yy} + q_{zz} \\)
!!
MODULE diffusion

IMPLICIT NONE

PUBLIC :: GetRHSDiffusion, InitializeDiffusion

TYPE diffusion_parameter
    INTEGER :: i_min, i_max, j_min, j_max, k_min, k_max
    DOUBLE PRECISION :: nu
END TYPE

TYPE(diffusion_parameter) :: param

!> Finite difference weight in the fourth order centred stencil
DOUBLE PRECISION, PARAMETER :: w1 = 1.0_8/12.0_8

!> Finite difference weight in the fourth order centred stencil
DOUBLE PRECISION, PARAMETER :: w2 = 4.0_8/3.0_8

!> Finite difference weight in the fourth order centred stencil
DOUBLE PRECISION, PARAMETER :: w3 = 5.0_8/2.0_8

CONTAINS

    !> Computes the contribution from the diffusive terms to the right hand side of the initial value problem.
    !! param[in] Q The solution
    !! param[out] RQ Discrete diffusion operator applied to Q
    !! param[in] order Either two or four.
    SUBROUTINE GetRHSDiffusion(Q, RQ, dx, dy, dz, i_start, i_end, j_start, j_end, k_start, k_end, order)
        DOUBLE PRECISION, DIMENSION(param%i_min:, param%j_min:, param%k_min:), INTENT(IN)  :: Q
        DOUBLE PRECISION, DIMENSION(param%i_min:, param%j_min:, param%k_min:), INTENT(OUT) :: RQ
        DOUBLE PRECISION,   INTENT(IN) :: dx, dy, dz
        INTEGER,            INTENT(IN) :: i_start, i_end, j_start, j_end, k_start, k_end, order
        
        DOUBLE PRECISION :: stencil_i, stencil_j, stencil_k, coeff_i, coeff_j, coeff_k
        
        INTEGER :: i, j, k
        
        SELECT CASE (order)
        
            CASE (2)
            
                coeff_i = 1.0/(dx*dx)
                coeff_j = 1.0/(dy*dy)
                coeff_k = 1.0/(dz*dz)
                
                DO k=k_start, k_end
                    DO j=j_start, j_end
                        DO i=i_start, i_end
                            stencil_i = coeff_i*(Q(i+1,j,k) - 2.0*Q(i,j,k) + Q(i-1,j,k))
                            stencil_j = coeff_j*(Q(i,j+1,k) - 2.0*Q(i,j,k) + Q(i,j-1,k))
                            stencil_k = coeff_k*(Q(i,j,k+1) - 2.0*Q(i,j,k) + Q(i,j,k-1))
                            RQ(i,j,k) = RQ(i,j,k) +  param%nu*(stencil_i + stencil_j + stencil_k)
                        END DO
                    END DO
                END DO  
                
            CASE (4)
            
                coeff_i = 1.0/(dx*dx)
                coeff_j = 1.0/(dy*dy)
                coeff_k = 1.0/(dz*dz)
                
                DO k=k_start, k_end
                    DO j=j_start, j_end
                        DO i=i_start, i_end
                            stencil_i = coeff_i*( -w1*Q(i+2,j,k) + w2*Q(i+1,j,k) &
                                - w3*Q(i,j,k) + w2*Q(i-1,j,k) - w1*Q(i-2,j,k) )
                                
                            stencil_j = coeff_j*( -w1*Q(i,j+2,k) + w2*Q(i,j+1,k) &
                                -w3*Q(i,j,k) + w2*Q(i,j-1,k) - w1*Q(i,j-2,k) )
                            
                            stencil_k = coeff_k*( -w1*Q(i,j,k+2) + w2*Q(i,j,k+1) &
                                -w3*Q(i,j,k) + w2*Q(i,j,k-1) - w1*Q(i,j,k-2) )
                            
                            RQ(i,j,k) = RQ(i,j,k) + param%nu*( stencil_i + stencil_j + stencil_k )    
                        END DO
                    END DO
                END DO
            
            CASE DEFAULT
                WRITE(*,'(A,I2,A)') 'No implementation available for order = ', order, ' .. now exiting'
                STOP                
        END SELECT 
        
    END SUBROUTINE GetRHSDiffusion

    !> Initialize the diffusion module
    SUBROUTINE InitializeDiffusion(i_min, i_max, j_min, j_max, k_min, k_max, nu)
    
        DOUBLE PRECISION, INTENT(IN) :: nu        
        INTEGER, INTENT(IN) :: i_min, i_max, j_min, j_max, k_min, k_max
        
        param%nu    = nu 
        param%i_min = i_min
        param%i_max = i_max
        param%j_min = j_min
        param%j_max = j_max
        param%k_min = k_min
        param%k_max = k_max
        
    END SUBROUTINE InitializeDiffusion
    
END MODULE diffusion