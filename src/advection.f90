!>
!! @todo docu
!!
MODULE advection

USE omp_lib, only : omp_get_thread_num
USE upwind,  only : UpwindFluxes
USE weno5,   only : WenoFluxes
USE fluxes,  only : InitializeFluxes, FinalizeFluxes, param, FluxI, FluxJ, FluxK

PRIVATE
PUBLIC :: GetRHSAdvection, InitializeAdvection, FinalizeAdvection

CONTAINS

    !> @todo docu
    SUBROUTINE GetRHSAdvection(Q, RQ, dx, dy, dz, i_start, i_end, j_start, j_end, k_start, k_end, order)

        DOUBLE PRECISION, DIMENSION(param%i_min:, param%j_min:, param%k_min:), INTENT(IN)  :: Q
        DOUBLE PRECISION, DIMENSION(param%i_min:, param%j_min:, param%k_min:), INTENT(OUT) :: RQ
        DOUBLE PRECISION, INTENT(IN)  :: dx, dy, dz
        INTEGER,          INTENT(IN)  :: i_start, i_end, j_start, j_end, k_start, k_end, order

        ! The current offset values better be in the range defined by the maximum offsets during initialization, otherwise
        ! you will get nice array out of bounds errors.
        param%i_start = i_start
        param%i_end   = i_end
        param%j_start = j_start
        param%j_end   = j_end
        param%k_start = k_start
        param%k_end   = k_end

        SELECT CASE (order)
        
            CASE (1)
                CALL UpwindFluxes(Q)
            CASE (5)
                CALL WenoFluxes(Q)
            CASE DEFAULT
                WRITE(*,'(A,I2,A)') 'Now implementation available for order = ', order, ' .. now exiting'
                STOP
        END SELECT
                
       ! Compute flux divergence
        CALL GetFluxDivergence(RQ, dx, dy, dz)

    END SUBROUTINE GetRHSAdvection

    !> @todo
    !! @param[in]
    !! @param[in]
    SUBROUTINE InitializeAdvection(i_min, i_max, j_min, j_max, k_min, k_max, Nthreads, echo_on)
    
        INTEGER, INTENT(IN) :: Nthreads, i_min, i_max, j_min, j_max, k_min, k_max
        LOGICAL, INTENT(IN) :: echo_on

        CALL InitializeFluxes(i_min, i_max, j_min, j_max, k_min, k_max, Nthreads, echo_on)

    END SUBROUTINE InitializeAdvection

    !> @todo
    SUBROUTINE FinalizeAdvection()
        CALL FinalizeFluxes
    END SUBROUTINE FinalizeAdvection

    !> For computed fluxes, compute the flux divergence @todo
    !! \\( f(q)_x + f(q)_y + f(q)_z \\)
    !! @param[in]
    !! @param[out]
    SUBROUTINE GetFluxDivergence(RQ, dx, dy, dz)
    
        DOUBLE PRECISION, DIMENSION(param%i_min:, param%j_min:, param%k_min:), INTENT(OUT) :: RQ
        DOUBLE PRECISION, INTENT(IN)  :: dx, dy, dz
        
        DOUBLE PRECISION :: coeff_x, coeff_y, coeff_z
        INTEGER          :: i, j, k, thread_nr
        
        coeff_x   = 1.0_8/dx
        coeff_y   = 1.0_8/dy
        coeff_z   = 1.0_8/dz
        
        thread_nr = omp_get_thread_num()
              
        DO k=param%k_start,param%k_end
            DO j=param%j_start,param%j_end
                DO i=param%i_start,param%i_end
                    RQ(i,j,k) = -coeff_x*( FluxI(i+1,j,   k,   thread_nr) - FluxI(i,j,k,thread_nr) ) &
                                -coeff_y*( FluxJ(i,  j+1, k,   thread_nr) - FluxJ(i,j,k,thread_nr) ) &
                                -coeff_z*( FluxK(i,  j,   k+1, thread_nr) - FluxK(i,j,k,thread_nr) )
                END DO
            END DO
        END DO
               
    END SUBROUTINE GetFluxDivergence

END MODULE advection