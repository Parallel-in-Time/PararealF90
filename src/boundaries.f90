!>
!! @todo docu
!!
MODULE boundaries

IMPLICIT NONE

PRIVATE
PUBLIC :: periodic, InitializeBoundaries, FinalizeBoundaries

!> @todo docu
TYPE boundaries_parameter
    INTEGER :: i_max, j_max, k_max
END TYPE

!> @todo docu
TYPE(boundaries_parameter) :: param

CONTAINS

    !> @todo docu
    !! @param[inout]
    SUBROUTINE periodic(Q)
        DOUBLE PRECISION, DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: Q
        
        INTEGER :: i, j, k

        ! The following indices indicate identical ghost-cells for periodic BC with a 3-cell halo
        ! -2 <-> Nx-2
        ! -1 <-> Nx-1
        !  0 <-> Nx
        ! Therefore, use the mapping i -> Nx+i, -2<=i<=0. 
        ! In the other direction it is
        ! Nx+1 <-> 1
        ! Nx+2 <-> 2
        ! Nx+3 <-> 3
        ! and therefore the mapping is i -> i-Nx
        
        ! Ghost-cells in i direction
        DO k=1,param%k_max
            DO j=1,param%j_max
                DO i=-2,0
                      Q(i,j,k) = Q(param%i_max+i,j,k)
                END DO
            END DO
        END DO
        
        DO k=1,param%k_max
            DO j=1,param%j_max
                 DO i=param%i_max+1,param%i_max+3
                      Q(i,j,k) = Q(i-param%i_max,j,k)
                 END DO
            END DO
        END DO
        
        ! Ghost-cells in j direction
        DO k=1,param%k_max
            DO j=-2,0
                DO i=1,param%i_max
                    Q(i,j,k) = Q(i,param%j_max+j,k)
                END DO
            END DO
        END DO
        
        DO k=1,param%k_max
            DO j=param%j_max+1,param%j_max+3
                DO i=1,param%i_max
                    Q(i,j,k) = Q(i,j-param%j_max,k)
                END DO
            END DO
        END DO
        
        ! Ghost-cell in k direction
        DO k=-2,0
            DO j=1,param%j_max
                DO i=1,param%i_max
                    Q(i,j,k) = Q(i,j,param%k_max+k)
                END DO
            END DO
        END DO                 
        
        DO k=param%k_max+1,param%k_max+3
            DO j=1,param%j_max
                DO i=1,param%i_max
                    Q(i,j,k) = Q(i,j,k-param%k_max)
                END DO
            END DO
        END DO  
                         
    END SUBROUTINE periodic

    !> @todo docu
    SUBROUTINE InitializeBoundaries(i_max, j_max, k_max)
        INTEGER, INTENT(IN) :: i_max, j_max, k_max
        
        param%i_max = i_max
        param%j_max = j_max
        param%k_max = k_max
        
    END SUBROUTINE InitializeBoundaries

    !> @todo docu
    SUBROUTINE FinalizeBoundaries
      ! Nothing to DEALLOCATE here.
    END SUBROUTINE FinalizeBoundaries
    
END MODULE boundaries