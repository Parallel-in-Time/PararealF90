PROGRAM diffusion_test

USE omp_lib
USE diffusion, only : GetRHSDiffusion, InitializeDiffusion

IMPLICIT NONE

INCLUDE 'mpif.h'

DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932_8 ! underscore indicates rounding to real(8) precision
INTEGER, PARAMETER :: Nthreads = 1

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, RQ, RQex
DOUBLE PRECISION :: dx, dy, dz, x, y, z
INTEGER :: mpi_thread_provided, ierr, Nx, Ny, Nz, i, j, k, nt, &
    i_start, j_start, k_start, i_end, j_end, k_end, order, nn

INTEGER, DIMENSION(6), PARAMETER :: N_v = (/ 10, 100, 200, 300, 400, 500 /)
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: err_v, convrate


CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)

DO order=2,4,2

    i_start = 1
    j_start = 1
    k_start = 1

    ALLOCATE(err_v(1:SIZE(N_v)))
    ALLOCATE(convrate(1:SIZE(N_v)-1))

    DO nn=1,SIZE(N_v)

        Nx = N_v(nn)
        Ny = N_v(nn)
        Nz = N_v(nn)

        i_end = i_start + Nx - 1
        j_end = j_start + Ny - 1
        k_end = k_start + Nz - 1

        dx = 1.0/DBLE(Nx)
        dy = 1.0/DBLE(Ny)
        dz = 1.0/DBLE(Nz)

        ALLOCATE(Q(   i_start-2:i_end+2, j_start-2:j_end+2, k_start-2:k_end+2, 0:Nthreads-1))
        ALLOCATE(RQ(  i_start-2:i_end+2, j_start-2:j_end+2, k_start-2:k_end+2, 0:Nthreads-1))
        ALLOCATE(RQex(i_start-2:i_end+2, j_start-2:j_end+2, k_start-2:k_end+2, 0:Nthreads-1))

        CALL InitializeDiffusion(i_start-2, i_end+2, j_start-2, j_end+2, k_start-2, k_end+2, DBLE(1.0))
        !$OMP PARALLEL DO schedule(static)
        DO nt=0,Nthreads-1
            Q(:,:,:,nt)  = 0.0
            RQ(:,:,:,nt) = 0.0
        END DO
        !$OMP END PARALLEL DO

        DO k=k_start-2,k_end+2
            DO j=j_start-2,j_end+2
                DO i=i_start-2,i_end+2
                    x = 0.5*dx + DBLE(i - i_start)*dx
                    y = 0.5*dy + DBLE(j - j_start)*dy
                    z = 0.5*dz + DBLE(k - k_start)*dz    
                    
                    Q(i,j,k,:)    = SIN(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z)
                    RQex(i,j,k,:) = -12.0*pi*pi*SIN(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z)
                            
                END DO
            END DO
        END DO

        !$OMP PARALLEL DO schedule(static)
        DO nt=0,Nthreads-1
            CALL GetRHSDiffusion(Q(:,:,:,nt), RQ(:,:,:,nt), dx, dy, dz, i_start, i_end, j_start, j_end, k_start, k_end, order)
        END DO
        !$OMP END PARALLEL DO

        err_v(nn) = MAXVAL(ABS(RQ(i_start:i_end,j_start:j_end,k_start:k_end,:) & 
            -RQex(i_start:i_end,j_start:j_end,k_start:k_end,:)))/MAXVAL(ABS(RQex(i_start:i_end,j_start:j_end,k_start:k_end,:)))
       
        DEALLOCATE(Q)
        DEALLOCATE(RQ)
        DEALLOCATE(RQex)

    END DO

    DO nn=2,SIZE(N_v)
        convrate(nn-1) = LOG10(err_v(nn)/err_v(nn-1)) / LOG10( DBLE(N_v(nn-1)) / DBLE(N_v(nn)) )
    END DO
    
    IF (MINVAL(convrate)<=0.95*order) THEN
        WRITE(*,'(A, I1)') 'Failed to verify convergence order for order = ', order    
        DO nn=1,SIZE(N_v)
            WRITE(*,'(F9.3)') convrate(nn)
        END DO
        STOP
    END IF
    
    DEALLOCATE(err_v)
    DEALLOCATE(convrate)

END DO

CALL MPI_FINALIZE(ierr)

WRITE(*,*) '**** Diffusion: Test successful ****'

END PROGRAM diffusion_test