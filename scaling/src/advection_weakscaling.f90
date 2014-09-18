PROGRAM advection_weakscaling

USE Advection, only : GetRHSAdvection, InitializeAdvection, FinalizeAdvection
USE OMP_LIB

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, PARAMETER :: Nthreads = 8, Nx_factor = 10, Nx_base = 20
            
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, RQ 
DOUBLE PRECISION                                  :: dx, dy, dz, T0, T1
DOUBLE PRECISION, DIMENSION(Nthreads)             :: runtimes, efficiency
INTEGER                                           :: i, r, nt, ierr, mpi_thread_provided

INTEGER :: Nx, Ny, Nz, i_end, j_end, i_start, j_start, k_start, k_end, order

CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)
i_start = 1
j_start = 1
k_start = 1

DO order=1,5,4

    DO r=1,Nx_factor

        Nx = r*Nx_base
        Ny = r*Nx_base
        Nz = r*Nx_base
        dx = 1.0/Nx
        dy = 1.0/Ny
        dz = 1.0/Nz

        DO nt=1,Nthreads

            i_end = i_start+Nx-1
            j_end = j_start+Ny-1
            k_end = k_start+Nz-1
    
            ALLOCATE(Q( j_start-3:j_end+3,   i_start-3:i_end+3,  k_start-3:k_end+3, 0:nt-1))
            ALLOCATE(RQ(j_start-3:j_end+3,   i_start-3:i_end+3,  k_start-3:k_end+3, 0:nt-1))

            CALL InitializeAdvection(i_start-3, i_end+3, j_start-3, j_end+3, k_start-3, k_end+3, nt, .true.)

            CALL OMP_SET_NUM_THREADS(nt)
    
            !$OMP PARALLEL DO schedule(static)
            DO i=0,nt-1
                Q(:,:,:,i)  = 3.33
        
                ! Initialize RQ with -1
                RQ(:,:,:,i) = -1.0
            END DO
            !$OMP END PARALLEL DO

            T0 = MPI_WTIME()
        
            !$OMP PARALLEL DO schedule(static)
            DO i=0,nt-1
                CALL GetRHSAdvection(Q(:,:,:,i), RQ(:,:,:,i), dx, dy, dz, i_start, i_end, j_start, j_end, k_start, k_end, order)
            END DO
            !$OMP END PARALLEL DO
            T1 = MPI_WTIME()
   
            runtimes(nt) = T1-T0

            DEALLOCATE(Q)
            DEALLOCATE(RQ)
            CALL FinalizeAdvection()

        END DO
    
        DO nt=1,Nthreads
            efficiency(nt) = (runtimes(1)/runtimes(nt))*100
        END DO
    
        WRITE(*,*) 
        WRITE(*,*) '========================================'
        WRITE(*,'(A, I4, A, I4, A, I4)') 'Problem size : ', Nx, ' x ', Ny, ' x ', Nz
       ! DO nt=1,Nthreads
        !        WRITE(*,'(A, I1, A, I2, A, F5.3)') 'Order ', order, ' -- Runtime for N=', nt , ' threads : ', runtimes(nt)
        !END DO
        WRITE(*,*)
        DO nt=1,Nthreads
            WRITE(*,'(A, I1, A, I2, A, F5.1)') 'Order ', order, ' -- Efficiency for N=', nt , ' threads : ', efficiency(nt)
        END DO
    END DO

END DO

CALL MPI_FINALIZE(ierr)

END PROGRAM advection_weakscaling