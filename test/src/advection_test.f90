!>
!
!
PROGRAM Upwind_Test

USE Advection, only : GetRHSAdvection, InitializeAdvection, FinalizeAdvection
USE omp_lib

IMPLICIT NONE

#include <preprocessor.f90>

INCLUDE 'mpif.h'

INTEGER, PARAMETER :: Nthreads = 8
INTEGER, PARAMETER, DIMENSION(11) :: Nx_v = (/ 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150 /)            
!INTEGER, PARAMETER, DIMENSION(1) :: Nx_v = (/ 25 /)            

DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932_8 ! underscore indicates rounding to real(8) precision

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)       :: err_v, convrate
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)   :: FQ_x_y, error
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, RQ 
DOUBLE PRECISION                                  :: dx, dy, dz, T0, T1, x, y, z, q_val
DOUBLE PRECISION, DIMENSION(Nthreads)             :: runtimes
INTEGER                                           :: i, j, k, r, nt, ierr, mpi_thread_provided, i_add, j_add, k_add

INTEGER :: Nx, Ny, Nz, i_start, i_end, j_start, j_end, k_start, k_end, seed, time1(8), order

LOGICAL, PARAMETER, DIMENSION(3) :: do_test = (/ .true. , .true., .true. /)

CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)

! Create two random value for the offsets
CALL DATE_AND_TIME(values=time1)
seed = 1000*time1(7)+time1(8)
CALL SRAND(seed)

IF (do_test(1)) THEN

    ! Test for order 1 and 5
    DO order=1,5,4

        Nx = 128 + INT( RAND()*128 )
        Ny = 128 + INT( RAND()*128 )
        Nz = 128 + INT( RAND()*128 )

        ! Generate three random integers between -256 and +256
        i_start = INT( RAND()*(256+1+256) )-256
        j_start = INT( RAND()*(256+1+256) )-256
        k_start = INT( RAND()*(256+1+256) )-256

        ! Then generate three more positive random integers between 3 and 32.
        ! In initialization, the used buffer is enlarged by these numbers to verify the routine
        ! can correctly work on a subset of the allocated buffer
        i_add = 3 + INT( RAND()*32 )
        j_add = 3 + INT( RAND()*32 )
        k_add = 3 + INT( RAND()*32 )

        ! The first test makes sure that a constant Q leads to a zero flux divergence. This is
        ! for running 1 to 8 concurrent calls to GetRHS. Also, the time is measured in some 
        ! form of basic weak scaling test and if 8 threads calling 8 instances of GetRHS take
        ! significantly longer than 1 thread calling 1 instance, a warning is issued.

        dx = 1.0/DBLE(Nx)
        dy = 1.0/DBLE(Ny)
        dz = 1.0/DBLE(Nz)

        DO nt=1,Nthreads

            i_end = i_start+Nx-1
            j_end = j_start+Ny-1
            k_end = k_start+Nz-1
    
            ALLOCATE(Q( i_start-i_add:i_end+i_add, j_start-j_add:j_end+j_add, k_start-k_add:k_end+k_add, 0:nt-1))
            ALLOCATE(RQ(i_start-i_add:i_end+i_add, j_start-j_add:j_end+j_add, k_start-k_add:k_end+k_add, 0:nt-1))

            CALL InitializeAdvection(i_start-i_add, i_end+i_add, j_start-j_add, j_end+j_add, k_start-k_add, k_end+k_add, nt, .true.)

            CALL OMP_SET_NUM_THREADS(nt)

            ! fill Q with constant random value
            q_val = RAND(1)
    
            !$OMP PARALLEL DO schedule(static)
            DO i=0,nt-1
                Q(:,:,:,i)  = q_val
        
                ! Initialize RQ with -1
                RQ(:,:,:,i) = -1.0
            END DO
            !$OMP END PARALLEL DO

            ! TEST 1  : Make sure that constant Q gives a zero flux divergence

            T0 = MPI_WTIME()
            !$OMP PARALLEL DO schedule(static)
            DO i=0,nt-1
                CALL GetRHSAdvection(Q(:,:,:,i), RQ(:,:,:,i), dx, dy, dz, i_start, i_end, j_start, j_end, k_start, k_end, order)
            END DO
            !$OMP END PARALLEL DO
            T1 = MPI_WTIME()
        
            IF (MAXVAL(ABS(RQ(i_start:i_end,j_start:j_end,k_start:k_end,0:nt-1)))>1e-12) THEN
                WRITE(*,'(A,I1,A)') 'ERROR: Test failed for order ', order, ' : constant Q does not result in zero flux divergence'
                OPEN(unit=1,file='rqfield.txt')
                    WRITE(1, '(F35.25)') RQ
                STOP
            END IF
    
            runtimes(nt) = T1-T0

            DEALLOCATE(Q)
            DEALLOCATE(RQ)
            CALL FinalizeAdvection

        END DO

        ! If max and min runtime vary by more than 25%, print out the timings and a warning
        IF (MAXVAL(runtimes)/MINVAL(runtimes) > 1.25) THEN
            WRITE(*,'(A, I1, A)') '**** WARNING Order ', order, ': Possible performance bug.... '
            DO nt=1,Nthreads
                WRITE(*,'(A, I2, A, F9.5, A, F5.1)') '#Threads = ', nt, ' -- Runtime: ', runtimes(nt), &
                        ' -- Efficiency: ', 100*runtimes(1)/runtimes(nt)
            END DO
        END IF
    
    END DO

ELSE
    WRITE(*,*) 'Warning: Test 1 is disabled.'
END IF

! TEST 2: Make sure that if Q is constant in one dimension, so is the resulting RQ
IF (do_test(2)) THEN

        Nx = 32 + INT( RAND()*32 )
        Ny = 32 + INT( RAND()*32 )
        Nz = 32 + INT( RAND()*32 )

        ! Generate three random integers between -256 and +256
        i_start = INT( RAND()*(256+1+256) )-256
        j_start = INT( RAND()*(256+1+256) )-256
        k_start = INT( RAND()*(256+1+256) )-256

        ! Then generate three more positive random integers between 3 and 32.
        ! In initialization, the used buffer is enlarged by these numbers to verify the routine
        ! can correctly work on a subset of the allocated buffer
        i_add = 3 + INT( RAND()*32 )
        j_add = 3 + INT( RAND()*32 )
        k_add = 3 + INT( RAND()*32 )

        ! The first test makes sure that a constant Q leads to a zero flux divergence. This is
        ! for running 1 to 8 concurrent calls to GetRHS. Also, the time is measured in some 
        ! form of basic weak scaling test and if 8 threads calling 8 instances of GetRHS take
        ! significantly longer than 1 thread calling 1 instance, a warning is issued.

        dx = 1.0/Nx
        dy = 1.0/Ny
        dz = 1.0/Nz
        
        i_end = i_start+Nx-1
        j_end = j_start+Ny-1
        k_end = k_start+Nz-1

        ALLOCATE(Q( i_start-i_add:i_end+i_add, j_start-j_add:j_end+j_add, k_start-k_add:k_end+k_add, 0:Nthreads-1))
        ALLOCATE(RQ(i_start-i_add:i_end+i_add, j_start-j_add:j_end+j_add, k_start-k_add:k_end+k_add, 0:Nthreads-1))

        CALL InitializeAdvection(i_start-i_add, i_end+i_add, j_start-j_add, j_end+j_add, k_start-k_add, k_end+k_add, Nthreads, .true.)
        
        ! First, test along i axis
        DO order=1,5,4
            ! Set Q
            DO k=k_start-k_add, k_end+k_add
                DO j=j_start-j_add,j_end+j_add
                    DO i=i_start-i_add, i_end+i_add
                        x = 0.5*dx + DBLE(i - i_start)*dx               
                        Q(i,j,k,:)    = SIN(2.0*pi*x)
                    END DO
                END DO
            END DO
            
            !$OMP PARALLEL DO schedule(static)
            DO nt=0,Nthreads-1
                CALL GetRHSAdvection(Q(:,:,:,nt), RQ(:,:,:,nt), dx, dy, dz, i_start, i_end, j_start, j_end, k_start, k_end, order)
            END DO
            !$OMP END PARALLEL DO
            
            DO j=j_start,j_end-1
                IF (MAXVAL(ABS(RQ(:,j+1,:,:)-RQ(:,j,:,:)))>1e-13) THEN
                    WRITE(*,'(A, I1, A)') 'ERROR: (IJ) -- Constant Q in j direction produced non-constant RQ for order = ', order, '. Now exiting...'
                    OPEN(unit=1,file='q.txt')
                    WRITE(1,'(F35.25)') Q(1:Nx, 1:Ny, 1:Nz, 0)
                    CLOSE(1)
                    STOP               
                END IF 
            END DO
            
            DO k=k_start,k_end-1
                IF (MAXVAL(ABS(RQ(:,:,k+1,:)-RQ(:,:,k,:)))>1e-13) THEN
                    WRITE(*,'(A, I1, A)') 'ERROR: (IK) -- Constant Q in k direction produced non-constant RQ for order = ', order, '. Now exiting...'
                    STOP               
                END IF 
            END DO                             
        END DO
        
        ! Second, test along j axis
        DO order=1,5,4
            ! Set Q
           DO k=k_start-k_add, k_end+k_add
                DO j=j_start-j_add,j_end+j_add
                    DO i=i_start-i_add, i_end+i_add
                        y = 0.5*dy + DBLE(j - j_start)*dy
                        Q(i,j,k,:) = SIN(2.0*pi*y)
                    END DO
                END DO
            END DO        
            
            !$OMP PARALLEL DO schedule(static)
            DO nt=0,Nthreads-1
                CALL GetRHSAdvection(Q(:,:,:,nt), RQ(:,:,:,nt), dx, dy, dz, i_start, i_end, j_start, j_end, k_start, k_end, order)
            END DO
            !$OMP END PARALLEL DO
            
            DO i=i_start,i_end-1
                IF (MAXVAL(ABS(RQ(i+1,:,:,:)-RQ(i,:,:,:)))>1e-13) THEN
                    WRITE(*,'(A, I1, A)') 'ERROR: (JI) -- Constant Q in i direction produced non-constant RQ for order = ', order, '. Now exiting...'
                END IF 
            END DO
            
            DO k=k_start,k_end-1
                IF (MAXVAL(ABS(RQ(:,:,k+1,:)-RQ(:,:,k,:)))>1e-13) THEN
                    WRITE(*,'(A, I1, A)') 'ERROR: (JK) -- Constant Q in k direction produced non-constant RQ for order = ', order, '. Now exiting...'
                END IF             
            END DO            
        END DO      
          
        ! Third, test along k axis
        DO order=1,5,4
        
            ! Set Q
           DO k=k_start-k_add, k_end+k_add
                DO j=j_start-j_add,j_end+j_add
                    DO i=i_start-i_add, i_end+i_add
                        z = 0.5*dz + DBLE(k - k_start)*dz
                        Q(i,j,k,:)    = SIN(2.0*pi*z)
                    END DO
                END DO
            END DO        
            
            !$OMP PARALLEL DO schedule(static)
            DO nt=0,Nthreads-1
                CALL GetRHSAdvection(Q(:,:,:,nt), RQ(:,:,:,nt), dx, dy, dz, i_start, i_end, j_start, j_end, k_start, k_end, order)
            END DO
            !$OMP END PARALLEL DO
            
            DO i=i_start,i_end-1
                IF (MAXVAL(ABS(RQ(i+1,:,:,:)-RQ(i,:,:,:)))>1e-13) THEN
                    WRITE(*,'(A, I1, A)') 'ERROR: (KI) -- Constant Q in i direction produced non-constant RQ for order = ', order, '. Now exiting...'
                END IF             
            END DO
            
            DO j=j_start,j_end-1
                IF (MAXVAL(ABS(RQ(:,j+1,:,:)-RQ(:,j,:,:)))>1e-13) THEN
                    WRITE(*,'(A, I1, A)') 'ERROR: (KJ) -- Constant Q in j direction produced non-constant RQ for order = ', order, '. Now exiting...'
                END IF             
            END DO            
        END DO
        CALL FinalizeAdvection
        
        DEALLOCATE(Q)
        DEALLOCATE(RQ)
        
ELSE
    WRITE(*,*) 'Warning: Test 2 is disabled.'
END IF

! TEST 3: Test convergence order
IF (do_test(3)) THEN

    ! Test order 1 and 5
    DO order=1,5,4
    
        ALLOCATE(err_v(SIZE(Nx_v)))
        ALLOCATE(convrate(SIZE(Nx_v)-1))
        DO r=1,SIZE(Nx_v)
            ! Fetch resolutions and add a small number of points to Nx and Ny, in order to
            ! have slightly different mesh widths along different axes
            Nx = Nx_v(r)+3   
            Ny = Nx_v(r)+2
            Nz = Nx_v(r)

            dx = 1.0/DBLE(Nx)
            dy = 1.0/DBLE(Ny)
            dz = 1.0/DBLE(Nz)

            ! Generate three random integers between -256 and +256
            i_start = INT( RAND()*(256+1+256) )-256
            j_start = INT( RAND()*(256+1+256) )-256
            k_start = INT( RAND()*(256+1+256) )-256

            ! Then generate three more positive random integers between 3 and 32.
            ! In initialization, the used buffer is enlarged by these numbers to verify the routine
            ! can correctly work on a subset of the allocated buffer
            i_add = 3 + INT( RAND()*8 )
            j_add = 3 + INT( RAND()*8 )
            k_add = 3 + INT( RAND()*8 )

            i_end = i_start + Nx - 1
            j_end = j_start + Ny - 1
            k_end = k_start + Nz - 1
    
            ALLOCATE (Q( i_start-i_add:i_end+i_add, j_start-j_add:j_end+j_add, k_start-k_add:k_end+k_add, 0:Nthreads-1))
            ALLOCATE(RQ( i_start-i_add:i_end+i_add, j_start-j_add:j_end+j_add, k_start-k_add:k_end+k_add, 0:Nthreads-1))
    
            ! Exact solution for d/dx f(q) + d/dy f(q) + d/dz f(q) with f(q) = q
            ALLOCATE(FQ_x_y( i_start:i_end, j_start:j_end, k_start:k_end))
            ALLOCATE(error(  i_start:i_end, j_start:j_end, k_start:k_end ))
    
            CALL InitializeAdvection(i_start-i_add, i_end+i_add, j_start-j_add, j_end+j_add, k_start-k_add, k_end+k_add, Nthreads, .true.)

            ! Write meaningful data only in a subset of the allocated Q
            DO k=k_start-k_add, k_end+k_add
                DO j=j_start-j_add,j_end+j_add
                    DO i=i_start-i_add, i_end+i_add
        
                        x = 0.5*dx + DBLE(i - i_start)*dx
                        y = 0.5*dy + DBLE(j - j_start)*dy
                        z = 0.5*dz + DBLE(k - k_start)*dz
                
                        Q(i,j,k,:)    = SIN(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z)
                    END DO
                END DO
            END DO

            DO k=k_start, k_end
                DO j=j_start, j_end
                    DO i=i_start, i_end
            
                        x = 0.5*dx + DBLE(i - i_start)*dx
                        y = 0.5*dy + DBLE(j - j_start)*dy
                        z = 0.5*dz + DBLE(k - k_start)*dz
#if(linear==1)                                   
                        FQ_x_y(i,j,k) = -2.0*pi*COS(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z) &
                                        -2.0*pi*COS(2.0*pi*y)*SIN(2.0*pi*x)*SIN(2.0*pi*z) &
                                        -2.0*pi*COS(2.0*pi*z)*SIN(2.0*pi*x)*SIN(2.0*pi*y)

#elif(linear==0)
                         FQ_x_y(i,j,k) = -2.0*pi*COS(2.0*pi*x)*SIN(2.0*pi*x) * SIN(2.0*pi*y)**2 * SIN(2.0*pi*z)**2  &
                                         -2.0*pi*COS(2.0*pi*y)*SIN(2.0*pi*y) * SIN(2.0*pi*x)**2 * SIN(2.0*pi*z)**2 &
                                         -2.0*pi*COS(2.0*pi*z)*SIN(2.0*pi*z) * SIN(2.0*pi*x)**2 * SIN(2.0*pi*y)**2 
#endif                                         
                    END DO
                END DO
            END DO   

            !$OMP PARALLEL DO schedule(static)
            DO nt=0,Nthreads-1
                CALL GetRHSAdvection(Q(:,:,:,nt), RQ(:,:,:,nt), dx, dy, dz, i_start, i_end, j_start, j_end, k_start, k_end, order)
            END DO
            !$OMP END PARALLEL DO

            DO k=k_start, k_end
                DO j=j_start, j_end
                    DO i=i_start, i_end
                        error(i,j,k) = ABS(RQ(i,j,k,Nthreads-1) - FQ_x_y(i,j,k))   
                    END DO
                END DO             
            END DO      
    
            err_v(r) = maxval(error)/maxval(abs(FQ_x_y))
        
            CALL FinalizeAdvection

            DEALLOCATE(Q)
            DEALLOCATE(RQ)
            DEALLOCATE(FQ_x_y)
            DEALLOCATE(error)
    
        END DO

        DO r=2,SIZE(Nx_v)
            convrate(r-1) = LOG10(err_v(r)/err_v(r-1)) / LOG10( DBLE(Nx_v(r-1)) / DBLE(Nx_v(r)) )
        END DO
        
        IF (minval(convrate)<=0.95*order) THEN
            WRITE(*,'(A, I2, A)') 'ERROR: Failed to verify convergence rate of ', order, ' ... printing errors and convergence rates, then stopping'
            DO r=1,SIZE(Nx_v)
                WRITE(*,'(ES9.3)') err_v(r)
            END DO
            DO r=1,SIZE(Nx_v)-1
                WRITE(*,'(A, F5.3)') 'Convergence rate: ', convrate(r)
            END DO          
            STOP
        END IF
        
        DEALLOCATE(err_v)
        DEALLOCATE(convrate)
        
    END DO    
        
ELSE
    WRITE(*,*) 'Warning: Test 3 is disabled.'

END IF

CALL MPI_FINALIZE(ierr)

IF(SIZE(Nx_v)==1) THEN
    WRITE(*,*) 'WARNING: No convergence test... output in ascii files'
ELSE
    WRITE(*,*) '[0] -- Successful: Advection module returns zero for constant input and produces expected rates of convergence.'
END IF

END PROGRAM Upwind_Test