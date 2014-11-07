!>
! Fills a 3D field including ghost-cells with values from some given function.
! Then ghost-cell values are also generated from inner values using the boundaries module
! and both results are compared to make sure they match.
!
PROGRAM advection_boundaries_test

USE Advection, only : GetRHSAdvection, InitializeAdvection, FinalizeAdvection
USE boundaries, only : periodic, InitializeBoundaries, FinalizeBoundaries

IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932_8 ! underscore indicates rounding to real(8) precision
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, RQ, Qbc, RQbc
DOUBLE PRECISION :: dx, dy, dz, x, y, z
INTEGER :: Nx, Ny, Nz, i, j, k, order, Nthreads, nt, seed, time1(8)

CALL DATE_AND_TIME(values=time1)
seed = 1000*time1(7)+time1(8)
CALL SRAND(seed)

! Generate random values of Nx, Ny, Nz that are at least 6, to exclude too small domains
Nx = 6 + INT( RAND()*128 )
Ny = 6 + INT( RAND()*128 )
Nz = 6 + INT( RAND()*128 )
Nthreads = 8

ALLOCATE(Q(   -2:Nx+3,-2:Ny+3,-2:Nz+3,0:Nthreads-1))
ALLOCATE(RQ(  -2:Nx+3,-2:Ny+3,-2:Nz+3,0:Nthreads-1))
ALLOCATE(Qbc( -2:Nx+3,-2:Ny+3,-2:Nz+3,0:Nthreads-1))
ALLOCATE(RQbc(-2:Nx+3,-2:Ny+3,-2:Nz+3,0:Nthreads-1))

CALL InitializeAdvection(-2, Nx+3, -2, Ny+3, -2, Nz+3, Nthreads, .false.)
CALL InitializeBoundaries(Nx, Ny, Nz)

DO order=1,5,4
    
    dx = 1.0/DBLE(Nx)
    dy = 1.0/DBLE(Ny)
    dz = 1.0/DBLE(Nz)
    
    ! Initialize Q and Qb with large negative values far away from actual q(x,y,z)
    Q   = -100.0
    Qbc = -100.0
    
    ! Fill all entries of Q
    DO k=-2,Nz+3
        DO j=-2,Ny+3
            DO i=-2,Nx+3
                x = 0.5*dx + DBLE(i - 1)*dx
                y = 0.5*dy + DBLE(j - 1)*dy
                z = 0.5*dz + DBLE(k - 1)*dz
                
                Q(i,j,k,:) = SIN(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z)
            END DO
        END DO
    END DO
    
    ! For Qbc, fill only inner entries
    DO k=1,Nz
        DO j=1,Ny
            DO i=1,Nx
                x = 0.5*dx + DBLE(i - 1)*dx
                y = 0.5*dy + DBLE(j - 1)*dy
                z = 0.5*dz + DBLE(k - 1)*dz
                     
                Qbc(i,j,k,:) = SIN(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z)                            
            END DO
        END DO
    END DO
    
    ! And then call BC routine
    !$OMP PARALLEL DO schedule(static)
    DO nt=0,Nthreads-1
        CALL periodic(Qbc(:,:,:,nt))    
    END DO
    !$OMP END PARALLEL DO
        
    ! Apply advection stencil to both Q and Qbc
    
    !$OMP PARALLEL DO schedule(static)
    DO nt=0,Nthreads-1
        CALL GetRHSAdvection(Q(:,:,:,nt), RQ(:,:,:,nt), dx, dy, dz, 1, Nx, 1, Ny, 1, Nz, order)
    END DO
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO schedule(static)
    DO nt=0,Nthreads-1
        CALL GetRHSAdvection(Qbc(:,:,:,nt), RQbc(:,:,:,nt), dx, dy, dz, 1, Nx, 1, Ny, 1, Nz, order)
    END DO
    !$OMP END PARALLEL DO    
    
    IF (MAXVAL(ABS(RQ(1:Nx,1:Ny,1:Nz,:)-RQbc(1:Nx,1:Ny,1:Nz,:)))>1e-13) THEN
        WRITE(*,'(A, I1, A)') 'ERROR: Mismatch in combination of boundaries and advection module for order=', order, '. Now exiting...'
        OPEN(unit=1,file='rq.txt')
        WRITE(1,'(F35.25)') RQ(1:Nx,1:Ny,1:Nz,0)
        CLOSE(1)
        OPEN(unit=1,file='rqbc.txt')
        WRITE(1,'(F35.25)') RQbc(1:Nx,1:Ny,1:Nz,0)
        CLOSE(1)
        STOP
    END IF
END DO

DEALLOCATE(Q)
DEALLOCATE(Qbc)
DEALLOCATE(RQ)
DEALLOCATE(RQbc)
CALL FinalizeAdvection
CALL FinalizeBoundaries

PRINT*, '\x1B[32m[0] -- Successful: GetRHSAdvection produces identical results when ghost-cells are filled manually or through the boundaries module.\x1B[0m'

END PROGRAM advection_boundaries_test