PROGRAM boundaries_test

USE boundaries, only : periodic, InitializeBoundaries, FinalizeBoundaries

IMPLICIT NONE

DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932_8 ! underscore indicates rounding to real(8) precision

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q, Qall
DOUBLE PRECISION :: dx, dy, dz, x, y, z
INTEGER :: Nx, Ny, Nz, seed, time1(8), i, j, k

CALL DATE_AND_TIME(values=time1)
seed = 1000*time1(7)+time1(8)
CALL SRAND(seed)

! Create random integer values for Nx, Ny, Nz
Nx = 25 + INT(RAND()*256)
Ny = 25 + INT(RAND()*256)
Nz = 25 + INT(RAND()*256)

dx = 1.0/DBLE(Nx)
dy = 1.0/DBLE(Ny)
dz = 1.0/DBLE(Nz)
ALLOCATE(Q(   -2:Nx+3, -2:Ny+3, -2:Nz+3))
ALLOCATE(Qall(-2:Nx+3, -2:Ny+3, -2:Nz+3))

CALL InitializeBoundaries(Nx, Ny, Nz)

! Initialize Q with -100
Q = DBLE(-100)

! Fill inner values of Q
DO k=1,Nz
    DO j=1,Ny
        DO i=1,Nx

            x = 0.5*dx + DBLE(i - 1)*dx
            y = 0.5*dy + DBLE(j - 1)*dy
            z = 0.5*dz + DBLE(k - 1)*dz
    
            Q(i,j,k) = SIN(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z)
        END DO
    END DO
END DO

! Fill Qall, including correctly set halo values
! Fill inner values of Q
DO k=-2,Nz+3
    DO j=-2,Ny+3
        DO i=-2,Nx+3

            x = 0.5*dx + DBLE(i - 1)*dx
            y = 0.5*dy + DBLE(j - 1)*dy
            z = 0.5*dz + DBLE(k - 1)*dz
    
            Qall(i,j,k) = SIN(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z)
        END DO
    END DO
END DO

CALL periodic(Q)

IF (MAXVAL(ABS(Q(-2:0,1:Ny,1:Nz) - Qall(-2:0,1:Ny,1:Nz)))>1e-14) THEN
    WRITE(*,*) 'Mismatch between Q and Qall in i-ghostcells with indices -2,-1,0'
    STOP
ELSE IF (MAXVAL(ABS(Q(Nx+1:Nx+3,1:Ny,1:Nz) - Qall(Nx+1:Nx+3,1:Ny,1:Nz)))>1e-14) THEN
    WRITE(*,*) 'Mismatch between Q and Qall in i-ghostcells with indices Nx+1, Nx+2, Nx+3'
    STOP
ELSE IF (MAXVAL(ABS(Q(1:Nx,-2:0,1:Nz) - Qall(1:Nx,-2:0,1:Nz)))>1e-14) THEN
    WRITE(*,*) 'Mismatch between Q and Qall in j-ghostcells with indices -2,-1,0'
    STOP
ELSE IF (MAXVAL(ABS(Q(1:Nx,Ny+1:Ny+3,1:Nz) - Qall(1:Nx,Ny+1:Ny+3,1:Nz)))>1e-14) THEN
    WRITE(*,*) 'Mismatch between Q and Qall in j-ghostcells with indices Ny+1, Ny+2, Ny+3'
    STOP   
ELSE IF (MAXVAL(ABS(Q(1:Nx,1:Ny,-2:0)-Qall(1:Nx,1:Ny,-2:0)))>1e-14) THEN 
    WRITE(*,*) 'Mismatch between Q and Qall in k-ghostcells with indices -2,-1,0'
    STOP   
ELSE IF (MAXVAL(ABS(Q(1:Nx,1:Ny,Nz+1:Nz+3)-Qall(1:Nx,1:Ny,Nz+1:Nz+3)))>1e-14) THEN 
    WRITE(*,*) 'Mismatch between Q and Qall in k-ghostcells with indices Nz+1, Nz+2, Nz+3'
    STOP      
END IF

IF (MAXVAL(ABS(Q(-2:0,1:Ny,1:Nz) - Q(Nx-2:Nx,1:Ny,1:Nz)))>1e-14) THEN
    WRITE(*,*) 'Mismatch in i-ghostcells with index -2,-1,0 and other side Q values'
    STOP
ELSE IF (MAXVAL(ABS(Q(Nx+1:Nx+3,1:Ny,1:Nz)-Q(1:3,1:Ny,1:Nz)))>1e-14) THEN
    WRITE(*,*) 'Mismatch in i-ghostcells with index Nx+1, Nx+2, Nx+3 and other side Q values'
    STOP  
ELSE IF (MAXVAL(ABS(Q(1:Nx,-2:0,1:Nz)-Q(1:Nx,Ny-2:Ny,1:Nz)))>1e-14) THEN    
    WRITE(*,*) 'Mismatch in j-ghostcells with index -2,-1,0 and other side Q values'
    STOP
ELSE IF (MAXVAL(ABS(Q(1:Nx,Ny+1:Ny+3,1:Nz)-Q(1:Nx,1:3,1:Nz)))>1e-14) THEN
    WRITE(*,*) 'Mismatch in j-ghostcells with index Nx+1, Nx+2, Nx+3 and other side Q values'
    STOP  
ELSE IF (MAXVAL(ABS(Q(1:Nx,1:Ny,-2:0)-Q(1:Nx,1:Ny,Nz-2:Nz)))>1e-14) THEN
    WRITE(*,*) 'Mismatch in k-ghostcells with index -2,-1,0 and other side Q values'
    STOP  
ELSE IF (MAXVAL(ABS(Q(1:Nx,1:Ny,Nz+1:Nz+3)-Q(1:Nx,1:Ny,1:3)))>1e-14) THEN
    WRITE(*,*) 'Mismatch in k-ghostcells with index Nx+1, Nx+2, Nx+3 and other side Q values'
    STOP  
END IF

DEALLOCATE(Q)
DEALLOCATE(Qall)

WRITE(*,*) '**** Boundaries: Test successful ****'

END PROGRAM boundaries_test