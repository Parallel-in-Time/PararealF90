!>
!
!
PROGRAM diffusion_threads_test

USE diffusion, only : GetRHSDiffusion, InitializeDiffusion

DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932_8 ! underscore indicates rounding to real(8) precision

INTEGER order, Nx, Ny, Nz, Nthreads, i, j, k, nt

INTEGER, DIMENSION(5) :: N_v = (/ 20, 30, 40, 50, 60 /)

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: Q, RQ

DOUBLE PRECISION :: dx, dy, dz, x, y, z

DO order=2,4,2
  DO nn=1,SIZE(N_v)
    Nx = N_v(nn)
    Ny = Nx
    Nz = Nx

    dx = 1.0/DBLE(Nx)
    dy = 1.0/DBLE(Ny)
    dz = 1.0/DBLE(Nz)

    DO Nthreads=1,2

      ALLOCATE(Q( -2:Nx+3,-2:Ny+3,-2:Nz+3,0:Nthreads-1))
      ALLOCATE(RQ(-2:Nx+3,-2:Ny+3,-2:Nz+3,0:Nthreads-1))


      CALL InitializeDiffusion(-2,Nx+3,-2,Ny+3,-2,Nz+3, DBLE(0.05))

      ! Fill Q
      DO i=-2,Nx+3
        DO j=-2,Ny+3
          DO k=-2,Nz+3

            x = 0.5*dx + DBLE(i - 1)*dx
            y = 0.5*dy + DBLE(j - 1)*dy
            z = 0.5*dz + DBLE(k - 1)*dz

            Q(i,j,k,:) = SIN(2.0*pi*x)*SIN(2.0*pi*y)*SIN(2.0*pi*z)
          END DO
        END DO
      END DO

      !$OMP PARALLEL DO schedule(static)
      DO nt=0,Nthreads-1
          CALL GetRHSDiffusion(Q(:,:,:,nt), RQ(:,:,:,nt), dx, dy, dz, 1, Nx, 1, Ny, 1, Nz, order)
      END DO
      !$OMP END PARALLEL DO

      DO nt=1,Nthreads-1
        IF ( MAXVAL(ABS(RQ(:,:,:,nt) - RQ(:,:,:,nt-1)))>1e-12) THEN
          WRITE(*,*) 'ERROR: '
          STOP
        END IF
        print *,MAXVAL(ABS(RQ(:,:,:,nt) - RQ(:,:,:,nt-1)))
      END DO


      DEALLOCATE(Q)
      DEALLOCATE(RQ)

    END DO ! Nthreads
  END DO ! N_v
END DO ! order

WRITE(*,*) '[0] -- Successful: GetRHSDiffusion on multiple threads on identical input yields identical output.'
END PROGRAM diffusion_threads_test