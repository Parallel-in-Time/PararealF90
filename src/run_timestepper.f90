!>
!! @todo docu
!!

!> Module: runt_timestepper
!> \author Daniel Ruprecht
!> \date 28 October, 2014
PROGRAM run_timestepper

USE omp_lib
USE timestepper, only : Euler, Rk3Ssp, InitializeTimestepper, FinalizeTimestepper

IMPLICIT NONE

INCLUDE 'mpif.h'

DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932_8 ! underscore indicates rounding to real(8) precision

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q

CHARACTER(len=32) :: arg
DOUBLE PRECISION :: nu, Tend, dx, dy, dz
INTEGER :: Nx, Ny, Nz, N_fine, N_coarse, Niter, method, order_adv, order_diff, &
    ierr, mpi_thread_provided, Nthreads
LOGICAL :: do_io, be_verbose

!>
NAMELIST /param/ nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose

! ***********************
! Here starts the program
! ***********************

Nthreads=1

OPEN(unit=20, FILE='parameter.in', ACTION='read', STATUS='old')
READ(20,NML=param)
CLOSE(20)

! Command line arguments: C = run coarse method, F = run fine method
CALL GET_COMMAND_ARGUMENT(1, arg)
IF (arg=='C') THEN
    method     = 1
    order_adv  = 1
    order_diff = 2
ELSE IF (arg=='F') THEN
    method     = 3
    order_adv  = 5
    order_diff = 4
ELSE
    WRITE(*,*) 'First command line argument has to be either C or F. Found other value, now exiting.'
    STOP
END IF

! Computational domain is unit cube
dy = 1.0/DBLE(Ny)
dz = 1.0/DBLE(Nz)
dx = 1.0/DBLE(Nx)

CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)

CALL InitializeTimestepper(nu, Nx, Ny, Nz, Nthreads)
ALLOCATE(Q(-2:Nx+3,-2:Ny+3,-2:Nz+3))

OPEN(unit=100, file='q0.dat', ACTION='read', STATUS='old')
READ(100,'(F35.25)') Q
CLOSE(100)

IF (be_verbose) THEN
    WRITE(*,'(A)') '--- Running serial timestepper ...'
    WRITE(*,'(A, F9.5)') 'Fine step length:   ', Tend/DBLE(N_fine)
    WRITE(*,'(A, F9.5)') 'Coarse step length: ', Tend/DBLE(N_coarse)
END IF
    
IF (method==1) THEN
    CALL Euler(Q, DBLE(0.0), Tend, N_coarse, dx, dy, dz, order_adv, order_diff)
ELSE IF (method==3) THEN
    CALL Rk3Ssp(Q, DBLE(0.0), Tend, N_fine, dx, dy, dz, order_adv, order_diff)
END IF

OPEN(unit=100, file='qend.dat')
WRITE(100, '(F35.25)') Q(1:Nx,1:Ny,1:Nz)
CLOSE(100)

CALL FinalizeTimestepper

CALL MPI_FINALIZE(ierr)

END PROGRAM run_timestepper