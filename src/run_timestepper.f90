!>
!! @todo docu
!!

!> Module: run_timestepper
!> \author Daniel Ruprecht
!> \date 28 October, 2014
PROGRAM run_timestepper

USE omp_lib, only : OMP_GET_WTIME
USE timestepper, only : Euler, Rk3Ssp, InitializeTimestepper, FinalizeTimestepper

IMPLICIT NONE

INCLUDE 'mpif.h'

DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932_8 ! underscore indicates rounding to real(8) precision

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Q

CHARACTER(len=32) :: arg
DOUBLE PRECISION :: nu, Tend, dx, dy, dz, T0, T1
INTEGER :: Nx, Ny, Nz, N_fine, N_coarse, Niter, method, order_adv, order_diff, &
    ierr, mpi_thread_provided, Nthreads
LOGICAL :: do_io, be_verbose

CHARACTER(len=64) :: filename

!>
NAMELIST /param/ nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose

! ***********************
! Here starts the program
! ***********************

Nthreads=1

! Read parameter
OPEN(UNIT=20, FILE='parameter.in', ACTION='read', STATUS='old')
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

! Initialize
CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED, mpi_thread_provided, ierr)
CALL InitializeTimestepper(nu, Nx, Ny, Nz, Nthreads)

! Load initial value
ALLOCATE(Q(-2:Nx+3,-2:Ny+3,-2:Nz+3))
OPEN(unit=100, file='q0.dat', ACTION='read', STATUS='old')
READ(100,'(F35.25)') Q
CLOSE(100)

! Output
IF (be_verbose) THEN
    WRITE(*,'(A)') '--- Running serial timestepper ...'
    WRITE(*,'(A, F9.5)') 'Fine step length:   ', Tend/DBLE(N_fine)
    WRITE(*,'(A, F9.5)') 'Coarse step length: ', Tend/DBLE(N_coarse)
END IF

T0 = OMP_GET_WTIME()
IF (method==1) THEN
    CALL Euler(Q, DBLE(0.0), Tend, N_coarse, dx, dy, dz, order_adv, order_diff)
ELSE IF (method==3) THEN
    CALL Rk3Ssp(Q, DBLE(0.0), Tend, N_fine, dx, dy, dz, order_adv, order_diff)
END IF
T1 = OMP_GET_WTIME()

IF (do_io) THEN
  IF (arg=='C') THEN
    OPEN(UNIT=10, FILE='q_final_coarse.dat', ACTION='write', STATUS='replace')
  ELSE IF (arg=='F') THEN
    OPEN(UNIT=10, FILE='q_final_fine.dat', ACTION='write', STATUS='replace')
  END IF
  WRITE(10, '(F35.25)') Q(1:Nx,1:Ny,1:Nz)
  CLOSE(10)
END IF

WRITE(filename, '(A)') 'timings_serial_fine.dat'
OPEN(UNIT=10, FILE=filename, ACTION='write', STATUS='replace')
WRITE(10, '(F8.2)') T1-T0

! Finalize
CALL FinalizeTimestepper
CALL MPI_FINALIZE(ierr)

END PROGRAM run_timestepper