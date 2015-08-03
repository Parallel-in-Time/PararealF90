!>
!! Module that reads in and provides parameters for the simulation.
!!
MODULE params

IMPLICIT NONE

!> If TRUE, the solution will be written out. Otherwise only timings are generated.
LOGICAL :: do_io

!> If TRUE, Parareal will generate a number of messages during runtime.
LOGICAL ::be_verbose

!> Diffusion coefficient in Burgers equation
DOUBLE PRECISION :: nu

!> Mesh width in x direction
DOUBLE PRECISION :: dx

!> Mesh width in y direction
DOUBLE PRECISION :: dy

!> Mesh width in z direction
DOUBLE PRECISION :: dz

!> Final simulation time
DOUBLE PRECISION :: Tend

!> Number of cells in x direction
INTEGER :: Nx

!> Number of cells in y direction
INTEGER :: Ny

!> Number of cells in z direction
INTEGER :: Nz

!> Time steps *per timeslice* for the fine integrator
INTEGER :: N_fine

!> Time steps *per timeslice* for the coarse integrator
INTEGER :: N_coarse

!> Number of Parareal iterations
INTEGER :: Niter


CHARACTER(len=32) :: param_file

NAMELIST /param/ nu, Nx, Ny, Nz, N_fine, N_coarse, Niter, Tend, do_io, be_verbose

CONTAINS

  !>
  !! Read parameter namelist from input file provided as first command line argument
  SUBROUTINE ReadParameter()

    CALL GET_COMMAND_ARGUMENT(1,param_file)

    ! Read parameters
    OPEN(unit=20, FILE=param_file, ACTION='read', STATUS='old')
    READ(20,NML=param)
    CLOSE(20)

    ! Computational domain is unit cube
    dy = DBLE(1.0)/DBLE(Ny)
    dz = DBLE(1.0)/DBLE(Nz)
    dx = DBLE(1.0)/DBLE(Nx)

  END SUBROUTINE ReadParameter

END MODULE params