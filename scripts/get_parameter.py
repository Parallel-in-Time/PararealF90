def get_parameter():
  nu = 0.01
  Nx = 32
  Ny = 32
  Nz = 32
  dt_fine   = 1.0/96.0
  dt_coarse = 1.0/64.0
  Niter = 2
  Tend  = 4.0
  do_io = False
  be_verbose = True
  return nu, Nx, Ny, Nz, dt_fine, dt_coarse, Niter, Tend, do_io, be_verbose