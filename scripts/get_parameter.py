def get_parameter():
  nu = 0.01
  Nx = 32
  Ny = 32
  Nz = 32
  dt_fine   = 1.0/200.0
  dt_coarse = 1.0/80.0
  Niter = 1
  Tend  = 2.0
  do_io = True
  be_verbose = False
  return nu, Nx, Ny, Nz, dt_fine, dt_coarse, Niter, Tend, do_io, be_verbose