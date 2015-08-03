def get_parameter():
  nu = 0.02
  Nx = 40
  Ny = 40
  Nz = 40
  dt_fine   = 1.0/240.0
  dt_coarse = 1.0/192.0
  Niter = 4
  Tend  = 1.0
  do_io = False
  be_verbose = False
  return nu, Nx, Ny, Nz, dt_fine, dt_coarse, Niter, Tend, do_io, be_verbose