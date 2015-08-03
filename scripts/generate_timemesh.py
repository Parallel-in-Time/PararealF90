import math
#
#
def generate_timemesh(T0, Tend, dt_fine, dt_coarse, Np):
  dt_slice = (Tend-T0)/Np
  Nfine    = int(math.ceil(dt_slice/dt_fine))
  Ncoarse  = int(math.ceil(dt_slice/dt_coarse))
  dt_fine_actual = dt_slice/Nfine
  dt_coarse_actual = dt_slice/Ncoarse
  if abs(dt_fine - dt_fine_actual)/abs(dt_fine) > 1e-12:
      print "WARNING: Actual fine timestep deviates from requested value"
      print ("Requested fine timestep: %9.7f" % dt_fine)
      print ("Acutal fine timestep:    %9.7f" % dt_fine_actual)
  if abs(dt_coarse - dt_coarse_actual)/abs(dt_coarse) > 1e-12:
      print "WARNING: Actual coarse timestep deviates from requested value"
      print ("Requested fine timestep: %9.7f" % dt_coarse)
      print ("Acutal fine timestep:    %9.7f" % dt_coarse_actual)      
  return {'Nfine':Nfine, 'Ncoarse':Ncoarse}