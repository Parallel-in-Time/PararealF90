import numpy as np

disc_err_fine   = 0.0
disc_err_coarse = 0.0
y_ref    = np.array([])
y_fine   = np.array([])
y_coarse = np.array([])

# load fine reference
filename = "q_final_fine_ref.dat"
with open(filename,'r') as fobj:
  while True:
    line = fobj.readline()
    if not line: break
    y_ref = np.append(y_ref, [float(line)])

# load fine serial
filename = "q_final_fine.dat"
with open(filename,'r') as fobj:
  while True:
    line = fobj.readline()
    if not line: break
    y_fine = np.append(y_fine, [float(line)])

# load coarse serial
filename = "q_final_coarse.dat"
with open(filename,'r') as fobj:
  while True:
    line = fobj.readline()
    if not line: break
    y_coarse = np.append(y_coarse, [float(line)])

assert np.size(y_ref)==np.size(y_fine), 'Mismatch in length of fine and reference solution'
assert np.size(y_fine)==np.size(y_coarse), 'Mismatch in length of fine and coarse solution'

err_fine   = y_fine - y_ref
err_coarse = y_coarse - y_ref

print ('Approximate discretization error of fine method:   %6.3e' % (np.linalg.norm(err_fine, np.inf)/np.linalg.norm(y_ref, np.inf)))
print ('Approximate discretization error of coarse method: %6.3e' % (np.linalg.norm(err_coarse, np.inf)/np.linalg.norm(y_ref, np.inf)))
