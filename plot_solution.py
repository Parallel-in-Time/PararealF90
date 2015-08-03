import sys
sys.path.append('./scripts')
from get_parameter import get_parameter

import numpy as np

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

nu, Nx, Ny, Nz, dt_fine, dt_coarse, Niter, Tend, do_io, be_verbose = get_parameter()
sol = np.array([])

filename = "q_final_fine.dat"
with open(filename,'r') as fobj:
  while True:
    line = fobj.readline()
    if not line: break
    sol = np.append(sol, [float(line)])

assert np.size(sol)==Nx*Ny*Nz, 'Length of solution does not match parameter... was probably generated with different setting.'

sol.shape = ((Nx, Ny, Nz))

x = np.linspace(0, 1, Nx)
y = np.linspace(0, 1, Ny)
xx, yy = np.meshgrid(x, y)

fig = plt.figure(figsize=(8,8))
ax = fig.gca(projection='3d')
ax.view_init(elev=0., azim=-90.)
surf = ax.plot_surface(xx, yy, sol[:,:,0], rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_xlim(left   = 0.0, right = 1.0)
ax.set_ylim(bottom = 0.0, top   = 1.0)
#ax.set_zlim(bottom = 0.0, top   = 1.0)
plt.xlabel('x')
plt.ylabel('y')
plt.show()
