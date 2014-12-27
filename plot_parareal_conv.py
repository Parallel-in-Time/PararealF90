import numpy as np
from matplotlib import pyplot as plt

fs = 18
Nproc = 8
type="mpi"
defect = np.zeros(Nproc-1)
y_fine = np.array([])

filename = "q_final_fine.dat"
with open(filename,'r') as fobj:
  while True:
    line = fobj.readline()
    if not line: break
    y_fine = np.append(y_fine, [float(line)])
fobj.close()
Nfine = np.size(y_fine)

for niter in range(1,Nproc):
  y_para = np.array([])
  niter_s    = '%0.2i' % niter
  nproc_s    = '%0.2i' % Nproc
  nproc_m1_s = '%0.2i' % (Nproc-1)
  filename = "q_final_"+niter_s+"_"+nproc_m1_s+"_"+nproc_s+"_"+type+".dat"
  with open(filename,'r') as fobj:
    while True:
      line = fobj.readline()
      if not line: break
      y_para = np.append(y_para, [float(line)])
  fobj.close()
  Npara = np.size(y_para)
  assert Npara == Nfine, 'Mismatch in length of serial and parallel solution.'
  diff = y_fine - y_para
  defect[niter-1] = np.linalg.norm(diff)

print defect