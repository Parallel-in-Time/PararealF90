Np = 2
Niter = 1

type = 'mpi'
nps     = '%0.2i' % Np
nps_m1  = '%0.2i' % (Np-1)
niter_s = '%0.2i' % Niter

fser = open('qend.dat')
fpar = open('q_final_'+niter_s+'_'+nps_m1+'_'+nps+'_'+type+'.dat')

max_err = 0.0
while True:
    line = fser.readline()
    if not line: break
    Qser = float(line)
    Qpar = float(fpar.readline())
    max_err = max( abs(Qser-Qpar), max_err )

print ('Maximum difference: %7.5e' % max_err)
