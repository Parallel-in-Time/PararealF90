import glob, os, sys
sys.path.append('./test/scripts')
sys.path.append('./scripts')
sys.path.append('./scripts/termcolor-1.1.0')
from para_compare_test import para_compare_test
from para_fine_test import para_fine_test
from para_coarse_test import para_coarse_test
import multiprocessing

with open("system.txt", "r") as rfile:
    system = rfile.readline()
    system = system.rstrip()
    runcmd = rfile.readline()
    runcmd = runcmd.rstrip()
    rfile.close()

be_verbose = True
Ntests = 1

para_coarse_test(system,runcmd,be_verbose,Ntests);
os.system('rm -f q_final*.dat')
os.system('rm -f parameter.in')

para_fine_test(system,runcmd,be_verbose,Ntests);
os.system('rm -f q_final*.dat')
os.system('rm -f parameter.in')

para_compare_test(system,runcmd,be_verbose, 50);
os.system('rm -f q_final*.dat')
os.system('rm -f parameter.in')

Np    = multiprocessing.cpu_count()
tests = glob.glob('test/bin/*.out')
for file in tests:
  os.system('OMP_NUM_THREADS='+str(Np)+' '+runcmd+' -n 1 '+file)
