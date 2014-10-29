import glob, os, sys
sys.path.append('./test/scripts')
sys.path.append('./scripts')
from para_compare_test import para_compare_test
from para_fine_test import para_fine_test
from para_coarse_test import para_coarse_test

para_coarse_test();
para_fine_test();
para_compare_test();

tests = glob.glob('test/bin/*.out')
for file in tests:
  os.system(file)
