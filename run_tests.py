import glob, os
tests = glob.glob('test/bin/*.out')
for file in tests:
    os.system(file)