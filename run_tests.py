import glob, os

#tests = glob.glob('test/scripts/*.py')
#for file in tests:
#  os.system('python '+file)

tests = glob.glob('test/bin/*.out')
for file in tests:
    os.system(file)
