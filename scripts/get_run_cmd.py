import sys
def get_run_cmd(system, runcmd, mpi, Np):
  if system=="mac":
    if mpi:
      return "OMP_NUM_THREADS=1 "+runcmd+" -n "+str(Np)
    else:
      return "OMP_NUM_THREADS="+str(Np)+" "+runcmd+" -n 1"
  elif system=="cub":
    if mpi:
      return "OMP_NUM_THREADS=1 "+runcmd+" -n "+str(Np)
    else:
      return "OMP_NUM_THREADS="+str(Np)+" "+runcmd+" -n 1"
    return ""
  elif system=="rosa":
    if mpi:
      return "OMP_NUM_THREADS=1 "+runcmd+" -n "+str(Np)+" -d 1 "
    else:
      return "OMP_NUM_THREADS="+str(Np)+" "+runcmd+" -n "+str(Np)+" -d "+str(Np)+" "
  else:
    sys.exit("Unknow system: Check value provided in system.txt and make sure a configuration is provided in get_run_cmd.py")