import os
def build_runscript(Nproc,jobname,type,system):
    with open("submit_"+type+"_Np"+str(Nproc)+".sh", "w") as myfile:

        myfile.write("#!/bin/bash \n")
        myfile.write("#SBATCH --job-name="+jobname+"\n")
        myfile.write("#SBATCH --nodes=1")
        if type in ("mpi", "serial_f", "serial_g", "serial_f_ref"):
            myfile.write("#SBATCH --ntasks="+str(Nproc)+"\n")
            myfile.write("#SBATCH --ntasks-per-node="+str(Nproc)+"\n")
            #myfile.write("#SBATCH --cpus-per-task=1 \n")
        elif type in ("openmp","openmp_pipe"):
            myfile.write("#SBATCH --ntasks=1 \n")
            myfile.write("#SBATCH --ntasks-per-node=1 \n")
            #myfile.write("#SBATCH --cpus-per-task="+str(Nproc)+"\n")
        else:
            print "build_runscript: Encountered unknown type: "+type
        
        myfile.write("#SBATCH --time=00:05:00 \n")
        myfile.write("#SBATCH --mail-user=daniel.ruprecht@usi.ch \n")
        myfile.write("#SBATCH --output="+jobname+".out \n")
        myfile.write("echo JobID $SLURM_JOB_ID \n")
        myfile.write("echo Nodes $SLURM_NNODES \n")
        myfile.write("echo Tasks $SLURM_NTASKS \n")

        cwd = os.getcwd()
        
        if system=="cub":
            if type=="mpi":
                myfile.write('OMP_NUM_THREADS=1 mpirun -n '+str(Nproc)+' '+cwd+'/bin/run_parareal_mpi.out \n')
            elif type=="openmp":
                myfile.write('OMP_NUM_THREADS='+str(Nproc)+' mpirun -n 1 '+cwd+'/bin/run_parareal_openmp.out \n')
            elif type=="openmp_pipe":
                myfile.write('OMP_NUM_THREADS='+str(Nproc)+' mpirun -n 1 '+cwd+'/bin/run_parareal_openmp_pipe.out \n')
            elif type in ("serial_f","serial_f_ref"):
              myfile.write('mpirun -n 1 '+cwd+'/bin/run_timestepper.out F \n')
            elif type=="serial_g":
              myfile.write('mpirun -n 1 '+cwd+'/bin/run_timestepper.out C \n')
    
        elif system=="rosa":
            if type=="mpi":
                myfile.write('OMP_NUM_THREADS=1 aprun -n '+str(Nproc)+' -d 1 '+cwd+'/bin/run_parareal_mpi.out \n')
            elif type=="openmp":
                myfile.write('OMP_NUM_THREADS='+str(Nproc)+' aprun -n 1 -d '+str(Nproc)+' -cc=0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30 '+cwd+'/bin/run_parareal_openmp.out \n')
            elif type=="openmp_pipe":
                myfile.write('OMP_NUM_THREADS='+str(Nproc)+' aprun -n 1 -d '+str(Nproc)+' -cc=0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30 '+cwd+'/bin/run_parareal_openmp_pipe.out \n')
            elif type in ("serial_f","serial_f_ref"):
                myfile.write('aprun -n 1 '+cwd+'/bin/run_timestepper.out F \n')
            elif type=="serial_g":
                myfile.write('aprun -n 1 '+cwd+'/bin/run_timestepper.out C \n')

        elif system=="daint":
            if type=="mpi":
                myfile.write('OMP_NUM_THREADS=1 aprun -n '+str(Nproc)+' -d 1 '+cwd+'/bin/run_parareal_mpi.out \n')
            elif type=="openmp":
                myfile.write('OMP_NUM_THREADS='+str(Nproc)+' aprun -n 1 -d '+str(Nproc)+' '+cwd+'/bin/run_parareal_openmp.out \n')
            elif type=="openmp_pipe":
                myfile.write('OMP_NUM_THREADS='+str(Nproc)+' aprun -n 1 -d '+str(Nproc)+' '+cwd+'/bin/run_parareal_openmp_pipe.out \n')
            elif type in ("serial_f","serial_f_ref"):
                myfile.write('aprun -n 1 '+cwd+'/bin/run_timestepper.out F \n')
            elif type=="serial_g":
                myfile.write('aprun -n 1 '+cwd+'/bin/run_timestepper.out C \n')

        elif system=="mac":
            print "No SLURM on Mac, no runscript needed..."

        else:
            print "build_runscript: Encountered unknown string for system: "+system
