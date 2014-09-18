def build_runscript(Nproc,jobname,type,system):
    with open("submit_"+type+"_Np"+str(Nproc)+".sh", "w") as myfile:
        myfile.write("#!/bin/bash \n")
        myfile.write("#SBATCH --job-name="+jobname+"\n")
        if type=="mpi":
            myfile.write("#SBATCH --ntasks="+str(Nproc)+"\n")
            myfile.write("#SBATCH --ntasks-per-node="+str(Nproc)+"\n")
            #myfile.write("#SBATCH --cpus-per-task=1 \n")
        elif type in ("openmp","openmp_pipe"):
            myfile.write("#SBATCH --ntasks=1 \n")
            myfile.write("#SBATCH --ntasks-per-node=1 \n")
            #myfile.write("#SBATCH --cpus-per-task="+str(Nproc)+"\n")
        else:
            print "Encountered unknown type: "+type
        myfile.write("#SBATCH --time=00:05:00 \n")
        myfile.write("#SBATCH --mail-user=daniel.ruprecht@usi.ch \n")
        myfile.write("#SBATCH --output="+jobname+".out \n")
        myfile.write("echo JobID $SLURM_JOB_ID \n")
        myfile.write("echo Nodes $SLURM_NNODES \n")
        myfile.write("echo Tasks $SLURM_NTASKS \n")
        if system=="cub":
            if type=="mpi":
                myfile.write('OMP_NUM_THREADS=1 mpirun -n '+str(Nproc)+' bin/parareal_mpi.out \n')
            elif type=="openmp":
                myfile.write("OMP_NUM_THREADS="+str(Nproc)+" mpirun -n 1 "+run_add+" bin/parareal_openmp.out \n")
            elif type=="openmp_pipe":
                myfile.write("OMP_NUM_THREADS="+str(Nproc)+" mpirun -n 1 "+run_add+" bin/parareal_openmp_pipe.out \n")   
        elif system=="rosa":
            if type=="mpi":
                myfile.write('OMP_NUM_THREADS=1 aprun -n '+str(Nproc)+' bin/parareal_mpi.out \n')
            elif type=="openmp":
                myfile.write("OMP_NUM_THREADS="+str(Nproc)+" aprun -n 1 -d "+str(Nproc)+" bin/parareal_openmp.out \n")
            elif type=="openmp_pipe":
                myfile.write("OMP_NUM_THREADS="+str(Nproc)+" aprun -n 1 -d "+str(Nproc)+" bin/parareal_openmp_pipe.out \n")           
        elif system=="mac":
            print "No SLURM on Mac, no runscript needed..."          
        else:
            print "Encountered unknown string for system: "+system     