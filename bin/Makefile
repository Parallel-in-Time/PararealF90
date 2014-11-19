include ../makefile.defs
include dependencies

all: objects run_timestepper.out run_parareal_mpi.out run_parareal_openmp.out parareal_openmp_pipe.out

# Update objects where necessary
objects:
	cd $(CURDIR)/../obj; make

run_timestepper.out: ../obj/run_timestepper.o ../obj/timestepper.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/run_timestepper.o -o run_timestepper.out

run_parareal_mpi.out: ../obj/run_parareal_mpi.o ../obj/parareal_mpi.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/run_parareal_mpi.o -o run_parareal_mpi.out

run_parareal_openmp.out: ../obj/run_parareal_openmp.o ../obj/parareal_openmp.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/run_parareal_openmp.o -o run_parareal_openmp.out

parareal_openmp_pipe.out: ../obj/parareal_openmp_pipe.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/parareal_openmp_pipe.o -o parareal_openmp_pipe.out

clean:
	rm -f *.out