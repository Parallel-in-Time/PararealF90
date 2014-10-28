include ../makefile.defs

all: objects run_timestepper.out run_parareal_mpi.out parareal_openmp.out parareal_openmp_pipe.out

DEP = ../obj/fluxes.o ../obj/upwind.o ../obj/weno5.o ../obj/advection.o ../obj/diffusion.o ../obj/boundaries.o ../obj/timestepper.o ../obj/parareal_mpi.o

objects:
	cd $(CURDIR)/../obj; make

run_timestepper.out: ../obj/run_timestepper.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/run_timestepper.o -o run_timestepper.out

run_parareal_mpi.out: ../obj/run_parareal_mpi.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/parareal_mpi.o ../obj/run_parareal_mpi.o -o parareal_mpi.out

parareal_openmp.out: ../obj/parareal_openmp.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/run_parareal_openmp.o -o run_parareal_openmp.out

parareal_openmp_pipe.out: ../obj/parareal_openmp_pipe.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/parareal_openmp_pipe.o -o parareal_openmp_pipe.out

clean:
	rm -rf *.out