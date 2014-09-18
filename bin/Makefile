include ../makefile.defs

all: objects run_timestepper.out parareal_mpi.out parareal_openmp.out parareal_openmp_pipe.out

DEP = ../obj/fluxes.o ../obj/upwind.o ../obj/weno5.o ../obj/advection.o ../obj/diffusion.o ../obj/boundaries.o ../obj/timestepper.o 

objects:
	cd $(CURDIR)/../obj; make

run_timestepper.out: ../obj/run_timestepper.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/run_timestepper.o -o run_timestepper.out

parareal_mpi.out: ../obj/parareal_mpi.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/parareal_mpi.o -o parareal_mpi.out

parareal_openmp.out: ../obj/parareal_openmp.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/parareal_openmp.o -o parareal_openmp.out

parareal_openmp_pipe.out: ../obj/parareal_openmp_pipe.o
	$(F90) $(FCFLAGS) $(DEP) ../obj/parareal_openmp_pipe.o -o parareal_openmp_pipe.out

clean:
	rm -rf *.out