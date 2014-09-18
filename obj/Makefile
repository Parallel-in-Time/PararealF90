include ../makefile.defs

INC=-I$(CURDIR)/..

all: boundaries.o diffusion.o advection.o timestepper.o parareal_mpi.o parareal_openmp.o run_timestepper.o parareal_openmp_pipe.o

fluxes.o: ../src/fluxes.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/fluxes.f90 -o fluxes.o

upwind.o: ../src/upwind.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/upwind.f90 -o upwind.o

weno5.o: ../src/weno5.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/weno5.f90 -o weno5.o

diffusion.o: fluxes.o ../src/diffusion.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/diffusion.f90 -o diffusion.o

advection.o: fluxes.o upwind.o weno5.o ../preprocessor.f90 ../src/advection.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/advection.f90 -o advection.o

boundaries.o: ../src/boundaries.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/boundaries.f90 -o boundaries.o

timestepper.o: advection.o diffusion.o boundaries.o ../src/timestepper.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/timestepper.f90 -o timestepper.o

parareal_mpi.o: timestepper.o ../src/parareal_mpi.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/parareal_mpi.f90 -o parareal_mpi.o

parareal_openmp.o: timestepper.o ../src/parareal_openmp.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/parareal_openmp.f90 -o parareal_openmp.o

parareal_openmp_pipe.o: timestepper.o ../src/parareal_openmp_pipe.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/parareal_openmp_pipe.f90 -o parareal_openmp_pipe.o

run_timestepper.o : timestepper.o ../src/run_timestepper.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/run_timestepper.f90 -o run_timestepper.o

clean:
	rm -rf *.o
	rm -rf *.mod