include ../makefile.defs

INC=-I$(CURDIR)/..

all: $(patsubst ../src/%.f90, %.o, $(wildcard ../src/*.f90))

%.o: ../src/%.f90 ../preprocessor.f90
	$(F90) $(FCFLAGS) $(INC) -c $< -o $@

# NOTE: To ensure that *.mod files exist when required, all modules that USE other modules
# must have their dependencies list here explicitly.

advection.o: fluxes.o upwind.o weno5.o ../preprocessor.f90 ../src/advection.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/advection.f90 -o advection.o

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
