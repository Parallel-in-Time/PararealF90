include ../makefile.defs

INC=-I$(CURDIR)/..

all: 
	$(F90) -c $(FCFLAGS) $(INC) ../src/fluxes.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/boundaries.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/params.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/upwind.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/weno5.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/advection.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/diffusion.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/timestepper.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/parareal_mpi.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/parareal_openmp.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/parareal_openmp_pipe.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/run_parareal_mpi.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/run_parareal_openmp.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/run_parareal_openmp_pipe.f90
	$(F90) -c $(FCFLAGS) $(INC) ../src/run_timestepper.f90

include dependencies

clean:
	rm -f *.o
	rm -f *.mod
