include ../../makefile.defs

INC:=-I$(CURDIR)/../../obj -I$(CURDIR)/../../

all: advection_weakscaling.o diffusion_weakscaling.o advection_weakscaling_mpi.o

advection_weakscaling.o: ../src/advection_weakscaling.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/advection_weakscaling.f90 -o advection_weakscaling.o

diffusion_weakscaling.o: ../src/diffusion_weakscaling.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/diffusion_weakscaling.f90 -o diffusion_weakscaling.o

advection_weakscaling_mpi.o: ../src/advection_weakscaling_mpi.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/advection_weakscaling_mpi.f90 -o advection_weakscaling_mpi.o

clean:
	rm -rf *.o
	rm -rf *.mod