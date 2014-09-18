include ../../makefile.defs

INC=-I$(CURDIR)/../../obj -I$(CURDIR)/../../

all: advection_test.o diffusion_test.o advection_diffusion_test.o boundaries_test.o advection_boundaries_test.o timestepper_test.o

advection_test.o: ../../obj/advection.o ../src/advection_test.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/advection_test.f90 -o advection_test.o

diffusion_test.o: ../../obj/diffusion.o ../src/diffusion_test.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/diffusion_test.f90 -o diffusion_test.o

advection_diffusion_test.o: ../../obj/advection.o ../../obj/diffusion.o ../src/advection_diffusion_test.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/advection_diffusion_test.f90 -o advection_diffusion_test.o

boundaries_test.o: ../../obj/boundaries.o ../src/boundaries_test.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/boundaries_test.f90 -o boundaries_test.o

advection_boundaries_test.o: ../../obj/advection.o ../../obj/boundaries.o ../src/advection_boundaries_test.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/advection_boundaries_test.f90 -o advection_boundaries_test.o

timestepper_test.o: ../../obj/advection.o ../../obj/diffusion.o ../../obj/timestepper.o ../../obj/boundaries.o ../src/timestepper_test.f90
	$(F90) $(FCFLAGS) $(INC) -c ../src/timestepper_test.f90 -o timestepper_test.o

clean:
	rm -rf *.o