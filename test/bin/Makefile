include ../../makefile.defs

INC=-I$(CURDIR)/../obj -I$(CURDIR)../../obj

#DEP:=$(wildcard ../../obj/*.o)
DEP:=../../obj/fluxes.o ../../obj/boundaries.o ../../obj/upwind.o ../../obj/weno5.o ../../obj/timestepper.o ../../obj/advection.o ../../obj/diffusion.o

all: boundaries_test.out advection_test.out diffusion_test.out advection_boundaries_test.out advection_diffusion_test.out timestepper_test.out

boundaries_test.out: ../obj/boundaries_test.o ../src/boundaries_test.f90
	$(F90) $(FCFLAGS) $(INC) $(DEP) ../obj/boundaries_test.o -o boundaries_test.out

advection_test.out: ../obj/advection_test.o ../src/advection_test.f90
	$(F90) $(FCFLAGS) $(INC) $(DEP) ../obj/advection_test.o -o advection_test.out

diffusion_test.out: ../obj/diffusion_test.o ../src/diffusion_test.f90
	$(F90) $(FCFLAGS) $(INC) $(DEP) ../obj/diffusion_test.o -o diffusion_test.out	

advection_boundaries_test.out: ../obj/advection_boundaries_test.o ../src/advection_boundaries_test.f90
	$(F90) $(FCFLAGS) $(INC) $(DEP) ../obj/advection_boundaries_test.o -o advection_boundaries_test.out

advection_diffusion_test.out: ../obj/advection_diffusion_test.o ../src/advection_diffusion_test.f90
	$(F90) $(FCFLAGS) $(INC) $(DEP) ../obj/advection_diffusion_test.o -o advection_diffusion_test.out

timestepper_test.out: ../obj/timestepper_test.o ../src/timestepper_test.f90
	$(F90) $(FCFLAGS) $(INC) $(DEP) ../obj/timestepper_test.o -o timestepper_test.out

clean:
	rm -rf *.out