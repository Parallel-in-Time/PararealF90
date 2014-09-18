include ../../makefile.defs

INC=-I$(CURDIR)/../obj -I$(CURDIR)../../obj

DEP:=../../obj/fluxes.o ../../obj/boundaries.o ../../obj/upwind.o ../../obj/weno5.o ../../obj/timestepper.o ../../obj/advection.o ../../obj/diffusion.o

all: advection_weakscaling.out diffusion_weakscaling.out advection_weakscaling_mpi.out

advection_weakscaling.out: ../obj/advection_weakscaling.o
	$(F90) $(FCFLAGS) $(INC) $(DEP) ../obj/advection_weakscaling.o -o advection_weakscaling.out

diffusion_weakscaling.out: ../obj/diffusion_weakscaling.o
	$(F90) $(FCFLAGS) $(INC) $(DEP) ../obj/diffusion_weakscaling.o -o diffusion_weakscaling.out

advection_weakscaling_mpi.out: ../obj/advection_weakscaling_mpi.o
	$(F90) $(FCFLAGS) $(INC) $(DEP) ../obj/advection_weakscaling_mpi.o -o advection_weakscaling_mpi.out

clean:
	rm -rf *.out