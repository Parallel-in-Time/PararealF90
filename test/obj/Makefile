include ../../makefile.defs

INC=-I$(CURDIR)/../../obj -I$(CURDIR)/../../

all: objects $(patsubst ../src/%.f90, %.o, $(wildcard ../src/*.f90))

objects:
	cd $(CURDIR)/../../obj; make

%.o: ../src/%.f90
	$(F90) -c $(FCFLAGS) $(INC) $< -o $@

clean:
	rm -f *.o
	rm -f *.mod 