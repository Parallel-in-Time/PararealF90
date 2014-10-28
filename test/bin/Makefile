include ../../makefile.defs
include dependencies

INC=-I$(CURDIR)/../obj -I$(CURDIR)/../../obj -I$(CURDIR)/../..

all: objects $(patsubst ../src/%.f90, %.o, $(wildcard ../src/*.f90))

objects:
	cd $(CURDIR)/../../obj; make

%.o: ../src/%.f90
	$(F90) $(FCFLAGS) $(INC) $(DEP) $< -o $@

clean:
	rm -rf *.out