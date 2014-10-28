include ../makefile.defs

INC=-I$(CURDIR)/..

all: $(patsubst ../src/%.f90, %.o, $(wildcard ../src/*.f90))

%.o: ../src/%.f90
	$(F90) -c $(FCFLAGS) $(INC) $< -o $@

include dependencies

clean:
	rm -f *.o
	rm -f *.mod