include makefile.defs

all: compile

compile:
	cd $(CURDIR)/obj; make
	cd $(CURDIR)/test/obj; make
	cd $(CURDIR)/bin; make
	cd $(CURDIR)/test/bin; make
	cd $(CURDIR)/scaling/obj; make
	cd $(CURDIR)/scaling/bin; make

clean:
	cd $(CURDIR)/obj; make clean
	cd $(CURDIR)/test/obj; make clean
	cd $(CURDIR)/bin; make clean
	cd $(CURDIR)/test/bin; make clean
	cd $(CURDIR)/test/scripts; make clean
	cd $(CURDIR)/scaling/obj; make clean
	cd $(CURDIR)/scaling/bin; make clean
	rm -f *.dat
	rm -f *.pyc
	rm -f scripts/*.pyc
	rm -f *.out
	rm -f submit_*.sh

cleandata:
	rm -f *.dat
	rm -f *.out
	rm -f *.sh
	rm -f *.in
	rm -f rur.*