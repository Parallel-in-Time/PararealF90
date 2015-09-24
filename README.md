[![][image]][link]

A description of the code and performance results can be found [here](http://arxiv.org/abs/1509.06935).

Parareal F90 {#mainpage}
============

This is a lightweight standalone implementation of Parareal solving 3D Burger's equation using a forward Euler and a RK3SSP method as coarse and fine integrators. It contains three different implementations of Parareal, one based on MPI, one using OpenMP without pipelining and one using OpenMP with pipelining. All three versions compute the same result, the purpose of the code is to compare different implement strategies with respect to speedup, memory footprint and energy consumption.

Documentation
-------------

If you have Doxygen installed, you can run *doxygen Doxyfile* in the base directory to generate HTML documentation of the code.

How do I get set up?
--------------------


To build, you will need the following files in the base directory:
  - preprocessor.f90
  - makefile.defs

In preprocessor.f90, a flag is defined that switch between linear advection and nonlinear advection (only the latter case is really Burger's equation, the other is a linear advection-diffusion problem). 

In makefile.defs, the compiler is specified and compiler flags are set. The code needs no libraries, only a MPI compiler with -enable-threads and the correct flag to allow for OpenMP directives (e.g. -fopenmp for GCC). Just specify the compiler and the flags in makefile.defs and type make in the base directory.

To run, need file
  - system.defs

The entries in system.defs are used to automatically build scripts to run the code on different platforms (see e.g. build_runscript.py in /scripts). There are several examples provided. If you want to include your own system, you will need to modify the following python scripts
  - get_run_cmd.py
  - build_runscript.py
  - run_parareal.py 
  - ...

Python runscripts
-----------------
There is a collection of python scripts that automatically run the code with different configurations to generate different sets of data. The script run_parareal_scaling.py for example runs Parareal with different numbers of processors plus the fine integrator as reference and generates all data necessary to plot speedup.

Test harness
------------
There is a set of tests, some written in Fortran and some in Python. Just type *python run_tests* in the base directory to start testing. If you add the flag *nof90* then only the high level python test scripts will be run.


Who do I talk to?
-----------------

This code is written by Daniel Ruprecht.

[image]:  https://zenodo.org/badge/doi/10.5281/zenodo.31288.svg
[link]:  http://dx.doi.org/10.5281/zenodo.31288