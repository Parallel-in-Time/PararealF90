# Puts a prefix in front of all timing files: This allows to quickly rename files from one scaling run
# for a later averaging of results in the plot routine
for f in timings_mpi*.dat;    do mv "$f" "`echo $f | sed s/timings/$1_timings/`"; done
for f in timings_openmp*.dat; do mv "$f" "`echo $f | sed s/timings/$1_timings/`"; done
for f in timings_serial*.dat; do mv "$f" "`echo $f | sed s/timings/$1_timings/`"; done
for f in openmp_*Np*.rur;     do mv "$f" "`echo $f | sed s/openmp/$1_openmp/`"; done
for f in mpi_Np*.rur;         do mv "$f" "`echo $f | sed s/mpi/$1_mpi/`"; done
for f in serial*.rur;         do mv "$f" "`echo $f | sed s/serial/$1_serial/`"; done
