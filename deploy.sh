#!/bin/bash
echo "purge garbage"
rm -f *.dat *.o *.x *.mod

echo "generate objects and executables"
make parameters.o
make heated_plate_mpi_subroutines.o
make heated_plate_mpi.o
make heated_plate_mpi.x

echo "running"
mpirun -np 36 ./heated_plate_mpi.x