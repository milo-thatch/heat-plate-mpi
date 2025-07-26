CC = mpif90
CCFLAGS = -ftree-vectorize

## objects
heated_plate_mpi.o: heated_plate_mpi.f90
	$(CC) $(CCFLAGS) -c heated_plate_mpi.f90 -o heated_plate_mpi.o
parameters.o: parameters.f90
	$(CC) $(CCFLAGS) -c parameters.f90 -o parameters.o

## executables
heated_plate_mpi.x: parameters.o heated_plate_mpi.o
	$(CC) $(CCFLAGS) parameters.o heated_plate_mpi.o -o heated_plate_mpi.x

