F90 = mpif90 -I/usr/local/Cluster-Apps/intel/2017.4/compilers_and_libraries_2017.4.196/linux/mpi/intel64/include

FFLAGS = -check -traceback
#FFLAGS = -O2
LDFLAGS = -lnetcdff -lnetcdf

SRC = \
	make_inputs.f90 \
	double.f90

OBJ = $(SRC:.f90=.o)

make_inputs.exe : $(OBJ)
	$(F90) $(FFLAGS) -o make_inputs.exe $(OBJ) $(LDFLAGS)

# Main routine.
make_inputs.o : double.o make_inputs.f90
	$(F90) $(FFLAGS) -c make_inputs.f90
	
# Subroutines.
	
# Modules.
double.o : double.f90
	$(F90) $(FFLAGS) -c double.f90
