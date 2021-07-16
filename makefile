F90 = mpif90 -I/usr/local/Cluster-Apps/intel/2017.4/compilers_and_libraries_2017.4.196/linux/mpi/intel64/include

#FFLAGS = -check -traceback
FFLAGS = -O2
LDFLAGS = -lnetcdff -lnetcdf

SRC = \
	HYBRID10.2.f90 \
	RSF_In.f90 \
	get_clm.f90 \
	advance.f90 \
	Diag_Global.f90 \
	double.f90 \
	shared.f90

OBJ = $(SRC:.f90=.o)

HYBRID10.2.exe : $(OBJ)
	$(F90) $(FFLAGS) -o HYBRID10.2.exe $(OBJ) $(LDFLAGS)

# Main routine.
HYBRID10.2.o : RSF_In.o get_clm.o advance.o Diag_Global.o double.o shared.o HYBRID10.2.f90
	$(F90) $(FFLAGS) -c HYBRID10.2.f90
	
# Subroutines.
RSF_In.o : shared.o RSF_In.f90
	$(F90) $(FFLAGS) -c RSF_In.f90
	
get_clm.o : double.o shared.o get_clm.f90
	$(F90) $(FFLAGS) -c get_clm.f90
	
advance.o : double.o shared.o advance.f90
	$(F90) $(FFLAGS) -c advance.f90
	
Diag_Global.o : double.o shared.o Diag_Global.f90
	$(F90) $(FFLAGS) -c Diag_Global.f90
	
# Modules.
double.o : double.f90
	$(F90) $(FFLAGS) -c double.f90

shared.o : double.o shared.f90
	$(F90) $(FFLAGS) -c shared.f90
