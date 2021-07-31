PROGRAM PROCESS

USE mpi

IMPLICIT NONE

INTEGER, PARAMETER :: root = 0
INTEGER :: nland = 67420
INTEGER :: myrank, nprocs, size, file_handle, kyr_clm, nland_chunk
INTEGER :: error, k
CHARACTER(LEN=200) :: var_name, file_name
REAL :: TB, TLA, summary
REAL, ALLOCATABLE, DIMENSION (:) :: B_k, larea_k

var_name = 'larea'
nprocs = 4
myrank = 0
nland_chunk = nland / nprocs
ALLOCATE (larea_k (nland_chunk))
WRITE (file_name, "(A,I0.4,A,A,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
 &adf10/TRENDY2021/input/LUH2_GCB_2021/static_",nprocs,&
 &"CPUs/",TRIM(var_name),"_",myrank,".bin"
WRITE (*,*) 'Reading from ', TRIM(file_name)
! Open the file for reading.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_read(file_handle, larea_k, size, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)

write (*,*) larea_k

stop

END PROGRAM PROCESS