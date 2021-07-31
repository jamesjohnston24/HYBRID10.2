PROGRAM PROCESS

USE mpi
use netCDF

IMPLICIT NONE

INTEGER, PARAMETER :: nland = 67420
INTEGER :: nprocs, error, myrank, nland_chunk, file_handle
REAL, ALLOCATABLE, DIMENSION (:) :: B_k
CHARACTER(LEN=20) :: var_name
CHARACTER(LEN=200) :: file_name

!----------------------------------------------------------------------!
CALL MPI_INIT ( error )
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_Comm_size (MPI_COMM_WORLD,nprocs,error)
CALL MPI_Comm_rank (MPI_COMM_WORLD,myrank,error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
nland_chunk = nland / nprocs
kyr_clm = 1901
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
ALLOCATE (B_k(nland_chunk))
B_k = 20.0
var_name = 'B'
WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
&adf10/TRENDY2021/output/HYBRID10.3_",nprocs,&
&"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
WRITE (*,*) 'Reading from ', TRIM(file_name)
! Open the file for reading.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_read(file_handle, B_k, nland_chunk, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

write (*,*) myrank, B_k (10), error

!----------------------------------------------------------------------!
CALL MPI_FINALIZE ( error )
!----------------------------------------------------------------------!

END PROGRAM PROCESS