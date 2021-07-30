PROGRAM PROCESS
IMPLICIT NONE

INTEGER :: myrank, nprocs, size, file_handle, kyr_clm
CHARACTER(LEN=200) :: var_name, file_name
REAL, ALLOCATABLE. DIMENSION (:) :: B

nprocs = 4
size = nland / nprocs
kyr_clm = 1920
ALLOCATE (B(size))

DO myrank = 1, nprocs

 var_name = 'B'
 WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
 &adf10/TRENDY2021/output/HYBRID10.3_",nprocs,&
 &"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
 WRITE (*,*) 'Reading from ', TRIM(file_name)
 ! Open the file for reading.
 CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
  MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
 ! MPI_IO is binary output format. Write using individual file pointer.
 CALL MPI_File_read(file_handle, B, size, &
  MPI_REAL, MPI_STATUS_IGNORE, error)
 ! Close the file.
 CALL MPI_File_Close(file_handle, error)

END DO ! myrank = 1, nprocs

END PROGRAM PROCESS