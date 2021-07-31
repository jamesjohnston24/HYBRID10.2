PROGRAM PROCESS

USE mpi

IMPLICIT NONE

INTEGER, PARAMETER :: root = 0
INTEGER :: nland = 67420, nlon = 720, nlat = 360
INTEGER :: myrank, nprocs, size, file_handle, kyr_clm, nland_chunk
INTEGER :: error, k
CHARACTER(LEN=200) :: var_name, file_name
REAL :: TB, TLA, summary
REAL, PARAMETER :: fillvalue = 1.0E20
REAL, ALLOCATABLE, DIMENSION (:) :: B_k, larea_k
REAL, ALLOCATABLE, DIMENSION (:,:) :: B_grid
INTEGER, ALLOCATABLE, DIMENSION (:) :: i_k, j_k

!----------------------------------------------------------------------!
CALL MPI_INIT ( error )
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
kyr_clm = 1901
nprocs = 4
myrank = 0
nland_chunk = nland / nprocs
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
var_name = 'larea'
ALLOCATE (larea_k (nland_chunk))
WRITE (file_name, "(A,I0.4,A,A,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
 &adf10/TRENDY2021/input/LUH2_GCB_2021/static_",nprocs,&
 &"CPUs/",TRIM(var_name),"_",myrank,".bin"
WRITE (*,*) 'Reading from ', TRIM(file_name)
! Open the file for reading.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_read(file_handle, larea_k, nland_chunk, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
var_name = 'i'
ALLOCATE (i_k (nland_chunk))
WRITE (file_name, "(A,I0.4,A,A,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
 &adf10/TRENDY2021/input/LUH2_GCB_2021/static_",nprocs,&
 &"CPUs/",TRIM(var_name),"_",myrank,".bin"
WRITE (*,*) 'Reading from ', TRIM(file_name)
! Open the file for reading.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_read(file_handle, i_k, nland_chunk, &
 MPI_INTEGER, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
var_name = 'j'
ALLOCATE (j_k (nland_chunk))
WRITE (file_name, "(A,I0.4,A,A,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
 &adf10/TRENDY2021/input/LUH2_GCB_2021/static_",nprocs,&
 &"CPUs/",TRIM(var_name),"_",myrank,".bin"
WRITE (*,*) 'Reading from ', TRIM(file_name)
! Open the file for reading.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_read(file_handle, j_k, nland_chunk, &
 MPI_INTEGER, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
ALLOCATE (B_k(nland_chunk))
var_name = 'B'
WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
&adf10/TRENDY2021/output/HYBRID10.3_",nprocs,&
&"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
WRITE (*,*) 'Reading from ', TRIM(file_name)
! Open the file for reading.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_read(file_handle, B_k, size, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
ALLOCATE (B_grid(nlon,nlat))
B_grid = fillvalue
DO k = 1, nland_chunk
 i = i_k (k)
 j = j_k (k)
 B_grid (i,j) = B_k (k)
END DO ! k
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_FINALIZE ( error )
!----------------------------------------------------------------------!

END PROGRAM PROCESS