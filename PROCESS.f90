PROGRAM PROCESS

USE mpi
use netCDF

IMPLICIT NONE

INTEGER, PARAMETER :: root = 0
INTEGER, PARAMETER :: nland = 67420, nlon = 720, nlat = 360
INTEGER :: myrank, nprocs, size, file_handle, kyr_clm, nland_chunk
INTEGER :: error, i, j, k
INTEGER :: lon_dimid, lat_dimid, lon_varid, lat_varid, ncid, varid
INTEGER :: varidW, varidB, varidSOM
INTEGER, DIMENSION (2) :: dimids_two
CHARACTER(LEN=300) :: var_name, file_name
CHARACTER(LEN=4) :: char_year
REAL :: TB, TLA, summary
REAL, PARAMETER :: fillvalue = 1.0E20
REAL, ALLOCATABLE, DIMENSION (:) :: B_k, larea_k
REAL, ALLOCATABLE, DIMENSION (:,:) :: soilW_grid
REAL, ALLOCATABLE, DIMENSION (:,:) :: B_grid
REAL, ALLOCATABLE, DIMENSION (:,:) :: SOM_grid
INTEGER, ALLOCATABLE, DIMENSION (:) :: i_k, j_k
REAL, DIMENSION (nlon) :: lon
REAL, DIMENSION (nlat) :: lat

!----------------------------------------------------------------------!
CALL MPI_INIT ( error )
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_Comm_size (MPI_COMM_WORLD,nprocs,error)
CALL MPI_Comm_rank (MPI_COMM_WORLD,myrank,error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
kyr_clm = 1901
nland_chunk = nland / nprocs
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
CALL MPI_File_read(file_handle, B_k, size, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

write (*,*) myrank, B_k (10), error

!----------------------------------------------------------------------!
CALL MPI_FINALIZE ( error )
!----------------------------------------------------------------------!

END PROGRAM PROCESS