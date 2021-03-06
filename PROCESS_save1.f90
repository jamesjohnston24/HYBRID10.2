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
WRITE (file_name, "(A,I0.4,A,A,A,I0.4,A)") "/home/jhj34/rds/rds-mb425-geogscratch/&
 &adf10/TRENDY2021/input/LUH2_GCB_2021/static_",nprocs,&
 &"CPUs/",TRIM(var_name),"_",myrank,".bin"
OPEN (10,FILE=file_name,STATUS='OLD',FORM='UNFORMATTED')
READ (10) larea_k
CLOSE (10)

write (*,*) larea_k

stop

!----------------------------------------------------------------------!
CALL MPI_INIT ( error )
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_Comm_size (MPI_COMM_WORLD,nprocs,error)
CALL MPI_Comm_rank (MPI_COMM_WORLD,myrank,error)
!----------------------------------------------------------------------!

size = nland / nprocs
kyr_clm = 1901

var_name = 'larea'
 WRITE (file_name, "(A,I0.4,A,A,A,I0.4,A)") "/home/jhj34/rds/rds-mb425-geogscratch/&
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

TLA = SUM (larea_k)
WRITE (*,*) 'TLA = ',TLA

ALLOCATE (B_k(nland_chunk))
var_name = 'B'
WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") "/home/jhj34/rds/rds-mb425-geogscratch/&
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

TB = 0.0
DO k = 1, nland_chunk
 ! (kg[DM] m^2) * (km^2) * (m^2 km^-2) * (g kg-1)
 TB = TB + B_k (k) * larea_k (k) * 1000.0 * 1000.0 * 1000.0
 !TB = TB + larea_k (k)
END DO ! k
WRITE (*,*) 'TB (Pg[DM]) = ', TB/1.0E15, myrank

! Combine and produce netCDF output for mapping.
CALL MPI_Reduce (TB, summary, 1, MPI_REAL, &
 MPI_SUM, root, MPI_COMM_WORLD, error)
!CALL MPI_Gather(B,nland_chunk,MPI_REAL, &
!                B_fin,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
IF (myrank == root) THEN
 write (*,*) 'Total biomass = ', summary/1.0e15, 'Pg[DM]'
 !file_name = 
END IF

!----------------------------------------------------------------------!
CALL MPI_FINALIZE ( error )
!----------------------------------------------------------------------!

END PROGRAM PROCESS