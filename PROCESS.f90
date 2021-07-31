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
CHARACTER(LEN=200) :: var_name, file_name
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
kyr_clm = 1901
nprocs = 4
nland_chunk = nland / nprocs
ALLOCATE (soilW_grid(nlon,nlat))
ALLOCATE (B_grid(nlon,nlat))
ALLOCATE (SOM_grid(nlon,nlat))
ALLOCATE (larea_k (nland_chunk))
ALLOCATE (i_k (nland_chunk))
ALLOCATE (j_k (nland_chunk))
ALLOCATE (B_k (nland_chunk))
B_grid = fillvalue
myrank=0
!----------------------------------------------------------------------!

!DO myrank = 0, 0

!----------------------------------------------------------------------!
var_name = 'larea'
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
DO k = 1, nland_chunk
 i = i_k (k)
 j = j_k (k)
 B_grid (i,j) = B_k (k)
 write(*,*) i,j,k,B_grid(i,j)
END DO ! k
!----------------------------------------------------------------------!

!END DO ! myrank

!----------------------------------------------------------------------!
var_name = 'tmp'
kyr_clm =1901
WRITE (char_year, '(I4)') kyr_clm
file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/TRENDY2021/&
 &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
 &char_year//'.365d.noc.nc'
WRITE (*,*) 'Opening file: ',file_name
CALL CHECK ( NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid ))
varid = 2
CALL CHECK ( NF90_GET_VAR ( ncid, varid, lon ))
varid = 3
CALL CHECK ( NF90_GET_VAR ( ncid, varid, lat ))
CALL CHECK ( NF90_CLOSE ( ncid ))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
file_name = "fields_grid.nc"
write (*, *) 'Writing to ', trim (file_name)
CALL CHECK (NF90_CREATE (trim (file_name), cmode = nf90_clobber, &
            ncid = ncid))
CALL CHECK (NF90_DEF_DIM (ncid, "longitude", nlon, lon_dimid))
CALL CHECK (NF90_DEF_DIM (ncid, "latitude" , nlat, lat_dimid))
CALL CHECK (NF90_DEF_VAR (ncid, "longitude", nf90_float, lon_dimid, &
            lon_varid))
CALL CHECK (NF90_DEF_VAR (ncid, "latitude" , nf90_float, lat_dimid, &
            lat_varid))
dimids_two = (/ lon_dimid, lat_dimid /)
CALL CHECK (NF90_PUT_ATT (ncid, lon_varid, "units", "degrees_east"))
CALL CHECK (NF90_PUT_ATT (ncid, lat_varid, "units", "degrees_north"))
CALL CHECK (NF90_DEF_VAR (ncid, "Soil_Water", nf90_float, &
            dimids_two, varidW))
CALL CHECK (NF90_DEF_VAR (ncid, "Living_Biomass", nf90_float, &
            dimids_two, varidB))
CALL CHECK (NF90_DEF_VAR (ncid, "Soil_Organic_Matter", nf90_float, &
            dimids_two, varidSOM))
CALL CHECK (NF90_PUT_ATT (ncid, varidW, "units", "m"))
CALL CHECK (NF90_PUT_ATT (ncid, varidB, "units", "kg[DM] m-2"))
CALL CHECK (NF90_PUT_ATT (ncid, varidSOM, "units", "kg[DM] m-2"))
CALL CHECK (NF90_PUT_ATT (ncid, varidW, "_FillValue", fillvalue))
CALL CHECK (NF90_PUT_ATT (ncid, varidB, "_FillValue", fillvalue))
CALL CHECK (NF90_PUT_ATT (ncid, varidSOM, "_FillValue", fillvalue))
CALL CHECK (NF90_ENDDEF (ncid))
CALL CHECK (NF90_PUT_VAR (ncid, lon_varid, lon))
CALL CHECK (NF90_PUT_VAR (ncid, lat_varid, lat))
CALL CHECK (NF90_PUT_VAR (ncid,     varidW, soilW_grid))
CALL CHECK (NF90_PUT_VAR (ncid,     varidB, B_grid))
CALL CHECK (NF90_PUT_VAR (ncid,     varidSOM, SOM_grid))
CALL CHECK (NF90_close (ncid))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_FINALIZE ( error )
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CONTAINS
 SUBROUTINE check ( status )

 INTEGER, INTENT ( in ) :: status
 IF (status /= nf90_noerr) THEN
  PRINT *, TRIM (NF90_STRERROR( status ))
  STOP  "Stopped"
 END IF
 END SUBROUTINE check
!----------------------------------------------------------------------!

END PROGRAM PROCESS