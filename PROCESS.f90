PROGRAM PROCESS

! 

USE mpi
use netCDF

IMPLICIT NONE

INTEGER, PARAMETER :: root = 0
INTEGER, PARAMETER :: nland = 67420, nlon = 720, nlat = 360
REAL, PARAMETER :: fillvalue = 1.0E20
INTEGER :: nprocs, error, myrank, nland_chunk, file_handle, kyr_clm
INTEGER :: i, j, k, nyr_spin, kyr_rsf
REAL :: TA, TB, TW, TSOM, TaNPP, TaRh, TaNBP
REAL :: Wmax, Bmax, SOMmax, aNPPmax, aRhmax, aNBPmax
REAL, ALLOCATABLE, DIMENSION (:) :: B_k, larea_k, B_k_all, larea_k_all
REAL, ALLOCATABLE, DIMENSION (:) :: soilW_k, soilW_k_all
REAL, ALLOCATABLE, DIMENSION (:) :: SOM_k, SOM_k_all
REAL, ALLOCATABLE, DIMENSION (:) :: aNPP_k, aNPP_k_all
REAL, ALLOCATABLE, DIMENSION (:) :: aRh_k, aRh_k_all
REAL, ALLOCATABLE, DIMENSION (:) :: aNBP_k, aNBP_k_all
REAL, ALLOCATABLE, DIMENSION (:,:) :: soilW_grid
REAL, ALLOCATABLE, DIMENSION (:,:) :: B_grid
REAL, ALLOCATABLE, DIMENSION (:,:) :: SOM_grid
REAL, ALLOCATABLE, DIMENSION (:,:) :: aNPP_grid
REAL, ALLOCATABLE, DIMENSION (:,:) :: aRh_grid
REAL, ALLOCATABLE, DIMENSION (:,:) :: aNBP_grid
INTEGER, ALLOCATABLE, DIMENSION (:) :: i_k, j_k, i_k_all, j_k_all
REAL, DIMENSION (nlon) :: lon
REAL, DIMENSION (nlat) :: lat
CHARACTER(LEN=20) :: var_name
CHARACTER(LEN=200) :: file_name
INTEGER :: lon_dimid, lat_dimid, lon_varid, lat_varid, ncid, varid
INTEGER :: varidW, varidB, varidSOM, varidaNPP, varidaRh, varidaNBP
INTEGER, DIMENSION (2) :: dimids_two
CHARACTER(LEN=4) :: char_year

!----------------------------------------------------------------------!
CALL MPI_INIT ( error )
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_Comm_size (MPI_COMM_WORLD,nprocs,error)
CALL MPI_Comm_rank (MPI_COMM_WORLD,myrank,error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
OPEN (10, FILE = 'driver.txt', STATUS = 'OLD')
READ (10,*) nyr_spin
READ (10,*) kyr_rsf
CLOSe (10)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
nland_chunk = nland / nprocs
kyr_clm = kyr_rsf + nyr_spin ! Set to year for input file-name.
ALLOCATE (B_k        (nland_chunk))
ALLOCATE (B_k_all    (nland))
ALLOCATE (B_grid     (nlon,nlat))
ALLOCATE (soilW_k    (nland_chunk))
ALLOCATE (soilW_k_all(nland))
ALLOCATE (soilW_grid (nlon,nlat))
ALLOCATE (SOM_k    (nland_chunk))
ALLOCATE (SOM_k_all(nland))
ALLOCATE (SOM_grid (nlon,nlat))
ALLOCATE (aNPP_k    (nland_chunk))
ALLOCATE (aNPP_k_all(nland))
ALLOCATE (aNPP_grid (nlon,nlat))
ALLOCATE (aRh_k    (nland_chunk))
ALLOCATE (aRh_k_all(nland))
ALLOCATE (aRh_grid (nlon,nlat))
ALLOCATE (aNBP_k    (nland_chunk))
ALLOCATE (aNBP_k_all(nland))
ALLOCATE (aNBP_grid (nlon,nlat))
ALLOCATE (larea_k    (nland_chunk))
ALLOCATE (larea_k_all(nland))
ALLOCATE (i_k_all    (nland))
ALLOCATE (j_k_all    (nland))
ALLOCATE (i_k        (nland_chunk))
ALLOCATE (j_k        (nland_chunk))
!----------------------------------------------------------------------!

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
CALL MPI_File_read(file_handle, B_k, nland_chunk, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
var_name = 'soilW'
WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
&adf10/TRENDY2021/output/HYBRID10.3_",nprocs,&
&"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
WRITE (*,*) 'Reading from ', TRIM(file_name)
! Open the file for reading.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_read(file_handle, soilW_k, nland_chunk, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
var_name = 'SOM'
WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
&adf10/TRENDY2021/output/HYBRID10.3_",nprocs,&
&"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
WRITE (*,*) 'Reading from ', TRIM(file_name)
! Open the file for reading.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_read(file_handle, SOM_k, nland_chunk, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
var_name = 'aNPP'
WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
&adf10/TRENDY2021/output/HYBRID10.3_",nprocs,&
&"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
WRITE (*,*) 'Reading from ', TRIM(file_name)
! Open the file for reading.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_read(file_handle, aNPP_k, nland_chunk, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
var_name = 'aRh'
WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
&adf10/TRENDY2021/output/HYBRID10.3_",nprocs,&
&"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
WRITE (*,*) 'Reading from ', TRIM(file_name)
! Open the file for reading.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_read(file_handle, aRh_k, nland_chunk, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
var_name = 'aNBP'
WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
&adf10/TRENDY2021/output/HYBRID10.3_",nprocs,&
&"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
WRITE (*,*) 'Reading from ', TRIM(file_name)
! Open the file for reading.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_read(file_handle, aNBP_k, nland_chunk, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

write (*,*) myrank, B_k (10), soilW_k (10), SOM_k (10), aNPP_k (10), &
 aRh_k (10), aNBP_k (10), larea_k (10), i_k (10), j_k (10)
!B_grid = fillvalue
!DO k = 1, nland_chunk
! i = i_k (k)
! j = j_k (k)
! B_grid (i,j) = B_k (k)
! !write(*,*)i,j,k,B_grid(i,j)
!END DO ! k

! How does this know which order to place stuff in B_k_all?
! Assume in processor myrank order.
CALL MPI_Gather (B_k, nland_chunk, MPI_REAL, B_k_all, nland_chunk, MPI_REAL, &
 root, MPI_COMM_WORLD, error)
CALL MPI_Gather (soilW_k, nland_chunk, MPI_REAL, soilW_k_all, nland_chunk, MPI_REAL, &
 root, MPI_COMM_WORLD, error)
CALL MPI_Gather (SOM_k, nland_chunk, MPI_REAL, SOM_k_all, nland_chunk, MPI_REAL, &
 root, MPI_COMM_WORLD, error)
CALL MPI_Gather (aNPP_k, nland_chunk, MPI_REAL, aNPP_k_all, nland_chunk, MPI_REAL, &
 root, MPI_COMM_WORLD, error)
CALL MPI_Gather (aRh_k, nland_chunk, MPI_REAL, aRh_k_all, nland_chunk, MPI_REAL, &
 root, MPI_COMM_WORLD, error)
CALL MPI_Gather (aNBP_k, nland_chunk, MPI_REAL, aNBP_k_all, nland_chunk, MPI_REAL, &
 root, MPI_COMM_WORLD, error)
CALL MPI_Gather (larea_k, nland_chunk, MPI_INTEGER, larea_k_all, nland_chunk, MPI_INTEGER, &
 root, MPI_COMM_WORLD, error)
CALL MPI_Gather (i_k, nland_chunk, MPI_INTEGER, i_k_all, nland_chunk, MPI_INTEGER, &
 root, MPI_COMM_WORLD, error)
CALL MPI_Gather (j_k, nland_chunk, MPI_INTEGER, j_k_all, nland_chunk, MPI_INTEGER, &
 root, MPI_COMM_WORLD, error)

IF (myrank == root) THEN

TA = 0.0
TB = 0.0
TW = 0.0
TSOM = 0.0
TaNPP = 0.0
TaRh = 0.0
TaNBP = 0.0
Wmax = 0.0
Bmax = 0.0
SOMmax = 0.0
aNPPmax = 0.0
aRhmax = 0.0
aNBPmax = 0.0
B_grid = fillvalue
soilW_grid = fillvalue
SOM_grid = fillvalue
aNPP_grid = fillvalue
aRh_grid = fillvalue
aNBP_grid = fillvalue
DO k = 1, nland
 i = i_k_all (k)
 j = j_k_all (k)
 B_grid (i,j) = B_k_all (k)
 soilW_grid (i,j) = soilW_k_all (k)
 SOM_grid (i,j) = SOM_k_all (k)
 aNPP_grid (i,j) = aNPP_k_all (k)
 aRh_grid (i,j) = aRh_k_all (k)
 aNBP_grid (i,j) = aNBP_k_all (k)
 TA = TA + larea_k_all (k)
 TB = TB + B_k_all (k) * larea_k_all (k)
 TW = TW + soilW_k_all (k) * larea_k_all (k)
 TSOM = TSOM + SOM_k_all (k) * larea_k_all (k)
 TaNPP = TaNPP + aNPP_k_all (k) * larea_k_all (k)
 TaRh = TaRh + aRh_k_all (k) * larea_k_all (k)
 TaNBP = TaNBP + aNBP_k_all (k) * larea_k_all (k)
 Wmax = MAX (soilW_k_all (k), Wmax)
 Bmax = MAX (B_k_all (k), Bmax)
 SOMmax = MAX (SOM_k_all (k), SOMmax)
 aNPPmax = MAX (aNPP_k_all (k), aNPPmax)
 aRhmax = MAX (aRh_k_all (k), aRhmax)
 aNBPmax = MAX (aNBP_k_all (k), aNBPmax)
 !write(*,*)i,j,k,B_grid(i,j)
END DO ! k
WRITE (*,*) 'Total land area = ',TA
WRITE (*,*) 'Total biomass = ',TB/1.0E6
WRITE (*,*) 'Total water = ',TW
WRITE (*,*) 'Total SOM = ',TSOM/1.0e6
WRITE (*,*) 'Total aNPP = ',TaNPP/1.0e6
WRITE (*,*) 'Total aRh = ',TaRh/1.0e6
WRITE (*,*) 'Total aNBP = ',TaNBP/1.0e6
WRITE (*,*) 'Wmax = ', Wmax
WRITE (*,*) 'Bmax = ', Bmax
WRITE (*,*) 'SOMmax = ', SOMmax
WRITE (*,*) 'aNPPmax = ', aNPPmax
WRITE (*,*) 'aRhmax = ', aRhmax
WRITE (*,*) 'aNBPmax = ', aNBPmax

!----------------------------------------------------------------------!
var_name = 'tmp'
kyr_clm = 1901
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
CALL CHECK (NF90_DEF_VAR (ncid, "aNPP", nf90_float, &
            dimids_two, varidaNPP))
CALL CHECK (NF90_DEF_VAR (ncid, "aRh", nf90_float, &
            dimids_two, varidaRh))
CALL CHECK (NF90_DEF_VAR (ncid, "aNBP", nf90_float, &
            dimids_two, varidaNBP))
CALL CHECK (NF90_PUT_ATT (ncid, varidW, "units", "m"))
CALL CHECK (NF90_PUT_ATT (ncid, varidB, "units", "kg[DM] m-2"))
CALL CHECK (NF90_PUT_ATT (ncid, varidSOM, "units", "kg[DM] m-2"))
CALL CHECK (NF90_PUT_ATT (ncid, varidaNPP, "units", "kg[DM] m-2 yr-1"))
CALL CHECK (NF90_PUT_ATT (ncid, varidaRh, "units", "kg[DM] m-2 yr-1"))
CALL CHECK (NF90_PUT_ATT (ncid, varidaNBP, "units", "kg[DM] m-2 yr-1"))
CALL CHECK (NF90_PUT_ATT (ncid, varidW, "_FillValue", fillvalue))
CALL CHECK (NF90_PUT_ATT (ncid, varidB, "_FillValue", fillvalue))
CALL CHECK (NF90_PUT_ATT (ncid, varidSOM, "_FillValue", fillvalue))
CALL CHECK (NF90_PUT_ATT (ncid, varidaNPP, "_FillValue", fillvalue))
CALL CHECK (NF90_PUT_ATT (ncid, varidaRh, "_FillValue", fillvalue))
CALL CHECK (NF90_PUT_ATT (ncid, varidaNBP, "_FillValue", fillvalue))
CALL CHECK (NF90_ENDDEF (ncid))
CALL CHECK (NF90_PUT_VAR (ncid, lon_varid, lon))
CALL CHECK (NF90_PUT_VAR (ncid, lat_varid, lat))
CALL CHECK (NF90_PUT_VAR (ncid,     varidW, soilW_grid))
CALL CHECK (NF90_PUT_VAR (ncid,     varidB, B_grid))
CALL CHECK (NF90_PUT_VAR (ncid,     varidSOM, SOM_grid))
CALL CHECK (NF90_PUT_VAR (ncid,     varidaNPP, aNPP_grid))
CALL CHECK (NF90_PUT_VAR (ncid,     varidaRh, aRh_grid))
CALL CHECK (NF90_PUT_VAR (ncid,     varidaNBP, aNBP_grid))
CALL CHECK (NF90_close (ncid))
!----------------------------------------------------------------------!

END IF ! myrank == root

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