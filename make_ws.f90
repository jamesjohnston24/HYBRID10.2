PROGRAM make_ws

USE double
USE netcdf
USE mpi
IMPLICIT NONE
INTEGER, PARAMETER :: ntasks = 1
INTEGER, PARAMETER :: nlon = 720, nlat = 360, ntimes = 1460
INTEGER, PARAMETER :: nland = 67420, nyr = 70
INTEGER, PARAMETER :: size = ntimes * nland / ntasks, root = 1
!INTEGER, PARAMETER :: size = nland / ntasks, root = 1
INTEGER :: nprocs, namelen, myrank, error
INTEGER :: ncid, varid
INTEGER :: syr, kyr, i, j, k, l
CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: procname
CHARACTER(LEN=200) :: file_name
CHARACTER(LEN=4) :: char_year
CHARACTER(LEN=30) :: var_name
REAL(KIND=DP) :: before, after
REAL, DIMENSION (nlon) :: lon
REAL, DIMENSION (nlat) :: lat
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: ugrd ! 1.5 GB
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: vgrd ! 1.5 GB
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: wsgrd ! 1.5 GB
INTEGER :: lon_dimid, lat_dimid, t_dimid, lon_varid, lat_varid, t_varid
INTEGER, DIMENSION (3) :: dimids_three
REAL, PARAMETER :: fillvalue = 1.0E20
REAL, DIMENSION (ntimes) :: t

ALLOCATE (ugrd(nlon,nlat,ntimes))
ALLOCATE (vgrd(nlon,nlat,ntimes))
ALLOCATE (wsgrd(nlon,nlat,ntimes))
t (1) = 1825.0
DO i = 2, ntimes
 t (i) = t (i-1) + 0.25
END DO

CALL MPI_Init ( error )
before = MPI_Wtime()

syr = 1951
DO kyr = syr, syr+nyr-1

 WRITE (char_year, '(I4)') kyr
 var_name = 'ugrd'
 file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/TRENDY2021/&
  &input/CRUJRA2021/'//'crujra.v2.2.5d.'&
  &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'
 WRITE (*,"('Reading from file ',A)") TRIM(file_name)
 CALL CHECK ( NF90_OPEN (trim (file_name), NF90_NOWRITE, ncid ))
 varid = 4
 ! Origin at IDL and SP.
 CALL CHECK ( NF90_GET_VAR ( ncid, varid, ugrd ))
 IF ((kyr == syr) .AND. (TRIM ( var_name ) == 'ugrd')) THEN
  varid = 2
  CALL CHECK ( NF90_GET_VAR ( ncid, varid, lon ))
  varid = 3
  CALL CHECK ( NF90_GET_VAR ( ncid, varid, lat ))
 END IF
 CALL CHECK ( NF90_CLOSE ( ncid ))
 var_name = 'vgrd'
 file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/TRENDY2021/&
  &input/CRUJRA2021/'//'crujra.v2.2.5d.'&
  &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'
 WRITE (*,"('Reading from file ',A)") TRIM(file_name)
 CALL CHECK ( NF90_OPEN (trim (file_name), NF90_NOWRITE, ncid ))
 varid = 4
 ! Origin at IDL and SP.
 CALL CHECK ( NF90_GET_VAR ( ncid, varid, vgrd ))
 CALL CHECK ( NF90_CLOSE ( ncid ))
 
 wsgrd = fillvalue
 DO j = 1, nlat
  DO i = 1, nlon
   IF (ugrd ( i, j, 1 ) < 1.0E6) THEN
    wsgrd ( i, j, : ) = SQRT ( ugrd ( i, j, : ) ** 2 + &
                               vgrd ( i, j, : ) ** 2 )
   END IF
  END DO
 END DO
 
 var_name = 'wsgrd'
 file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/TRENDY2021/&
  &input/CRUJRA2021/'//'crujra.v2.2.5d.'&
  &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'
 WRITE (*,"('Writing to file ',A)") TRIM(file_name)
 call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
             ncid = ncid))
 call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
 call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
 call check (nf90_def_dim (ncid, "time" , ntimes, t_dimid))
 call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
             lon_varid))
 call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
             lat_varid))
 call check (nf90_def_var (ncid, "time" , nf90_int, t_dimid, &
             t_varid))
 dimids_three = (/ lon_dimid, lat_dimid, t_dimid /)
 call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
 call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
 call check (nf90_put_att (ncid, t_varid, "units", "days since 1901-01-01 00:00:00"))
 call check (nf90_def_var (ncid, "Ground wind speed", nf90_float, &
             dimids_three, varid))
 call check (nf90_put_att (ncid, varid, "units", "m s-1"))
 call check (nf90_put_att (ncid, varid, "_FillValue", fillvalue))
 call check (nf90_enddef (ncid))
 call check (nf90_put_var (ncid, lon_varid, lon))
 call check (nf90_put_var (ncid, lat_varid, lat))
 call check (nf90_put_var (ncid, t_varid, t))
 call check (nf90_put_var (ncid,     varid, wsgrd))
 call check (nf90_close (ncid))
 
END DO ! kyr
after = MPI_Wtime()
WRITE (*,"('Took ',F0.4,' seconds')") after-before

CALL MPI_Finalize ( error )

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

END PROGRAM make_ws
