SUBROUTINE get_clm (kyr_clm, kyr)

! CRUJRA climate data are SP real, so read in as such.

use netcdf
USE mpi
USE shared
IMPLICIT NONE
INTEGER :: ncid, varid
INTEGER :: kyr_clm, kyr, l, i, j, k
CHARACTER(LEN=200) :: file_name
CHARACTER(LEN=4) :: char_year
CHARACTER(LEN=30) :: var_name

DO l = 1, 9
 IF (myrank == root) THEN
  clm_in = 0.0
  IF (l == 1) var_name = 'tmp'
  IF (l == 2) var_name = 'pre'
  IF (l == 3) var_name = 'spfh'
  IF (l == 4) var_name = 'dswrf'
  IF (l == 5) var_name = 'dlwrf'
  IF (l == 6) var_name = 'pres'
  IF (l == 7) var_name = 'tmax'
  IF (l == 8) var_name = 'tmin'
  IF (l == 9) var_name = 'wsgrd'
  WRITE (char_year, '(I4)') kyr_clm
  IF (ntasks == 4) WRITE (*,"('Reading ',I5,' ',A)") kyr_clm, TRIM(var_name)
  file_name = '/rds/user/jhj34/rds-mb425-geogscratch/adf10/FORCINGS/&
   &CRUJRA_2.1/CRUJRA2020/'//TRIM(var_name)//'/crujra.v2.1.5d.'&
   &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'
  CALL CHECK ( NF90_OPEN (trim (file_name), NF90_NOWRITE, ncid ))
  varid = 4
  ! Origin at IDL and SP.
  CALL CHECK ( NF90_GET_VAR ( ncid, varid, clm_in ))
  CALL CHECK ( NF90_CLOSE ( ncid ))
  k = 1
  DO j = 1, nlat
   DO i = 1, nlon
    IF (clm_in (i,j,1) < 1.0E10) THEN
     source (:,k) = clm_in (i,j,:)
     k = k + 1
    END IF
   END DO ! i
  END DO ! j
 END IF ! root
 CALL MPI_Scatter (source,size,MPI_REAL, &
                   result,size,MPI_REAL,root,MPI_COMM_WORLD,error)
 IF (l == 1) tmp   (:,kyr,:) = result (:,:)
 IF (l == 2) pre   (:,kyr,:) = result (:,:)
 IF (l == 3) spfh  (:,kyr,:) = result (:,:)
 IF (l == 4) dswrf (:,kyr,:) = result (:,:)
 IF (l == 5) dlwrf (:,kyr,:) = result (:,:)
 IF (l == 6) pres  (:,kyr,:) = result (:,:)
 IF (l == 7) tmax  (:,kyr,:) = result (:,:)
 IF (l == 8) tmin  (:,kyr,:) = result (:,:)
 IF (l == 9) ws    (:,kyr,:) = result (:,:)
END DO ! l

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

END SUBROUTINE get_clm
