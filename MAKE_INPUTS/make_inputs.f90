!======================================================================!
PROGRAM MAKE_INPUTS
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Create TRENDY input files for each processor.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE netcdf
USE double
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
INTEGER, PARAMETER :: nlon = 720, nlat = 360, ntimes = 1460
INTEGER, PARAMETER :: nland = 67420
REAL, PARAMETER :: tf = 273.15
REAL, PARAMETER :: clm_fill = 1.0E20
INTEGER :: kyr_clm, ncid, varid, i, j, k, ii, jj
REAL :: Aland ! Total land area (km^2)
REAL :: Tmean ! Global mean annual surface temperature (oC)
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: clm_in
REAL, ALLOCATABLE, DIMENSION (:,:) :: source
REAL, ALLOCATABLE, DIMENSION (:,:) :: carea ! QD
REAL, ALLOCATABLE, DIMENSION (:,:) :: icwtr ! QD
REAL, ALLOCATABLE, DIMENSION (:,:) :: larea ! HD
REAL, ALLOCATABLE, DIMENSION (:,:) :: fwice ! HD
CHARACTER(LEN=200) :: file_name, var_name
CHARACTER(LEN=4) :: char_year
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
ALLOCATE (clm_in (nlon, nlat, ntimes)) ! Reverse order from netCDF file
ALLOCATE (source (ntimes, nland))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
ALLOCATE (carea (2*nlon, 2*nlat)) ! Reverse order from net CDF file
ALLOCATE (icwtr (2*nlon, 2*nlat)) ! Reverse order from net CDF file
file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/TRENDY2021/&
 &input/LUH2_GCB_2021/staticData_quarterdeg.nc'
CALL CHECK (NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid))
varid = 5 ! QD area, km^2
CALL CHECK (NF90_GET_VAR (ncid, varid, carea))
varid = 6 ! QD ice/water fraction, area fraction
CALL CHECK (NF90_GET_VAR (ncid, varid, icwtr))
CALL CHECK (NF90_CLOSE (ncid))
! Aggregate from QD to HD.
ALLOCATE (larea(nlon,nlat))
jj = 1
DO j = 1, nlat
 ii = 1
 DO i = 1, nlon
  larea (i,j) = SUM ( carea (ii:ii+1,jj:jj+1) ) ! km^2
  fwice (i,j) = SUM ( icwtr (ii:ii+1,jj:jj+1) ) / 4.0 ! area fraction
  ii = ii + 2
 END DO
 jj = jj + 2
END DO ! j
DEALLOCATE (carea)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
kyr_clm = 2020
var_name = 'tmp'
WRITE (char_year, '(I4)') kyr_clm
file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/TRENDY2021/&
 &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
 &char_year//'.365d.noc.nc'
WRITE (*,*) 'Opening file: ',file_name
CALL CHECK ( NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid ))
varid = 4 ! Temperature (K)
! Origin at IDL and SP.
CALL CHECK ( NF90_GET_VAR ( ncid, varid, clm_in ))
CALL CHECK ( NF90_CLOSE ( ncid ))
k = 1
DO j = 1, nlat
 DO i = 1, nlon
  IF (clm_in (i,j,1) /= clm_fill) THEN
   source (:,k) = clm_in (i,j,:)
   k = k + 1
  END IF
 END DO ! i
END DO ! j
!----------------------------------------------------------------------!

WRITE (*,*) source (1, 1)
!DEALLOCATE (clm_in)
DEALLOCATE (source)

!----------------------------------------------------------------------!
! Compute global mean annual land surface temperature (oC).
!----------------------------------------------------------------------!
Aland = 0.0 ! Total land area (km^2)
Tmean = 0.0 ! Mean land temperature (oC)
DO j = 1, nlat
 DO i = 1, nlon
  IF (clm_in (i,j,1) /= clm_fill) THEN
   Tmean = Tmean + SUM (clm_in (i,j,:)) * larea (i,j)
   Aland = Aland + larea (i,j)
  END IF
 END DO ! i
END DO ! j
Tmean = Tmean / (FLOAT (ntimes) * Aland) - tf
WRITE (*,"('Total land area = ',F0.4,' km^2')") Aland
WRITE (*,"('Land temperature = ',F0.4,' degC')") Tmean
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

!----------------------------------------------------------------------!
END PROGRAM MAKE_INPUTS
!======================================================================!
