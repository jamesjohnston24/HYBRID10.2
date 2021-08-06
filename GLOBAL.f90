!======================================================================!
PROGRAM GLOBAL
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read all temperature files from TRENDY20021 and create simple text
! file of global annual values. Keep as simple as possible.
! Then check values are good.
! Then do same for all climate variables.
! The use these global land means to drive a global NBP model.
! See if more like HYBRID10.2 or 10.3.
! See which might be wrong and why? Fix!
! Also, use GLOBAL model to look at source/sink effects...
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE netcdf
IMPLICIT NONE
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! First read land areas in each grid-box.
! staticData_quarterdeg.nc contains quarter-degree areas of grid-boxes.
! Read and see if as expected from cosine rule.
! Use unidata.ucar.edu method (simple_xy_rd.f90)
!----------------------------------------------------------------------!

! This is the name of the data file we will read. 
  character (len = *), parameter :: FILE_NAME = "/home/adf10/rds/&
  &rds-mb425-geogscratch/adf10/TRENDY2021/input/LUH2_GCB_2021/&
  &staticData_quarterdeg.nc"
  character (len = *), parameter :: FILE_NAME_tmp = "/home/adf10/rds/&
  &rds-mb425-geogscratch/adf10/TRENDY2021/input/CRUJRA2021/&
  &crujra.v2.2.5d.tmp.1901.365d.noc.nc"
character (len = *), parameter :: FILE_NAME_ptbio = "ptbio.nc"

INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12)
INTEGER, PARAMETER :: sp = KIND (1E0)

  ! We are reading 2D data, a 1440 x 720 grid.
INTEGER, PARAMETER :: NDIMS = 2
  integer, parameter :: NX = 1440, NY = 720
  integer, parameter :: NX_tmp = 720, NY_tmp = 360, NTIMES = 1460
  REAL (KIND=dp) :: data_in_lon (NX)
  REAL (KIND=dp) :: data_in_lat (NY)
  REAL (KIND=sp) :: data_in_lon_tmp (NX_tmp)
  REAL (KIND=sp) :: data_in_lat_tmp (NY_tmp)
  REAL (KIND=sp) :: data_in_tmp(NX_tmp, NY_tmp, NTIMES)
  REAL (KIND=sp) :: data_in_ptbio(NX, NY)
  REAL (KIND=dp) :: data_in_carea(NX, NY)
REAL (KIND=DP) :: sum_carea, tmp_mean, carea_land
REAL (KIND=DP) :: carea_tmp (NX_tmp, NY_tmp)
REAL (KIND=SP), PARAMETER :: tmp_fill = 1.0E+20

  ! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid_lon, varid_lat, varid_ptbio, varid_carea
INTEGER :: varid_tmp
INTEGER :: x_dimid, y_dimid, dimids (NDIMS), dimid_lon, dimid_lat

  ! Loop indexes, and error handling.
  integer :: x, y, i, j
INTEGER :: ntmp, it

  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
  ! the file.
  call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )

  ! Get the varid of the data variable, based on its name.
  ! Data starts at (-179.875; 89.875).
  call check( nf90_inq_varid(ncid, "lon", varid_lon) )
  call check( nf90_inq_varid(ncid, "lat", varid_lat) )
  call check( nf90_inq_varid(ncid, "ptbio", varid_ptbio) )
  call check( nf90_inq_varid(ncid, "carea", varid_carea) )

  ! Read the data.
  call check( nf90_get_var(ncid, varid_lon, data_in_lon) )
  call check( nf90_get_var(ncid, varid_lat, data_in_lat) )
  call check( nf90_get_var(ncid, varid_ptbio, data_in_ptbio) )
  call check( nf90_get_var(ncid, varid_carea, data_in_carea) )

  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )

  print *,"*** SUCCESS reading file ", FILE_NAME, "! "

!----------------------------------------------------------------------!
! Total area seems to take into account bulge.
sum_carea = SUM (data_in_carea)
PRINT *, "Sum of carea = ", sum_carea
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Write potential vegetation carbon for looking at map.
!----------------------------------------------------------------------!

  ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  ! overwrite this file, if it already exists.
  call check( nf90_create(FILE_NAME_ptbio, NF90_CLOBBER, ncid) )

  ! Define the dimensions. NetCDF will hand back an ID for each. 
  call check( nf90_def_dim(ncid, "lon", NX, x_dimid) )
  call check( nf90_def_dim(ncid, "lat", NY, y_dimid) )

  ! The dimids array is used to pass the IDs of the dimensions of
  ! the variables. Note that in fortran arrays are stored in
  ! column-major format.
  dimid_lon = x_dimid
  dimid_lat = y_dimid
  dimids =  (/ x_dimid, y_dimid /)

  ! Define the variables.
  call check( nf90_def_var(ncid, "lon", NF90_DOUBLE, dimid_lon, varid_lon) )
  call check( nf90_def_var(ncid, "lat", NF90_DOUBLE, dimid_lat, varid_lat) )
  call check( nf90_def_var(ncid, "ptbio", NF90_FLOAT, dimids, varid_ptbio) )

  ! End define mode. This tells netCDF we are done defining metadata.
  call check( nf90_enddef(ncid) )

  ! Write the data to the file.
  call check( nf90_put_var(ncid, varid_lon, data_in_lon) )
  call check( nf90_put_var(ncid, varid_lat, data_in_lat) )
  call check( nf90_put_var(ncid, varid_ptbio, data_in_ptbio) )

  ! Close the file. This frees up any internal netCDF resources
  ! associated with the file, and flushes any buffers.
  call check( nf90_close(ncid) )

  print *,"*** SUCCESS writing file ", FILE_NAME_ptbio, "! "
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
  ! the file.
  call check( nf90_open(FILE_NAME_tmp, NF90_NOWRITE, ncid) )

  ! Get the varid of the data variable, based on its name.
  ! Data starts at (-179.75; -89.75).
  call check( nf90_inq_varid(ncid, "lon", varid_lon) )
  call check( nf90_inq_varid(ncid, "lat", varid_lat) )
  call check( nf90_inq_varid(ncid, "tmp", varid_tmp) )

  ! Read the data.
  call check( nf90_get_var(ncid, varid_lon, data_in_lon_tmp) )
  call check( nf90_get_var(ncid, varid_lat, data_in_lat_tmp) )
  call check( nf90_get_var(ncid, varid_tmp, data_in_tmp) )

  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )

  print *,"*** SUCCESS reading file ", FILE_NAME_tmp, "! "
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Create file of HD areas in same format as the climate data.
!----------------------------------------------------------------------!
carea_tmp = 0.0
! Aggregate QD areas to HD, check sum to same.
j = 720
DO y = 1, NY_tmp
 i = 1
 DO x = 1, NX_tmp
  carea_tmp (x,y) = carea_tmp (x,y) + data_in_carea (i,j)
  carea_tmp (x,y) = carea_tmp (x,y) + data_in_carea (i+1,j)
  carea_tmp (x,y) = carea_tmp (x,y) + data_in_carea (i,j-1)
  carea_tmp (x,y) = carea_tmp (x,y) + data_in_carea (i+1,j-1)
  i = i + 2
 END DO
 j = j - 2
END DO
PRINT *, "Sum of carea_tmp = ", SUM (carea_tmp)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Compute annual global land temperatures.
!----------------------------------------------------------------------!
carea_land = 0.0_DP
tmp_mean = 0.0_DP
ntmp = 0
it = 1
DO y = 1, NY_tmp
 DO x = 1, NX_tmp
  IF (data_in_tmp (x,y,1) /= tmp_fill) THEN
   carea_land = carea_land + carea_tmp (x,y)
   tmp_mean = tmp_mean + SUM (data_in_tmp (x,y,1:NTIMES))
   ntmp = ntmp + 1
  END IF
 END DO
END DO
tmp_mean = tmp_mean / DBLE (ntmp*NTIMES)
PRINT *, "carea_land = ", carea_land, carea_land/SUM (carea_tmp)
PRINT *, "mean tmp = ", tmp_mean, ntmp
!----------------------------------------------------------------------!

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
END PROGRAM GLOBAL
!======================================================================!
