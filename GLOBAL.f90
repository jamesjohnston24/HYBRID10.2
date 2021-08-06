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

  ! We are reading 2D data, a 1440 x 720 grid. 
  integer, parameter :: NX = 1440, NY = 720
  real :: data_in(NX, NY)

  ! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid

  ! Loop indexes, and error handling.
  integer :: x, y

  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
  ! the file.
  call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )

  ! Get the varid of the data variable, based on its name.
  call check( nf90_inq_varid(ncid, "carea", varid) )

  ! Read the data.
  call check( nf90_get_var(ncid, varid, data_in) )

  ! Check the data.
  do x = 1, NX
     do y = 1, NY
        if (data_in(y, x) /= (x - 1) * NY + (y - 1)) then
           print *, "data_in(", y, ", ", x, ") = ", data_in(y, x)
           stop "Stopped"
        end if
     end do
  end do

  ! Close the file, freeing all resources.
  call check( nf90_close(ncid) )

  print *,"*** SUCCESS reading example file ", FILE_NAME, "! "

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
