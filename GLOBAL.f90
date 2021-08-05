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

!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
END PROGRAM GLOBAL
!======================================================================!
