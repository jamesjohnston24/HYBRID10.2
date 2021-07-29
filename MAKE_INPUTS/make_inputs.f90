!======================================================================!
PROGRAM MAKE_INPUTS
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Create TRENDY input files for each processor.
!----------------------------------------------------------------------!

IMPLICIT NONE
INTEGER :: kyr_clm
CHARACTER(LEN=200) :: file_name, var_name
CHARACTER(LEN=4) :: char_year

kyr_clm = 2020
var_name = 'tmp'
WRITE (char_year, '(I4)') kyr_clm
file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/TRENDY2021/&
 &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//char_year//&
 &'.365d.noc.nc'

!TRIM(var_name)//'/crujra.v2.1.5d.'&
! &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'

!----------------------------------------------------------------------!
END PROGRAM MAKE_INPUTS
!======================================================================!
