!======================================================================!
PROGRAM NBP
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Global model of terrestrial carbon fluxes. Runs using TRENDY climate
! and land cover.
!----------------------------------------------------------------------!


implicit none
integer, parameter :: nyr = 2020 - 1901 + 1
integer :: kyr, iyr
real, dimension (nyr) :: t
real :: B, SOM, fT

open (10,file="tmp_mean.txt",status="old")
do kyr = 1, nyr
 read (10,*) iyr, t (kyr)
end do
close (10)

B = 0.0
SOM = 0.0
iyr = 1
do kyr = 1, 100
 iyr = iyr + 1
 if (iyr > 20) iyr = 1
 fT = 2.0 ** (0.1 * (t (iyr) - 25.0)) / ((1.0 + EXP (0.3 * (t (iyr) - 36.0))) * &
      (1.0 + EXP (0.3 * (0.0 - t (iyr)))))
 write (*,*) iyr,t(iyr),fT
end do

!----------------------------------------------------------------------!
END PROGRAM NBP
!======================================================================!
