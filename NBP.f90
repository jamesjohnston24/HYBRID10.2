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
real :: B, SOM, fT, NPP, BL, dB, dSOM, ET_SOIL, EV, Rh

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
 NPP = fT * 3.0
 BL = B / 12.5
 dB = NPP - BL
 ET_SOIL = 0.0326 + 0.00351 * t (iyr) ** 1.652 - (0.023953 * t (iyr)) ** 7.19
 ET_SOIL = MAX (0.0, ET_SOIL)
 ET_SOIL = MIN (1.0, ET_SOIL)
 EV = ET_SOIL
 Rh = EV * SOM / 6.25
 dSOM = BL - Rh
 B = B + dB
 SOM = SOM + dSOM
 write (*,*) kyr,B,SOM
end do

!----------------------------------------------------------------------!
END PROGRAM NBP
!======================================================================!
