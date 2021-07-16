MODULE shared

USE double
IMPLICIT NONE
LOGICAL :: RSF_In, RSF_Out
INTEGER, PARAMETER :: nyr_spin_clm = 1 ! 10
INTEGER, PARAMETER :: nyr_spin = 1 ! 120
INTEGER, PARAMETER :: syr_trans = 1901
INTEGER, PARAMETER :: eyr_trans = 1900 ! 2019
INTEGER, PARAMETER :: nlon = 720, nlat = 360, ntimes = 1460
INTEGER, PARAMETER :: ntasks = 4 ! 4 ! 32
INTEGER, PARAMETER :: nland = 80000 ! nland = 67420
INTEGER, PARAMETER :: nland_chunk = nland / ntasks
INTEGER, PARAMETER :: size = ntimes * nland / ntasks
INTEGER, PARAMETER :: nplots = 1 ! 1 ! 50
INTEGER, PARAMETER :: root = 1
INTEGER :: myrank
INTEGER :: error
INTEGER :: nyr_run = nyr_spin + (eyr_trans - syr_trans) + 1
INTEGER :: kyr_off
REAL, DIMENSION (ntimes,nland) :: source ! 0.4 GB
REAL, DIMENSION (ntimes,nland/ntasks) :: result ! 0.1 GB
REAL, DIMENSION (nland) :: source_lat
REAL, DIMENSION (nland) :: source_larea ! km2
REAL, DIMENSION (nlon,nlat,ntimes) :: clm_in ! 1.5 GB
REAL, DIMENSION (nland/ntasks) :: lon_chunk ! 
REAL, DIMENSION (nland/ntasks) :: lat_chunk !
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: tmp ! K ! 0.59 (30 yr over 20 processes)
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: pre ! mm 6-hr-1
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: spfh ! kg kg-1
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: dswrf ! J m-2 6-hr-1
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: dlwrf ! W m-2
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: pres ! Pa
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: tmax ! K
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: tmin ! K
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: ws ! m s-1
REAL, ALLOCATABLE, DIMENSION (:,:) :: B_plot
REAL, ALLOCATABLE, DIMENSION (:,:) :: SOM_plot
REAL, ALLOCATABLE, DIMENSION (:,:) :: soilW_plot
REAL, ALLOCATABLE, DIMENSION (:) :: NPP_gbox
REAL, ALLOCATABLE, DIMENSION (:) :: NPP_fin
REAL, ALLOCATABLE, DIMENSION (:) :: Rh_gbox
REAL, ALLOCATABLE, DIMENSION (:) :: Rh_fin
REAL, ALLOCATABLE, DIMENSION (:) :: NEE_gbox
REAL, ALLOCATABLE, DIMENSION (:) :: NEE_fin
REAL, ALLOCATABLE, DIMENSION (:) :: B_gbox
REAL, ALLOCATABLE, DIMENSION (:) :: SOM_gbox
REAL, ALLOCATABLE, DIMENSION (:) :: B_fin
REAL, ALLOCATABLE, DIMENSION (:) :: SOM_fin

END MODULE shared
