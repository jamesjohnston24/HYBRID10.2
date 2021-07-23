MODULE shared

USE double
IMPLICIT NONE
CHARACTER(LEN=200) :: RSF_Out_file_name
CHARACTER(LEN=200) :: RSF_In_file_name
LOGICAL :: RSF_In, RSF_Out
INTEGER, PARAMETER :: nlon = 720, nlat = 360, ntimes = 1460
INTEGER, PARAMETER :: nland = 80000 ! nland = 67420
INTEGER, PARAMETER :: root = 1
INTEGER :: myrank
INTEGER :: error
INTEGER :: nyr_spin_clm
INTEGER :: nyr_spin
INTEGER :: syr_trans
INTEGER :: eyr_trans
INTEGER :: ntasks
INTEGER :: nplots
INTEGER :: nland_chunk
INTEGER :: size
INTEGER :: nyr_run
INTEGER :: kyr_off
REAL (DP), DIMENSION (ntimes,nland) :: source ! 0.4 GB
REAL (DP), DIMENSION (nland) :: source_lat
REAL (DP), DIMENSION (nland) :: source_larea ! km2
REAL (DP), DIMENSION (nlon,nlat,ntimes) :: clm_in ! 1.5 GB
REAL (DP), ALLOCATABLE, DIMENSION (:,:) :: result ! 0.1 GB
REAL, ALLOCATABLE, DIMENSION (:) :: lon_chunk
REAL, ALLOCATABLE, DIMENSION (:) :: lat_chunk
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
