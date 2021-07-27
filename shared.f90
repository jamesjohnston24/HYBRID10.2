MODULE shared

USE double
IMPLICIT NONE
CHARACTER(LEN=200) :: RSF_Out_file_name
CHARACTER(LEN=200) :: RSF_In_file_name
LOGICAL :: RSF_In, RSF_Out
INTEGER, PARAMETER :: nlon = 720, nlat = 360, ntimes = 1460
INTEGER, PARAMETER :: nland = 67424 ! nland = 67420
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
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:) :: source ! 0.4 GB
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: source_lat
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: source_larea ! km2
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: clm_in ! 1.5 GB
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:) :: result ! 0.1 GB
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: lon_chunk
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: lat_chunk
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: tmp ! K ! 0.59 (30 yr over 20 processes)
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: pre ! mm 6-hr-1
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: spfh ! kg kg-1
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: dswrf ! J m-2 6-hr-1
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: dlwrf ! W m-2
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: pres ! Pa
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: tmax ! K
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: tmin ! K
REAL (KIND=SP), ALLOCATABLE, DIMENSION (:,:,:) :: ws ! m s-1
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:,:) :: B_plot
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:,:) :: SOM_plot
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:,:) :: soilW_plot
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: NPP_gbox
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: NPP_fin
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: Rh_gbox
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: Rh_fin
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: NEE_gbox
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: NEE_fin
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: B_gbox
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: SOM_gbox
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: B_fin
REAL (KIND=DP), ALLOCATABLE, DIMENSION (:) :: SOM_fin

END MODULE shared
