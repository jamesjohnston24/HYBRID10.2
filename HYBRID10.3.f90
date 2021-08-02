!======================================================================!
PROGRAM HYBRID10_3
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Global model of terrestrial carbon fluxes. Runs using TRENDY climate
! and land cover.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE mpi
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
INTEGER, PARAMETER :: ntimes = 1460, nland = 67420
INTEGER :: t, k, nland_chunk
INTEGER :: error, nprocs, myrank, file_handle, size, kyr_clm
REAL :: dB, NPP, BL, fT, Tc, ro, win, eas, ea, evap
REAL, PARAMETER :: dt = 21600.0
REAL, PARAMETER :: tf = 273.15
REAL, PARAMETER :: swc = 0.5
REAL, PARAMETER :: R = 8.3144
REAL, ALLOCATABLE, DIMENSION (:,:) :: tmp ! K
REAL, ALLOCATABLE, DIMENSION (:,:) :: pre ! mm/6h
REAL, ALLOCATABLE, DIMENSION (:,:) :: spfh ! kg/kg
REAL, ALLOCATABLE, DIMENSION (:,:) :: pres ! 
REAL, ALLOCATABLE, DIMENSION (:) :: B
REAL, ALLOCATABLE, DIMENSION (:) :: soilW
CHARACTER(LEN=200) :: file_name, var_name
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_INIT ( error )
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_Comm_size (MPI_COMM_WORLD,nprocs,error)
CALL MPI_Comm_rank (MPI_COMM_WORLD,myrank,error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
nland_chunk = nland / nprocs
ALLOCATE (B(nland_chunk))
ALLOCATE (soilW(nland_chunk))
B = 0.0
soilW = 0.0
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read input data for this processor.
!----------------------------------------------------------------------!
size = ntimes * nland / nprocs
ALLOCATE (tmp(ntimes,nland/nprocs))
ALLOCATE (pre(ntimes,nland/nprocs))
ALLOCATE (spfh(ntimes,nland/nprocs))
ALLOCATE (pres(ntimes,nland/nprocs))
DO kyr_clm = 1901, 1901

 !---------------------------------------------------------------------!
 var_name = 'tmp'
 WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") &
 &"/home/adf10/rds/rds-mb425-geogscratch/&
 &adf10/TRENDY2021/input/CRUJRA2021/CRUJRA2021_",nprocs,&
 &"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
 ! Open the file for reading.
 write (*,*) 'reading from ',trim(file_name)
 CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
  MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error)
 ! MPI_IO is binary output format.
 CALL MPI_File_read(file_handle, tmp, size, &
  MPI_REAL, MPI_STATUS_IGNORE, error)
 ! Close the file.
 CALL MPI_File_Close(file_handle, error)
 !---------------------------------------------------------------------!
 var_name = 'pre'
 WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") &
 &"/home/adf10/rds/rds-mb425-geogscratch/&
 &adf10/TRENDY2021/input/CRUJRA2021/CRUJRA2021_",nprocs,&
 &"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
 ! Open the file for reading.
 CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
  MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error)
 ! MPI_IO is binary output format.
 CALL MPI_File_read(file_handle, pre, size, &
  MPI_REAL, MPI_STATUS_IGNORE, error)
 ! Close the file.
 CALL MPI_File_Close(file_handle, error)
 !---------------------------------------------------------------------!
 var_name = 'spfh'
 WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") &
 &"/home/adf10/rds/rds-mb425-geogscratch/&
 &adf10/TRENDY2021/input/CRUJRA2021/CRUJRA2021_",nprocs,&
 &"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
 ! Open the file for reading.
 CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
  MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error)
 ! MPI_IO is binary output format.
 CALL MPI_File_read(file_handle, spfh, size, &
  MPI_REAL, MPI_STATUS_IGNORE, error)
 ! Close the file.
 CALL MPI_File_Close(file_handle, error)
 !---------------------------------------------------------------------!
 var_name = 'pres'
 WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") &
 &"/home/adf10/rds/rds-mb425-geogscratch/&
 &adf10/TRENDY2021/input/CRUJRA2021/CRUJRA2021_",nprocs,&
 &"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
 ! Open the file for reading.
 CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
  MPI_MODE_RDONLY, MPI_INFO_NULL, file_handle, error)
 ! MPI_IO is binary output format.
 CALL MPI_File_read(file_handle, pres, size, &
  MPI_REAL, MPI_STATUS_IGNORE, error)
 ! Close the file.
 CALL MPI_File_Close(file_handle, error)

 !---------------------------------------------------------------------!
 DO t = 1, ntimes
  DO k = 1, nland_chunk
   ro = soilW (k) + pre (t,k) / 1.0E3 - swc
   ro = MAX (0.0, ro)
   win = (pre (t,k) / 1.0e3 - ro) / dt
   ! Pa.
   eas = 611.0 * EXP (17.27 * (tmp (t,k) - 273.15) / &
         (237.3 + tmp (t,k) - 273.15))
   ! Pa.
   ea = spfh (t,k) * pres (t,k) * 29.0E-3 / 18.0E-3
   ! Potential (aerodynamic) evaporation (m s-1).
   ! http://mgebrekiros.github.io/IntroductoryHydrology/EvaporationAndTranspiration.pdf
   evap = (eas - ea) * 0.622 * 0.4 ** 2 * &
          (29.0E-3 / (R * tmp (t,k))) * ws (t,k) / &
          (997.0_DP * (log (2.0_DP / 0.0003_DP)) ** 2)
   evap = MIN (evap, soilW_plot (kp,k) / dt)
   Tc = tmp (t,k) - tf
   fT = 2.0 ** (0.1 * (Tc - 25.0)) / ((1.0 + EXP (0.3 * (Tc - 36.0))) * &
        (1.0 + EXP (0.3 * (0.0 - Tc))))
   NPP = fT * 3.0 / (1460.0 * dt)
   BL = B (k) / (12.5 * 365.0 * 86400.0)
   dB = NPP - BL
   B (k) = B (k) + dt * dB
  END DO ! k = 1, nland_chunk
 END DO ! t = 1, ntimes
 !---------------------------------------------------------------------!

 write (*,*) kyr_clm, myrank, tmp (1,1), pre (1,1), B(100)

!----------------------------------------------------------------------!
! Write output files for each processor.
!----------------------------------------------------------------------!
var_name = 'B'
WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
&adf10/TRENDY2021/output/HYBRID10.3_",nprocs,&
&"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
! Delete existing file.
CALL MPI_File_delete(file_name, MPI_INFO_NULL, error)
WRITE (*,*) 'Writing to ', TRIM(file_name)
! Open the file for writing.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_write(file_handle, B, size/ntimes, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!

END DO ! kyr_clm

!----------------------------------------------------------------------!
CALL MPI_FINALIZE ( error )
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
END PROGRAM HYBRID10_3
!======================================================================!
