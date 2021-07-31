!======================================================================!
PROGRAM MAKE_INPUTS
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Create TRENDY input files for each processor.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE netcdf
USE double
USE mpi
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
INTEGER, PARAMETER :: nlon = 720, nlat = 360, ntimes = 1460
INTEGER, PARAMETER :: nland = 67420
INTEGER, PARAMETER :: root = 0
REAL, PARAMETER :: tf = 273.15
REAL, PARAMETER :: clm_fill = 1.0E20
INTEGER :: kyr_clm, ncid, varid, i, j, k, ii, jj
INTEGER :: error, nprocs, myrank, file_handle, dest, size, errcode
REAL :: Aland ! Total land area (km^2)
REAL :: Tmean ! Global mean annual surface temperature (oC)
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: clm_in
REAL, ALLOCATABLE, DIMENSION (:,:) :: clm_k
REAL, ALLOCATABLE, DIMENSION (:,:) :: clm_buffer
REAL, ALLOCATABLE, DIMENSION (:,:) :: carea ! QD
REAL, ALLOCATABLE, DIMENSION (:,:) :: icwtr ! QD
REAL, ALLOCATABLE, DIMENSION (:,:) :: larea ! HD
REAL, ALLOCATABLE, DIMENSION (:,:) :: fwice ! HD
REAL, ALLOCATABLE, DIMENSION (:) :: larea_k
REAL, ALLOCATABLE, DIMENSION (:) :: larea_buffer
CHARACTER(LEN=200) :: file_name, var_name
CHARACTER(LEN=4) :: char_year, char_nprocs, char_myrank
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_INIT ( error )
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_Comm_size (MPI_COMM_WORLD,nprocs,error)
CALL MPI_Comm_rank (MPI_COMM_WORLD,myrank,error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
ALLOCATE (clm_in (nlon, nlat, ntimes)) ! Reverse order from netCDF file
ALLOCATE (clm_k (ntimes, nland))
ALLOCATE (larea_k (nland))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
ALLOCATE (carea (2*nlon, 2*nlat)) ! Reverse order from net CDF file
ALLOCATE (icwtr (2*nlon, 2*nlat)) ! Reverse order from net CDF file
file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/TRENDY2021/&
 &input/LUH2_GCB_2021/staticData_quarterdeg.nc'
CALL CHECK (NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid))
varid = 5 ! QD area, km^2
CALL CHECK (NF90_GET_VAR (ncid, varid, carea))
varid = 6 ! QD ice/water fraction, area fraction
CALL CHECK (NF90_GET_VAR (ncid, varid, icwtr))
CALL CHECK (NF90_CLOSE (ncid))
! Aggregate from QD to HD.
ALLOCATE (larea(nlon,nlat))
ALLOCATE (fwice(nlon,nlat))
jj = 1
DO j = 1, nlat
 ii = 1
 DO i = 1, nlon
  ! Invert as these start at NP, I presume.
  larea (i,nlat-j+1) = SUM ( carea (ii:ii+1,jj:jj+1) ) ! km^2
  fwice (i,nlat-j+1) = SUM ( icwtr (ii:ii+1,jj:jj+1) ) / 4.0 ! area fraction
  ii = ii + 2
 END DO
 jj = jj + 2
END DO ! j
DEALLOCATE (carea)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
size = ntimes * nland / nprocs
ALLOCATE (clm_buffer (ntimes,nland/nprocs))
ALLOCATE (larea_buffer (nland/nprocs))
DO kyr_clm = 1901, 1910

!kyr_clm = 2020
 var_name = 'tmp'

 IF (myrank == root) THEN

  WRITE (char_year, '(I4)') kyr_clm
  file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/TRENDY2021/&
   &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
   &char_year//'.365d.noc.nc'
  WRITE (*,*) 'Opening file: ',file_name
  CALL CHECK ( NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid ))
  varid = 4 ! Temperature (K)
  ! Origin at IDL and SP.
  CALL CHECK ( NF90_GET_VAR ( ncid, varid, clm_in ))
  CALL CHECK ( NF90_CLOSE ( ncid ))
  k = 1
  DO j = 1, nlat
   DO i = 1, nlon
    IF (clm_in (i,j,1) /= clm_fill) THEN
     clm_k (:,k) = clm_in (i,j,:)
     larea_k (k) = larea (i,j)
     k = k + 1
    END IF
   END DO ! i
  END DO ! j
  !--------------------------------------------------------------------!

  WRITE (*,*) kyr_clm, source (1, 1)
  !DEALLOCATE (clm_in)
  !DEALLOCATE (source)

  !--------------------------------------------------------------------!
  ! Compute global mean annual land surface temperature (oC).
  !--------------------------------------------------------------------!
  Aland = 0.0 ! Total land area (km^2)
  Tmean = 0.0 ! Mean land temperature (oC)
  DO j = 1, nlat
   DO i = 1, nlon
    IF (clm_in (i,j,1) /= clm_fill) THEN
     Tmean = Tmean + SUM (clm_in (i,j,:)) * larea (i,j) * &
             (1.0 - fwice (i,j))
     Aland = Aland + larea (i,j) * (1.0 - fwice (i,j))
     END IF
   END DO ! i
  END DO ! j
  Tmean = Tmean / (FLOAT (ntimes) * Aland) - tf
  WRITE (*,"('Total land area = ',F0.4,' km^2')") Aland
  WRITE (*,"('Land temperature = ',F0.4,' degC')") Tmean
  !--------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 END IF ! myrank == root
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 ! Send data to processors.
 ! Play with this if want different grid-box distributions.
 ! Here just in order of k, evenly.
 !---------------------------------------------------------------------!
 IF (MOD (nland, nprocs) /= 0.0) THEN
  WRITE (*,*) 'Problem, stopping, nland, nprocs = ', nland, nprocs
  CALL MPI_ABORT (MPI_COMM_WORLD, errcode, error)
 END IF
 IF (myrank == root) THEN
  ! Send data to each processor as 'buffer'.
  DO dest = 1, nprocs-1
    i = dest * nland / nprocs + 1
    clm_buffer (:,:) = clm_k (:,i:i+nland/nprocs-1)
    CALL MPI_SEND ( clm_buffer, size, MPI_REAL, dest, 1, MPI_COMM_WORLD, error)
    larea_buffer (:) = larea_k (i:i+nland/nprocs-1)
    CALL MPI_SEND ( larea_buffer, size/ntimes, MPI_REAL, dest, 2, MPI_COMM_WORLD, error)
  END DO
  ! Set 'clm_buffer' for root as well.
  clm_buffer (:,:) = clm_k (:,1:size)
  larea_buffer (:) = larea_k (1:size)
 ELSE
  WRITE (*,*) 'Receiving by myrank = ',myrank
  CALL MPI_RECV ( clm_buffer, size, MPI_REAL, 0, 1, MPI_COMM_WORLD, &
                  MPI_STATUS_IGNORE, error)
  CALL MPI_RECV ( larea_buffer, size/ntimes, MPI_REAL, 0, 2, MPI_COMM_WORLD, &
                  MPI_STATUS_IGNORE, error)
 END IF
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 ! Write input files for each processor.
 !---------------------------------------------------------------------!
 WRITE (char_nprocs, '(I4)') nprocs
 WRITE (char_myrank, '(I4)') myrank
 WRITE (file_name, "(A,I0.4,A,A,I0.4,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
 &adf10/TRENDY2021/input/CRUJRA2021/CRUJRA2021_",nprocs,&
 &"CPUs/",TRIM(var_name),kyr_clm,"_",myrank,".bin"
 ! Delete existing file.
 CALL MPI_File_delete(file_name, MPI_INFO_NULL, error)
 WRITE (*,*) 'Writing to ', TRIM(file_name)
 ! Open the file for writing.
 CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
  MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file_handle, error) 
 ! MPI_IO is binary output format. Write using individual file pointer.
 CALL MPI_File_write(file_handle, clm_buffer, size, &
  MPI_REAL, MPI_STATUS_IGNORE, error)
 ! Close the file.
 CALL MPI_File_Close(file_handle, error)
 !---------------------------------------------------------------------!

!----------------------------------------------------------------------!
END DO ! kyr_clm = 1901, 1910
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_FINALIZE ( error )
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CONTAINS
 SUBROUTINE check ( status )

 INTEGER, INTENT ( in ) :: status
 IF (status /= nf90_noerr) THEN
  PRINT *, TRIM (NF90_STRERROR( status ))
  STOP  "Stopped"
 END IF
 END SUBROUTINE check
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
END PROGRAM MAKE_INPUTS
!======================================================================!
