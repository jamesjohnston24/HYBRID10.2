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
REAL, PARAMETER :: tmp_fill = 1.e+20
INTEGER :: kyr_clm, ncid, varid, i, j, k, ii, jj
INTEGER :: error, nprocs, myrank, file_handle, dest, size, errcode
REAL :: Aland ! Total land area (km^2)
REAL :: Tmean ! Global mean annual surface temperature (oC)
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: clm_in
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: tmp_in
REAL, ALLOCATABLE, DIMENSION (:,:) :: clm_k
REAL, ALLOCATABLE, DIMENSION (:,:) :: clm_buffer
REAL, ALLOCATABLE, DIMENSION (:,:) :: carea ! QD
REAL, ALLOCATABLE, DIMENSION (:,:) :: icwtr ! QD
REAL, ALLOCATABLE, DIMENSION (:,:) :: larea ! HD
REAL, ALLOCATABLE, DIMENSION (:,:) :: fwice ! HD
REAL, ALLOCATABLE, DIMENSION (:) :: larea_k
INTEGER, ALLOCATABLE, DIMENSION (:) :: i_k, j_k
REAL, ALLOCATABLE, DIMENSION (:) :: larea_buffer
INTEGER, ALLOCATABLE, DIMENSION (:) :: i_buffer, j_buffer
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
ALLOCATE (tmp_in (nlon, nlat, ntimes)) ! Reverse order from netCDF file
ALLOCATE (clm_k (ntimes, nland))
ALLOCATE (larea_k (nland))
ALLOCATE (i_k (nland))
ALLOCATE (j_k (nland))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Get tmp for land mask.
!----------------------------------------------------------------------!
IF (myrank == root) THEN
 kyr_clm = 1901
 var_name = 'tmp'
 WRITE (char_year, '(I4)') kyr_clm
 file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/TRENDY2021/&
  &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
  &char_year//'.365d.noc.nc'
 WRITE (*,*) 'Opening file: ',file_name
 CALL CHECK ( NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid ))
 varid = 4 ! Temperature (K)
 ! Origin at IDL and SP.
 CALL CHECK ( NF90_GET_VAR ( ncid, varid, tmp_in ))
 CALL CHECK ( NF90_CLOSE ( ncid ))
END IF ! myrank == root
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
ALLOCATE (i_buffer (nland/nprocs))
ALLOCATE (j_buffer (nland/nprocs))
DO kyr_clm = 1901, 1902

 var_name = 'wsgrd' ! change as wish

 IF (myrank == root) THEN

  WRITE (char_year, '(I4)') kyr_clm
  file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/TRENDY2021/&
   &input/CRUJRA2021/'//'crujra.v2.2.5d.'//TRIM(var_name)//'.'//&
   &char_year//'.365d.noc.nc'
  WRITE (*,*) 'Opening file: ',file_name
  CALL CHECK ( NF90_OPEN (TRIM (file_name), NF90_NOWRITE, ncid ))
  varid = 4 ! Temperature (K)
  !varid = Precipitation (mm/6h)
  ! Origin at IDL and SP.
  CALL CHECK ( NF90_GET_VAR ( ncid, varid, clm_in ))
  CALL CHECK ( NF90_CLOSE ( ncid ))
  k = 1
  DO j = 1, nlat
   DO i = 1, nlon
    ! Need to use tmp here to ensure common mask.
    ! pre has one more box than tmp, for example.
    IF (tmp_in (i,j,1) /= tmp_fill) THEN
     clm_k (:,k) = clm_in (i,j,:)
     larea_k (k) = larea (i,j)
     i_k (k) = i
     j_k (k) = j
     k = k + 1
    END IF
   END DO ! i
  END DO ! j
  !--------------------------------------------------------------------!

  WRITE (*,*) kyr_clm, clm_k (1, 1), larea_k (1), i_k (1), j_k (1)
  !DEALLOCATE (clm_in)
  !DEALLOCATE (source)

  !--------------------------------------------------------------------!
  ! Compute global mean annual land surface temperature (oC),
  ! or precipitation (mm/6h).
  !--------------------------------------------------------------------!
  Aland = 0.0 ! Total land area (km^2)
  Tmean = 0.0 ! Mean land temperature (oC)
  DO j = 1, nlat
   DO i = 1, nlon
    IF (tmp_in (i,j,1) /= tmp_fill) THEN
     Tmean = Tmean + SUM (clm_in (i,j,:)) * larea (i,j) * &
             (1.0 - fwice (i,j))
     Aland = Aland + larea (i,j) * (1.0 - fwice (i,j))
     END IF
   END DO ! i
  END DO ! j
  Tmean = Tmean / (FLOAT (ntimes) * Aland)
  WRITE (*,"('Total land area = ',F0.4,' km^2')") Aland
  WRITE (*,"('Land climate mean = ',F0.4,A)") Tmean, TRIM(var_name)
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
  END DO
  ! Set 'clm_buffer' for root as well.
  clm_buffer (:,:) = clm_k (:,1:size)
 ELSE
  WRITE (*,*) 'Receiving by myrank = ',myrank
  CALL MPI_RECV ( clm_buffer, size, MPI_REAL, 0, 1, MPI_COMM_WORLD, &
                  MPI_STATUS_IGNORE, error)
 END IF
 !---------------------------------------------------------------------!

 !---------------------------------------------------------------------!
 ! Write input files for each processor.
 !---------------------------------------------------------------------!
 WRITE (char_nprocs, '(I4)') nprocs
 WRITE (char_myrank, '(I4)') myrank
 ! Climate
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
END DO ! kyr_clm
!----------------------------------------------------------------------!

IF (myrank == root) THEN
 ! Send data to each processor as 'buffer'.
 DO dest = 1, nprocs-1
   i = dest * nland / nprocs + 1
   larea_buffer (:) = larea_k (i:i+nland/nprocs-1)
   i_buffer (:) = i_k (i:i+nland/nprocs-1)
   j_buffer (:) = j_k (i:i+nland/nprocs-1)
   CALL MPI_SEND ( larea_buffer, size/ntimes, MPI_REAL, dest, 2, MPI_COMM_WORLD, error)
   CALL MPI_SEND ( i_buffer, size/ntimes, MPI_INTEGER, dest, 3, MPI_COMM_WORLD, error)
   CALL MPI_SEND ( j_buffer, size/ntimes, MPI_INTEGER, dest, 4, MPI_COMM_WORLD, error)
 END DO
 ! Set 'larea_buffer' for root as well.
 larea_buffer (:) = larea_k (1:size/ntimes)
 i_buffer (:) = i_k (1:size/ntimes)
 j_buffer (:) = j_k (1:size/ntimes)
ELSE
 WRITE (*,*) 'Receiving by myrank = ',myrank
 CALL MPI_RECV ( larea_buffer, size/ntimes, MPI_REAL, 0, 2, MPI_COMM_WORLD, &
                 MPI_STATUS_IGNORE, error)
 CALL MPI_RECV ( i_buffer, size/ntimes, MPI_INTEGER, 0, 3, MPI_COMM_WORLD, &
                 MPI_STATUS_IGNORE, error)
 CALL MPI_RECV ( j_buffer, size/ntimes, MPI_INTEGER, 0, 4, MPI_COMM_WORLD, &
                 MPI_STATUS_IGNORE, error)
END IF ! myrank
!----------------------------------------------------------------------!
! larea
var_name = 'larea'
WRITE (file_name, "(A,I0.4,A,A,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
&adf10/TRENDY2021/input/LUH2_GCB_2021/static_",nprocs,&
&"CPUs/",TRIM(var_name),"_",myrank,".bin"
! Delete existing file.
CALL MPI_File_delete(file_name, MPI_INFO_NULL, error)
WRITE (*,*) 'Writing to ', TRIM(file_name)
! Open the file for writing.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_write(file_handle, larea_buffer, size/ntimes, &
 MPI_REAL, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!
! i
var_name = 'i'
WRITE (file_name, "(A,I0.4,A,A,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
&adf10/TRENDY2021/input/LUH2_GCB_2021/static_",nprocs,&
&"CPUs/",TRIM(var_name),"_",myrank,".bin"
! Delete existing file.
CALL MPI_File_delete(file_name, MPI_INFO_NULL, error)
WRITE (*,*) 'Writing to ', TRIM(file_name)
! Open the file for writing.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_write(file_handle, i_buffer, size/ntimes, &
 MPI_INTEGER, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
!----------------------------------------------------------------------!
! j
var_name = 'j'
WRITE (file_name, "(A,I0.4,A,A,A,I0.4,A)") "/home/adf10/rds/rds-mb425-geogscratch/&
&adf10/TRENDY2021/input/LUH2_GCB_2021/static_",nprocs,&
&"CPUs/",TRIM(var_name),"_",myrank,".bin"
! Delete existing file.
CALL MPI_File_delete(file_name, MPI_INFO_NULL, error)
WRITE (*,*) 'Writing to ', TRIM(file_name)
! Open the file for writing.
CALL MPI_File_open(MPI_COMM_WORLD, file_name, &
 MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, file_handle, error) 
! MPI_IO is binary output format. Write using individual file pointer.
CALL MPI_File_write(file_handle, j_buffer, size/ntimes, &
 MPI_INTEGER, MPI_STATUS_IGNORE, error)
! Close the file.
CALL MPI_File_Close(file_handle, error)
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
