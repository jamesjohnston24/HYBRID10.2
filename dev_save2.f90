PROGRAM dev

USE double
USE netcdf
USE mpi
IMPLICIT NONE
INTEGER, PARAMETER :: ntasks = 4
INTEGER, PARAMETER :: nlon = 720, nlat = 360, ntimes = 1460, nland = 67420
INTEGER, PARAMETER :: size = ntimes * nland / ntasks, root = 1, nyr = 10
INTEGER, PARAMETER :: nland_chunk = nland / ntasks, nplots = 50, nyr_spin = 10
INTEGER :: nprocs, namelen, myrank, error
INTEGER :: ncid, varid, varidW, varidB, varidSOM
INTEGER :: syr, kyr, i, j, k, l, t, z, zt, kp, ii, jj, kyr_spin
INTEGER :: lon_dimid, lat_dimid, lon_varid, lat_varid
INTEGER, DIMENSION (2) :: dimids_two
CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: procname
CHARACTER(LEN=200) :: file_name
CHARACTER(LEN=4) :: char_year
CHARACTER(LEN=30) :: var_name
REAL(KIND=DP) :: before_all, after_all, before_in, after_in
REAL(KIND=DP) :: before_scatter, after_scatter
REAL(KIND=DP), PARAMETER :: dt = 21600.0
REAL(KIND=DP) :: dSOM, NPP, BL, dB, dsoilW, ro, win, evap, eas, ea
REAL, PARAMETER :: fillvalue = 1.0E20
REAL, PARAMETER :: R = 8.3144
REAL, DIMENSION (nlon) :: lon
REAL, DIMENSION (nlat) :: lat
REAL, DIMENSION (nlon,nlat,ntimes) :: clm_in ! 1.5 GB
REAL, DIMENSION (ntimes,nland) :: source ! 0.4 GB
REAL, DIMENSION (nland) :: source_lon ! 
REAL, DIMENSION (nland) :: source_lat ! 
REAL, DIMENSION (ntimes,nland/ntasks) :: result ! 0.1 GB
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
REAL, ALLOCATABLE, DIMENSION (:,:) :: soilW_plot
REAL, ALLOCATABLE, DIMENSION (:) :: soilW_gbox
REAL, ALLOCATABLE, DIMENSION (:) :: soilW_fin
REAL, ALLOCATABLE, DIMENSION (:,:) :: soilW_grid
REAL, ALLOCATABLE, DIMENSION (:,:) :: B_plot
REAL, ALLOCATABLE, DIMENSION (:) :: B_gbox
REAL, ALLOCATABLE, DIMENSION (:) :: B_fin
REAL, ALLOCATABLE, DIMENSION (:,:) :: B_grid
REAL, ALLOCATABLE, DIMENSION (:,:) :: SOM_plot
REAL, ALLOCATABLE, DIMENSION (:) :: SOM_gbox
REAL, ALLOCATABLE, DIMENSION (:) :: SOM_fin
REAL, ALLOCATABLE, DIMENSION (:,:) :: SOM_grid
REAL, ALLOCATABLE, DIMENSION (:,:) :: land
REAL, ALLOCATABLE, DIMENSION (:,:) :: icwtr_qd
REAL, ALLOCATABLE, DIMENSION (:,:) :: larea_qd
REAL, ALLOCATABLE, DIMENSION (:,:) :: icwtr
REAL, ALLOCATABLE, DIMENSION (:,:) :: larea
REAL, ALLOCATABLE, DIMENSION (:,:) :: fcover_in
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: fcover
REAL, ALLOCATABLE, DIMENSION (:,:) :: trans_in
REAL, ALLOCATABLE, DIMENSION (:,:,:) :: trans

source = 0.0
result = 0.0

CALL MPI_Init ( error )
CALL MPI_Comm_size(MPI_COMM_WORLD,nprocs,error)
before_all = MPI_Wtime()
IF (size*nprocs /= ntimes*nland) THEN
 PRINT *, 'Invalid number of processors',nprocs,size*nprocs,nland
 CALL MPI_Abort(MPI_COMM_WORLD,1,error)
END IF
CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank,error)
CALL MPI_Get_processor_name(procname,namelen,error)
WRITE (*,"('This is process ',I1,' in a communicator of ',    &
 &    I1,' processes running on processor ',A)")    &
 myrank,nprocs,procname(1:namelen)

ALLOCATE (tmp(ntimes,nyr,nland/ntasks))
ALLOCATE (pre(ntimes,nyr,nland/ntasks))
ALLOCATE (spfh(ntimes,nyr,nland/ntasks))
ALLOCATE (dswrf(ntimes,nyr,nland/ntasks))
ALLOCATE (dlwrf(ntimes,nyr,nland/ntasks))
ALLOCATE (pres(ntimes,nyr,nland/ntasks))
ALLOCATE (tmax(ntimes,nyr,nland/ntasks))
ALLOCATE (tmin(ntimes,nyr,nland/ntasks))
ALLOCATE (ws(ntimes,nyr,nland/ntasks))
ALLOCATE (soilW_plot(nplots,nland_chunk))
ALLOCATE (soilW_gbox(nland_chunk))
ALLOCATE (soilW_fin(nland))
ALLOCATE (B_plot(nplots,nland_chunk))
ALLOCATE (B_gbox(nland_chunk))
ALLOCATE (B_fin(nland))
ALLOCATE (SOM_plot(nplots,nland_chunk))
ALLOCATE (SOM_gbox(nland_chunk))
ALLOCATE (SOM_fin(nland))

IF (myrank == root) THEN
 ALLOCATE (fcover_in(nlon,nlat))
 ALLOCATE (fcover(14,nlon,nlat))
 file_name = '/home/adf10/rds/rds-mb425-geogscratch/adf10/FORCINGS/&
             &LUH2_new/states_half.nc'
 CALL CHECK(NF90_OPEN (TRIM(file_name), NF90_NOWRITE, ncid))
 DO z = 1, 14
  CALL check(nf90_inquire_variable(ncid,z+3,var_name))
  !WRITE (*,*) z,var_name
  CALL check(nf90_inq_varid(ncid, TRIM(var_name), varid))
  CALL CHECK(nf90_get_var(ncid, varid, fcover_in, &
                      !start = (/   1,    1, 1700-850+1/), &
                      start = (/   1,    1, 1/), &
                      count = (/nlon, nlat,        1/)))
  do j = 1, nlat
   fcover (z,:,j) = fcover_in (:,nlat-j+1)
  END DO
 END DO
 CALL check(nf90_close(ncid))
 ALLOCATE (trans_in(nlon,nlat))
 ALLOCATE (trans(118,nlon,nlat))
 file_name = '/home/adf10/rds/rds-mb425-geogscratch/adf10/FORCINGS/&
             &LUH2_new/transitions_half.nc'
 CALL CHECK(NF90_OPEN (TRIM(file_name), NF90_NOWRITE, ncid))
 DO zt = 1, 118
  CALL check(nf90_inquire_variable(ncid,zt+3,var_name))
  !WRITE (*,*) zt,var_name
  CALL check(nf90_inq_varid(ncid, TRIM(var_name), varid))
  CALL CHECK(nf90_get_var(ncid, varid, trans_in, &
                      !start = (/   1,    1, 1700-850+1/), &
                      start = (/   1,    1, 1/), &
                      count = (/nlon, nlat,        1/)))
  do j = 1, nlat
   trans (zt,:,j) = trans_in (:,nlat-j+1)
  END DO
 END DO
 CALL check(nf90_close(ncid))
 i = NINT (0.0 + 360.0)
 j = NINT (90.0+52.0)*2
 !write(*,*)trans(:,i,j)
END IF

syr = 1901
DO kyr = syr, syr+nyr-1
 DO l = 1, 9

 IF (myrank == root) THEN
  clm_in = 0.0
  IF (l == 1) var_name = 'tmp'
  IF (l == 2) var_name = 'pre'
  IF (l == 3) var_name = 'spfh'
  IF (l == 4) var_name = 'dswrf'
  IF (l == 5) var_name = 'dlwrf'
  IF (l == 6) var_name = 'pres'
  IF (l == 7) var_name = 'tmax'
  IF (l == 8) var_name = 'tmin'
  IF (l == 9) var_name = 'wsgrd'
  before_in = MPI_Wtime()
  WRITE (char_year, '(I4)') kyr
  file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/FORCINGS/&
   &CRUJRA_2.1/CRUJRA2020/'//TRIM(var_name)//'/crujra.v2.1.5d.'&
   &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'
  WRITE (*,"('Reading from file ',A)") TRIM(file_name)
  CALL CHECK ( NF90_OPEN (trim (file_name), NF90_NOWRITE, ncid ))
  varid = 4
  ! Origin at IDL and SP.
  CALL CHECK ( NF90_GET_VAR ( ncid, varid, clm_in ))
  IF ((kyr == syr) .AND. (TRIM ( var_name ) == 'tmp')) THEN
   varid = 2
   CALL CHECK ( NF90_GET_VAR ( ncid, varid, lon ))
   varid = 3
   CALL CHECK ( NF90_GET_VAR ( ncid, varid, lat ))
   ALLOCATE (land(nlon,nlat))
   DO j = 1, nlat
    DO i = 1, nlon
     IF (clm_in (i,j,1) /= fillvalue) THEN
      land (i,j) = 1.0
     ELSE
      land (i,j) = 0.0
     ENDIF
    END DO
   END DO
  END IF
  CALL CHECK ( NF90_CLOSE ( ncid ))
  k = 1
  DO j = 1, nlat
   DO i = 1, nlon
    IF (clm_in (i,j,1) < 1.0E10) THEN
     source (:,k) = clm_in (i,j,:)
     IF ((l == 1) .AND. (kyr == syr)) THEN
      source_lon (k) = lon (i)
      source_lat (k) = lat (j)
     END IF
     k = k + 1
    END IF
   END DO
  END DO
  after_in = MPI_Wtime()
  WRITE (*,"('In took ',F0.4,' seconds')") after_in-before_in
 END IF ! root

 before_scatter = MPI_Wtime()
 CALL MPI_Scatter (source,size,MPI_REAL, &
                   result,size,MPI_REAL,root,MPI_COMM_WORLD,error)

 IF (l == 1) tmp   (:,kyr-syr+1,:) = result (:,:)
 IF (l == 2) pre   (:,kyr-syr+1,:) = result (:,:)
 IF (l == 3) spfh  (:,kyr-syr+1,:) = result (:,:)
 IF (l == 4) dswrf (:,kyr-syr+1,:) = result (:,:)
 IF (l == 5) dlwrf (:,kyr-syr+1,:) = result (:,:)
 IF (l == 6) pres  (:,kyr-syr+1,:) = result (:,:)
 IF (l == 7) tmax  (:,kyr-syr+1,:) = result (:,:)
 IF (l == 8) tmin  (:,kyr-syr+1,:) = result (:,:)
 IF (l == 9) ws    (:,kyr-syr+1,:) = result (:,:)
 after_scatter = MPI_Wtime()
 WRITE (*,"('Scatter took ',F0.4,' seconds on ',I0)") &
  after_scatter-before_scatter,myrank
 END DO ! l
 
 CALL MPI_Scatter (source_lon,nland_chunk,MPI_REAL, &
                   lon_chunk,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
 CALL MPI_Scatter (source_lat,nland_chunk,MPI_REAL, &
                   lat_chunk,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
 
END DO ! kyr

soilW_plot = 0.0
B_plot = 0.0
SOM_plot = 0.0
kyr = 1
DO kyr_spin = 1, nyr_spin
 WRITE (*,"('Running ',I0,' of ',I0,' spin-up years on processor ', I0)") &
  kyr_spin, nyr_spin, myrank
 DO t = 1, ntimes
  DO k = 1, nland_chunk
   IF ((lon_chunk (k) == -60.25) .AND. (lat_chunk (k) == -3.25)) THEN
    IF ((kyr_spin == 1) .AND. (t == 1)) THEN
     OPEN (20, FILE='climate_out.txt', STATUS='UNKNOWN')
    END IF
    WRITE (20, "(2(1X, i0),9(1X,f0.4))") kyr, t, tmp (t,kyr,k), pre (t,kyr,k), &
     spfh (t,kyr,k), dswrf (t,kyr,k), dlwrf (t,kyr,k), pres (t,kyr,k), &
     tmax (t,kyr,k), tmin (t,kyr,k), ws (t,kyr,k)
    IF ((kyr_spin == nyr_spin) .AND. (t == ntimes) .AND. (k == nland_chunk)) THEN
     CLOSE (20)
    END IF
   END IF
   DO kp = 1, nplots
    ro = soilW_plot (kp,k) + pre (t,kyr,k) / 1.0D3 - 0.5_DP
    ro = MAX (0.0_DP, ro)
    win = (pre (t,kyr,k) / 1.0D3 - ro) / dt
    ! Pa.
    eas = 611.0_DP * EXP (17.27_DP * (tmp (t,kyr,k) - 273.15_DP) / &
           (237.3_DP + tmp (t,kyr,k) - 273.15_DP))
    ! Pa.
    ea = spfh (t,kyr,k) * pres (t,kyr,k) * 29.0D-3 / 18.0D-3
    ! Potential (aerodynamic) evaporation (m s-1).
    ! http://mgebrekiros.github.io/IntroductoryHydrology/EvaporationAndTranspiration.pdf
    evap = (eas - ea) * 0.622_DP * 0.4_DP ** 2 * &
           (29.0D-3 / (R * tmp (t,kyr,k))) * ws (t,kyr,k) / &
           (997.0 * (log (2.0 / 0.0003)) ** 2)
    evap = MIN (evap, soilW_plot (kp,k) / dt)
    dsoilW = win - evap
    BL = B_plot (kp,k) / (12.5_DP * 365.0_DP * 86400.0_DP)
    NPP = (soilW_plot (kp,k) / 0.5) * (tmp (t,kyr,k) / 298.0_DP) * &
          3.0_DP / (1460.0_DP * dt)
    dB = NPP - BL
    dSOM = BL - (2.0_DP ** ((tmp (t,kyr,k) - 273.15_DP) / 10.0_DP)) * SOM_plot (kp,k) * &
           (10.0 / (12.5_DP * 365.0_DP * 86400.0_DP))
    soilW_plot (kp,k) = soilW_plot (kp,k) + dt * dsoilW
    B_plot (kp,k) = B_plot (kp,k) + dt * dB
    SOM_plot (kp,k) = SOM_plot (kp,k) + dt * dSOM
   END DO
  END DO
 END DO
 kyr = kyr + 1
 IF (kyr > nyr) kyr = 1
END DO ! kyr_spin
DO k = 1, nland_chunk
 soilW_gbox (k) = SUM ( soilW_plot(:,k) ) / FLOAT (nplots)
 B_gbox (k) = SUM ( B_plot(:,k) ) / FLOAT (nplots)
 SOM_gbox (k) = SUM ( SOM_plot(:,k) ) / FLOAT (nplots)
END DO
CALL MPI_Gather(soilW_gbox,nland_chunk,MPI_REAL, &
                soilW_fin,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
CALL MPI_Gather(B_gbox,nland_chunk,MPI_REAL, &
                B_fin,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
CALL MPI_Gather(SOM_gbox,nland_chunk,MPI_REAL, &
                SOM_fin,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
after_all = MPI_Wtime()
WRITE (*,"('All took ',F0.4,' seconds on ',I0)") &
 after_all-before_all,myrank
 
IF (myrank == root) THEN
 allocate (icwtr_qd (2*nlon,2*nlat)) ! QD ice/water            (fraction)
 allocate (larea_qd (2*nlon,2*nlat)) ! QD grid-box area             (km2)
 file_name = '/home/adf10/rds/rds-mb425-geogscratch/adf10/FORCINGS/&
 &LUH2_new/staticData_quarterdeg.nc'
 write (*, *) 'Reading from ', trim (file_name)
 call check (nf90_open (trim (file_name), nf90_nowrite, ncid))
 varid = 7
 call check (nf90_get_var (ncid, varid, icwtr_qd))
 varid = 9
 call check (nf90_get_var (ncid, varid, larea_qd))
 call check (nf90_close (ncid))
 ALLOCATE (icwtr(nlon,nlat))
 ALLOCATE (larea(nlon,nlat))
 jj = 1
 do j = 1, nlat
  ii = 1
  do i = 1, nlon
   icwtr (i, nlat-j+1) = sum (icwtr_qd (ii:ii+1, jj:jj+1)) / 4.0
   larea (i, nlat-j+1) = sum (larea_qd (ii:ii+1, jj:jj+1))
   ii = ii + 2
  end do
  jj = jj + 2
 end do
 ALLOCATE (soilW_grid(nlon,nlat))
 ALLOCATE (B_grid(nlon,nlat))
 ALLOCATE (SOM_grid(nlon,nlat))
 k = 1
 DO j = 1, nlat
  DO i = 1, nlon
   IF (land (i,j) > 0.0) THEN
    soilW_grid (i,j) = (1.0 - icwtr (i,j)) * soilW_fin (k)
    B_grid (i,j) = (1.0 - icwtr (i,j)) * B_fin (k)
    SOM_grid (i,j) = (1.0 - icwtr (i,j)) * SOM_fin (k)
    k = k + 1
   ELSE
    soilW_grid (i,j) = fillvalue
    B_grid (i,j) = fillvalue
    SOM_grid (i,j) = fillvalue
   END IF
  END DO
 END DO
 file_name = "fields_grid.nc"
 write (*, *) 'Writing to ', trim (file_name)
 call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
             ncid = ncid))
 call check (nf90_def_dim (ncid, "longitude", nlon, lon_dimid))
 call check (nf90_def_dim (ncid, "latitude" , nlat, lat_dimid))
 call check (nf90_def_var (ncid, "longitude", nf90_float, lon_dimid, &
             lon_varid))
 call check (nf90_def_var (ncid, "latitude" , nf90_float, lat_dimid, &
             lat_varid))
 dimids_two = (/ lon_dimid, lat_dimid /)
 call check (nf90_put_att (ncid, lon_varid, "units", "degrees_east"))
 call check (nf90_put_att (ncid, lat_varid, "units", "degrees_north"))
 call check (nf90_def_var (ncid, "Soil_Water", nf90_float, &
             dimids_two, varidW))
 call check (nf90_def_var (ncid, "Living_Biomass", nf90_float, &
             dimids_two, varidB))
 call check (nf90_def_var (ncid, "Soil_Organic_Matter", nf90_float, &
             dimids_two, varidSOM))
 call check (nf90_put_att (ncid, varidW, "units", "m"))
 call check (nf90_put_att (ncid, varidB, "units", "kg[DM] m-2"))
 call check (nf90_put_att (ncid, varidSOM, "units", "kg[DM] m-2"))
 call check (nf90_put_att (ncid, varidW, "_FillValue", fillvalue))
 call check (nf90_put_att (ncid, varidB, "_FillValue", fillvalue))
 call check (nf90_put_att (ncid, varidSOM, "_FillValue", fillvalue))
 call check (nf90_enddef (ncid))
 call check (nf90_put_var (ncid, lon_varid, lon))
 call check (nf90_put_var (ncid, lat_varid, lat))
 call check (nf90_put_var (ncid,     varidW, soilW_grid))
 call check (nf90_put_var (ncid,     varidB, B_grid))
 call check (nf90_put_var (ncid,     varidSOM, SOM_grid))
 call check (nf90_close (ncid))
END IF

CALL MPI_Finalize ( error )

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

END PROGRAM dev
