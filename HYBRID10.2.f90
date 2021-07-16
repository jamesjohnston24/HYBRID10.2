!======================================================================!
PROGRAM HYBRID10_2
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Global model of terrestrial carbon fluxes. Runs using TRENDY climate
! and land cover.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE double
USE netcdf
USE mpi
USE shared
!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
INTEGER :: nprocs, namelen
INTEGER :: ncid, varid, varidW, varidB, varidSOM
INTEGER :: syr, kyr, i, j, k, l, z, zt, ii, jj
INTEGER :: kyr_spin, kyr_clm, nyr_clm_alloc
INTEGER :: lon_varid, lat_varid
INTEGER :: lon_dimid, lat_dimid, vector_dimid, plot_dimid
INTEGER, DIMENSION (2) :: dimids_plots
INTEGER, DIMENSION (2) :: dimids_two
CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME) :: procname
CHARACTER(LEN=200) :: file_name
CHARACTER(LEN=4) :: char_year
CHARACTER(LEN=30) :: var_name
REAL(KIND=DP) :: before_all, after_all, before_in, after_in
REAL(KIND=DP) :: before_scatter, after_scatter
REAL, PARAMETER :: fillvalue = 1.0E20
REAL, DIMENSION (nlon) :: lon
REAL, DIMENSION (nlat) :: lat
REAL, DIMENSION (nland) :: source_lon !
REAL, ALLOCATABLE, DIMENSION (:) :: soilW_gbox
REAL, ALLOCATABLE, DIMENSION (:) :: soilW_fin
REAL, ALLOCATABLE, DIMENSION (:,:) :: soilW_grid
REAL, ALLOCATABLE, DIMENSION (:,:) :: B_grid
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
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
OPEN (10, FILE = 'driver.txt', STATUS = 'OLD')
READ (10,*) RSF ! Input from restart files?
CLOSE (10)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IF (.NOT. (RSF)) kyr_off = 0
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
source = 0.0
result = 0.0
source_lat = 0.0
source_larea = 0.0
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_Init ( error )
CALL MPI_Comm_size(MPI_COMM_WORLD,nprocs,error)
before_all = MPI_Wtime()
IF (size*nprocs /= ntimes*nland) THEN
 PRINT *, 'Invalid number of processors',nprocs,size*nprocs,nland
 CALL MPI_Abort(MPI_COMM_WORLD,1,error)
END IF
CALL MPI_Comm_rank(MPI_COMM_WORLD,myrank,error)
CALL MPI_Get_processor_name(procname,namelen,error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! For diagnostics on root:
! Read quarter-degree areas and static cover fractions of ice/water and
! aggregate both to half-degree grid-boxes. Also allocate all grid-box
! diagnostic variables on root.
!----------------------------------------------------------------------!
IF (myrank == root) THEN
 allocate (icwtr_qd (2*nlon,2*nlat)) ! QD ice/water           (fraction)
 allocate (larea_qd (2*nlon,2*nlat)) ! QD grid-box area            (km2)
 file_name = '/home/adf10/rds/rds-mb425-geogscratch/adf10/FORCINGS/&
 &LUH2_new/staticData_quarterdeg.nc'
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
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Set number of climate years for climate on each processor.
!----------------------------------------------------------------------!
nyr_clm_alloc = MAX (1, nyr_spin_clm)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Allocate climate variables on each processor.
!----------------------------------------------------------------------!
ALLOCATE (tmp(ntimes,nyr_clm_alloc,nland/ntasks))
ALLOCATE (pre(ntimes,nyr_clm_alloc,nland/ntasks))
ALLOCATE (spfh(ntimes,nyr_clm_alloc,nland/ntasks))
ALLOCATE (dswrf(ntimes,nyr_clm_alloc,nland/ntasks))
ALLOCATE (dlwrf(ntimes,nyr_clm_alloc,nland/ntasks))
ALLOCATE (pres(ntimes,nyr_clm_alloc,nland/ntasks))
ALLOCATE (tmax(ntimes,nyr_clm_alloc,nland/ntasks))
ALLOCATE (tmin(ntimes,nyr_clm_alloc,nland/ntasks))
ALLOCATE (ws(ntimes,nyr_clm_alloc,nland/ntasks))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Allocate state variables on each processor.
!----------------------------------------------------------------------!
ALLOCATE (soilW_plot(nplots,nland_chunk)) ! Soil water               (m)
ALLOCATE (B_plot    (nplots,nland_chunk)) ! Biomass         (kg[DM] m-2)
ALLOCATE (SOM_plot  (nplots,nland_chunk)) ! SOM             (kg[DM] m-2)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Allocate grid-box diagnostic variables on each processor.
!----------------------------------------------------------------------!
ALLOCATE (soilW_gbox(nland_chunk)) ! Soil water                      (m)
ALLOCATE (B_gbox    (nland_chunk)) ! Biomass                (kg[DM] m-2)
ALLOCATE (SOM_gbox  (nland_chunk)) ! SOM                    (kg[DM] m-2)
ALLOCATE (NPP_gbox  (nland_chunk)) ! NPP               (kg[DM] m-2 yr-1)
ALLOCATE (Rh_gbox   (nland_chunk)) ! Rh                (kg[DM] m-2 yr-1)
ALLOCATE (NEE_gbox  (nland_chunk)) ! NEE               (kg[DM] m-2 yr-1)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Allocate global diagnostic vectors.
!----------------------------------------------------------------------!
ALLOCATE (soilW_fin(nland)) ! Soil water                             (m)
ALLOCATE (B_fin    (nland)) ! Biomass                       (kg[DM] m-2)
ALLOCATE (SOM_fin  (nland)) ! SOM                           (kg[DM] m-2)
ALLOCATE (NPP_fin  (nland)) ! NPP                      (kg[DM] m-2 yr-1)
ALLOCATE (Rh_fin   (nland)) ! Rh                       (kg[DM] m-2 yr-1)
ALLOCATE (NEE_fin  (nland)) ! NEE                      (kg[DM] m-2 yr-1)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read cover fractions and transitions in CE 850.
! Currently not used.
!----------------------------------------------------------------------!
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
 CALL CHECK(NF90_CLOSE(ncid))
 i = NINT (0.0 + 360.0)
 j = NINT (90.0+52.0)*2
 !write(*,*)trans(:,i,j)
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Set global land array.
! Read lon and lat, and set source_lon and source_lon global vectors.
! Set source_larea global vector.
! Could perhaps get these from files above?
!----------------------------------------------------------------------!
IF (myrank == root) THEN
 clm_in = 0.0
 kyr_clm = 1901
 var_name = 'tmp'
 WRITE (char_year, '(I4)') kyr_clm
 file_name = '/rds/user/adf10/rds-mb425-geogscratch/adf10/FORCINGS/&
  &CRUJRA_2.1/CRUJRA2020/'//TRIM(var_name)//'/crujra.v2.1.5d.'&
  &//TRIM(var_name)//'.'//char_year//'.365d.noc.nc'
! WRITE (*,"('Reading from file ',A)") TRIM(file_name)
 CALL CHECK ( NF90_OPEN (trim (file_name), NF90_NOWRITE, ncid ))
 varid = 4
 ! Origin at IDL and SP.
 CALL CHECK ( NF90_GET_VAR ( ncid, varid, clm_in ))
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
 CALL CHECK ( NF90_CLOSE ( ncid ))
 k = 1
 DO j = 1, nlat
  DO i = 1, nlon
   IF (clm_in (i,j,1) < 1.0E10) THEN
    source_lon (k) = lon (i)
    source_lat (k) = lat (j)
    source_larea (k) = larea (i, j)
    k = k + 1
   END IF
  END DO ! i
 END DO ! j
END IF ! root
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Scatter lon and lat to all processors.
!----------------------------------------------------------------------!
CALL MPI_Scatter (source_lon,nland_chunk,MPI_REAL, &
               lon_chunk,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
CALL MPI_Scatter (source_lat,nland_chunk,MPI_REAL, &
               lat_chunk,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! If required, read state variables from restart file for each
! processor, else initialise them at zero.
!----------------------------------------------------------------------!
IF (RSF) THEN
 CALL RSF_In
ELSE
 soilW_plot = 0.0
 B_plot     = 0.0
 SOM_plot   = 0.0
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Initialise diagnostic file of global mean annual values.
!----------------------------------------------------------------------!
IF (myrank == root) THEN
 WRITE (file_name, "(A12,I0.5,A4)") "global_means", kyr_off, ".txt"
 OPEN (21, FILE = file_name, STATUS = 'UNKNOWN')
 WRITE (21,"('kyr kyr_clm NPP Rh NEE B SOM')")
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read all spin-up climate and scatter to all processors.
!----------------------------------------------------------------------!
syr = 1901
kyr = 1
DO kyr_clm = syr, syr + nyr_spin_clm - 1
 CALL get_clm (kyr_clm, kyr)
 kyr = kyr + 1
END DO ! kyr
!----------------------------------------------------------------------!

! Perform spin-up.
kyr = 1
kyr_clm = syr
DO kyr_spin = 1, nyr_spin
 NPP_gbox = 0.0
 Rh_gbox = 0.0
 NEE_gbox = 0.0
 IF (ntasks == 4) THEN
  WRITE (*,"('Running ',I0,' of ',I0,' spin-up years on processor ', I0)") &
   kyr_spin, nyr_spin, myrank
 END IF
 ! Advance state variables by one year.
 CALL advance (kyr)
 ! Write global fields each 10 yr.
 !IF (MOD (kyr_spin, 10) == 0) CALL Diag_Global (kyr, kyr_clm)
 CALL Diag_Global (kyr_spin+kyr_off, kyr_clm)
 kyr = kyr + 1
 kyr_clm = kyr_clm + 1
 IF (kyr > nyr_spin_clm) THEN
  kyr = 1
  kyr_clm = syr
 END IF
END DO ! kyr_spin

! Perform transient.
kyr = 1
DO kyr_clm = syr_trans, eyr_trans
 NPP_gbox = 0.0
 Rh_gbox = 0.0
 NEE_gbox = 0.0
 IF (kyr_clm > (1901+nyr_spin_clm-1)) CALL get_clm (kyr_clm, kyr)
 ! Advance state variables by one year.
 CALL advance (kyr)
 CALL Diag_Global (kyr+kyr_off+nyr_spin, kyr_clm)
 kyr = kyr + 1
END DO

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
!WRITE (*,"('All took ',F0.4,' seconds on ',I0)") &
! after_all-before_all,myrank

!IF (myrank <10) THEN
 WRITE (file_name, "(A9,I0.4,A4)") "RSF_plots", myrank, ".bin"
!ELSE
! WRITE (file_name, "(A10,I2,A4)") "RSF_plots_", myrank, ".bin"
!END IF
OPEN (22,FILE=file_name,FORM='UNFORMATTED')
WRITE (22) kyr_off + nyr_run
WRITE (22) soilW_plot
WRITE (22) B_plot
WRITE (22) SOM_plot
CLOSe (22)

CALL MPI_Barrier ( MPI_COMM_WORLD, error )
IF (myrank == root) THEN
! file_name = "RSF_plots.bin"
! call check (nf90_create (trim (file_name), cmode = nf90_clobber, &
!             ncid = ncid))
! call check (nf90_def_dim (ncid, "vector", nland, vector_dimid))
! call check (nf90_def_dim (ncid, "plot", nplots, plot_dimid))
! dimids_plots = (/ plot_dimid, vector_dimid /)
! call check (nf90_def_var (ncid, "Soil_Water", nf90_float, &
!             dimids_plot, varidW))
! call check (nf90_def_var (ncid, "Living_Biomass", nf90_float, &
!             dimids_plot, varidB))
! call check (nf90_def_var (ncid, "Soil_Organic_Matter", nf90_float, &
!             dimids_plot, varidSOM))
! call check (nf90_enddef (ncid))
! call check (nf90_put_var (ncid, varidW, soilW_plot))
! call check (nf90_put_var (ncid, varidB, B_plot))
! call check (nf90_put_var (ncid, varidSOM, SOM_plot))
! call check (nf90_close (ncid))

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
! write (*, *) 'Writing to ', trim (file_name)
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
 CLOSE (21) ! global_means.txt
END IF ! myrank == root

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

END PROGRAM HYBRID10_2
