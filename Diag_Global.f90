SUBROUTINE Diag_Global (kyr, kyr_clm)

USE mpi
USE shared
IMPLICIT NONE
INTEGER :: k, kyr, kyr_clm
REAL :: NPP_total
REAL :: Rh_total
REAL :: NEE_total
REAL :: B_total
REAL :: SOM_total

DO k = 1, nland_chunk
 NPP_gbox (k) = NPP_gbox (k) / FLOAT (nplots)
 Rh_gbox (k) = Rh_gbox (k) / FLOAT (nplots)
 NEE_gbox (k) = NEE_gbox (k) / FLOAT (nplots)
 ! Mean biomasss of each grid-box over plots (kg[DM] m-2).
 B_gbox (k) = SUM ( B_plot(:,k) ) / FLOAT (nplots)
 ! Mean SOM of each grid-box over plots (kg[DM] m-2).
 SOM_gbox (k) = SUM ( SOM_plot(:,k) ) / FLOAT (nplots)
END DO ! k = 1, nland_chunk
CALL MPI_Gather(NPP_gbox,nland_chunk,MPI_REAL, &
                NPP_fin,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
CALL MPI_Gather(Rh_gbox,nland_chunk,MPI_REAL, &
                Rh_fin,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
CALL MPI_Gather(NEE_gbox,nland_chunk,MPI_REAL, &
                NEE_fin,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
CALL MPI_Gather(B_gbox,nland_chunk,MPI_REAL, &
                B_fin,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
CALL MPI_Gather(SOM_gbox,nland_chunk,MPI_REAL, &
                SOM_fin,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
CALL MPI_Barrier ( MPI_COMM_WORLD, error )

IF (myrank == root) THEN
 NPP_total = 0.0
 Rh_total = 0.0
 NEE_total = 0.0
 B_total = 0.0
 SOM_total = 0.0
 DO k = 1, nland
  NPP_total = NPP_total + source_larea (k) * NPP_fin (k)
  Rh_total = Rh_total + source_larea (k) * Rh_fin (k)
  NEE_total = NEE_total + source_larea (k) * NEE_fin(k)
  B_total = B_total + source_larea (k) * B_fin (k)
  SOM_total = SOM_total + source_larea (k) * SOM_fin (k)
 END DO ! k = 1, nland
 WRITE (21,"(2(I0,' '), 5(F0.5,' '))") kyr, kyr_clm, NPP_total / 1.0E6, &
  Rh_total / 1.0E6, NEE_total / 1.0E6, B_total / 1.0E6, SOM_total / 1.0E6
 IF (ntasks == 4) &
  WRITE (*,"(2(I0,' '), 5(F0.5,' '))") kyr, kyr_clm, NPP_total / 1.0E6, &
  Rh_total / 1.0E6, NEE_total / 1.0E6, B_total / 1.0E6, SOM_total / 1.0E6
END IF

END SUBROUTINE Diag_Global
