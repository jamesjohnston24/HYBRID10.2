SUBROUTINE Diag_Global (kyr, kyr_clm)

USE mpi
USE shared
IMPLICIT NONE
INTEGER :: k, kyr, kyr_clm
REAL (KIND=DP) :: NPP_total
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
REAL (KIND=DP) :: tmp_total, TA
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
REAL (KIND=DP) :: Rh_total
REAL (KIND=DP) :: NEE_total
REAL (KIND=DP) :: B_total
REAL (KIND=DP) :: SOM_total

DO k = 1, nland_chunk
 NPP_gbox (k) = NPP_gbox (k) / FLOAT (nplots)
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
 tmp_gbox (k) = tmp_gbox (k) / FLOAT (nplots)
if (tmp_gbox(k)<0.0)then
 write (*,*)myrank,k,tmp_gbox(k)
 stop
end if
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
 Rh_gbox (k) = Rh_gbox (k) / FLOAT (nplots)
 NEE_gbox (k) = NEE_gbox (k) / FLOAT (nplots)
 ! Mean biomasss of each grid-box over plots (kg[DM] m-2).
 B_gbox (k) = SUM ( B_plot(:,k) ) / FLOAT (nplots)
 ! Mean SOM of each grid-box over plots (kg[DM] m-2).
 SOM_gbox (k) = SUM ( SOM_plot(:,k) ) / FLOAT (nplots)
END DO ! k = 1, nland_chunk
CALL MPI_Gather(NPP_gbox,nland_chunk,MPI_REAL, &
                NPP_fin,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
CALL MPI_Gather(tmp_gbox,nland_chunk,MPI_REAL, &
                tmp_fin,nland_chunk,MPI_REAL,root,MPI_COMM_WORLD,error)
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
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
 NPP_total = 0.0_DP
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
 tmp_total = 0.0_DP
 TA = 0.0_DP
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
 Rh_total = 0.0_DP
 NEE_total = 0.0_DP
 B_total = 0.0_DP
 SOM_total = 0.0_DP
 DO k = 1, nland
  NPP_total = NPP_total + source_larea (k) * NPP_fin (k)
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
if((k>8426).and.(k<8430))write(*,*)myrank,k,tmp_total,source_larea(k),tmp_fin(k)
  tmp_total = tmp_total + source_larea (k) * tmp_fin (k)
  TA = TA + source_larea (k)
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
  Rh_total = Rh_total + source_larea (k) * Rh_fin (k)
  NEE_total = NEE_total + source_larea (k) * NEE_fin(k)
!if(B_fin(k) /= 0.0_DP) write(*,*)k,source_larea(k),B_fin(k)
  B_total = B_total + source_larea (k) * B_fin (k)
  SOM_total = SOM_total + source_larea (k) * SOM_fin (k)
 END DO ! k = 1, nland
 WRITE (21,"(2(I0,' '), 5(F0.5,' '))") kyr, kyr_clm, NPP_total / 1.0D6, &
  Rh_total / 1.0D6, NEE_total / 1.0D6, B_total / 1.0E6, SOM_total / 1.0D6
 IF (ntasks == 4) &
  WRITE (*,"(2(I0,' '), 5(F0.5,' '))") kyr, kyr_clm, NPP_total / 1.0D6, &
  Rh_total / 1.0D6, NEE_total / 1.0D6, B_total / 1.0D6, SOM_total / 1.0D6
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
WRITE (*,*) 'Global mean land tmp = ', tmp_total / (TA * REAL (ntimes))
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""!
END IF

END SUBROUTINE Diag_Global
