SUBROUTINE advance (kyr)

USE double
USE shared
IMPLICIT NONE
INTEGER :: t, k, kp, kyr
REAL(KIND=DP) :: ro, win, eas, ea, evap, dsoilW, BL, Tc, fT, NPP, Rh, NEE, dB, dSOM
REAL(KIND=DP) :: Tsoil, ET_SOIL, WFPS, EM, EV
REAL(KIND=DP), PARAMETER :: dt = 21600.0
REAL(KIND=DP), PARAMETER :: R = 8.3144
REAL(KIND=DP), PARAMETER :: tf = 273.15
REAL(KIND=DP), PARAMETER :: EPS = 1.0D-8 !?
REAL(KIND=DP), PARAMETER :: swc = 0.5_DP

DO t = 1, ntimes
 DO k = 1, nland_chunk
!  IF ((lon_chunk (k) == -60.25) .AND. (lat_chunk (k) == -3.25)) THEN
!   IF ((kyr_spin == 1) .AND. (t == 1)) THEN
!    OPEN (20, FILE='climate_out.txt', STATUS='UNKNOWN')
!   END IF
!   WRITE (20, "(2(1X, i0),9(1X,f0.4))") kyr_clm, t, tmp (t,kyr,k), pre (t,kyr,k), &
!    spfh (t,kyr,k), dswrf (t,kyr,k), dlwrf (t,kyr,k), pres (t,kyr,k), &
!    tmax (t,kyr,k), tmin (t,kyr,k), ws (t,kyr,k)
!   IF ((kyr_spin == nyr_spin) .AND. (t == ntimes) .AND. (k == nland_chunk)) THEN
!    CLOSE (20)
!   END IF
!  END IF
  DO kp = 1, nplots
   ro = soilW_plot (kp,k) + pre (t,kyr,k) / 1.0D3 - swc
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
   Tc = tmp (t,kyr,k) - tf ! degC
   fT = 2.0 ** (0.1 * (Tc - 25.0)) / ((1.0 + EXP (0.3 * (Tc - 36.0))) * &
        (1.0 + EXP (0.3 * (0.0 - Tc))))
   NPP = (soilW_plot (kp,k) / 0.5) * fT * 3.0_DP / (1460.0_DP * dt)
   dB = NPP - BL
   Tsoil = Tc
   IF (Tsoil > EPS) then
    ET_SOIL = 0.0326_DP + 0.00351_DP * Tsoil ** 1.652_DP - &
              (0.023953_DP * Tsoil) ** 7.19_DP
   ELSE
    ET_SOIL = 0.0326_DP
   END IF
   ET_SOIL = MAX (0.0_dp, ET_SOIL)
   ET_SOIL = MIN (1.0_dp, ET_SOIL)
   ! Convert water to %age water-filled pore space from Williams et al.
   ! Assumes micro-pore space = swc and macro-pore space = 42% of saturation
   ! content (from TEM for loam).
   WFPS = 100.0_DP * soilW_plot (kp,k) / swc
   WFPS = MIN (100.0_DP, WFPS)
   IF (WFPS < 60.0_DP) THEN
    EM = EXP ((WFPS - 60.0_DP) ** 2 / (-800.0_DP))
   ELSE
    EM = 0.000371_DP * WFPS ** 2 - 0.0748_DP * WFPS + 4.13_DP
   END IF
   EM = MAX (0.0_DP, EM)
   EM = MIN (1.0_DP, EM)
   EV = ET_SOIL * EM
   Rh = EV * SOM_plot (kp,k) * (1.0_DP / (12.5_DP * 365.0_DP * 86400.0_DP))
   dSOM = BL - Rh
   NEE = NPP - Rh ! For now!
   soilW_plot (kp,k) = soilW_plot (kp,k) + dt * dsoilW
   B_plot (kp,k) = B_plot (kp,k) + dt * dB
   SOM_plot (kp,k) = SOM_plot (kp,k) + dt * dSOM
   NPP_gbox (k) = NPP_gbox (k) + dt * NPP
   Rh_gbox (k) = Rh_gbox (k) + dt * Rh
   NEE_gbox (k) = NEE_gbox (k) + dt * NEE
  END DO ! nplots
 END DO ! k
END DO ! t

END SUBROUTINE advance
