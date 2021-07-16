SUBROUTINE Get_RSF

USE mpi
USE shared
IMPLICIT NONE
CHARACTER(LEN=200) :: file_name

!WRITE (file_name, "(A15,I0.4,A4)") "RSF_plots", myrank, ".bin"
!WRITE (file_name, "(A20,I0.4,A4)") "RUN2/SAVE5/RSF_plots", myrank, ".bin"
WRITE (file_name, "(A0,I0.4,A4)") TRIM(RSF_In_file_name), myrank, ".bin"
write (*,*)file_name
stop
OPEN (10,FILE=file_name,FORM='UNFORMATTED')
READ (10) kyr_off
READ (10) soilW_plot
READ (10) B_plot
READ (10) SOM_plot
CLOSe (10)

END SUBROUTINE Get_RSF
