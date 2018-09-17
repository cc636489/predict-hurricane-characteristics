        PROGRAM MODEL_HEADING 
        
        IMPLICIT NONE
        INTEGER(8) :: NR
        DOUBLE PRECISION :: YEAR
        DOUBLE PRECISION :: DISTANCE
        DOUBLE PRECISION :: HEADING
        DOUBLE PRECISION :: PRESSURE
        DOUBLE PRECISION :: RADIUS
        DOUBLE PRECISION :: VFSPEED

        OPEN(UNIT=1,FILE="DATASET.CSV")
        CALL READRC(NR)

        !INPUT DATAD DATA FOR EACH PARAMETERS
        DO I = 1,NR
           READ(1,1000) YEAR(I),DISRANCE(I),HEADING(I),PRESSURE(I),&
                        RADIUS(I),VFSPEED(I)
        ENDDO
        REWIND(1)
        CLOSE(1)

        !DOING OPTIMIZATION FOR HD
        


        END PROGRAM MODEL_HEADING
        

        SUBROUTINE READRC(NR)

        INTEGER,INTENT(OUT)  :: NR
        INTEGER              :: I,STAT

        !read the row of the sample
        NR = 0
        do while(.true.)
               read(1,*,IOSTAT=STAT)
               if(STAT==0) then
                 NR = NR + 1
               else if(STAT>0) then
                 write(*,*) "SOMETHING WRONG DURING READING FILE!"
                 exit
               else
                 exit
               endif
        enddo
        rewind(1)

        END SUBROUTINE  READRC




