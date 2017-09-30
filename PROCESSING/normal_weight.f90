        PROGRAM CALWEIGHT        
        
        IMPLICIT NONE        
        
        INTEGER                      :: I,QQ          !FOR LOOP
        INTEGER                      :: NR            !FOR READ FILE
        INTEGER                      :: N             !FOR OPTIM
        INTEGER*8                    :: OPT           !FOR OPTIM
        INTEGER                      :: IRES          !FOR OPTIM
        DOUBLE PRECISION             :: MINF          !FOR OPTIM
        INTEGER,PARAMETER            :: PP=100        !FOR RANDOM SEED
        DOUBLE PRECISION             :: RANDOM(PP)    !FOR RANDOM SEED 
        DOUBLE PRECISION,ALLOCATABLE :: X(:)          !WEIGHT(NR),FOR OPTIM
        DOUBLE PRECISION,ALLOCATABLE :: NODAL(:)      !NODAL(NR*4),FOR OPTIM

        EXTERNAL MYFUNC, MYCONSTRAINT
        INCLUDE 'nlopt.f'        

        CALL READRC(NR)

        N = NR*4
        ALLOCATE(X(NR))
        ALLOCATE(NODAL(N))
        
        OPEN(UNIT=1,FILE="DATA_ORIGINAL.CSV")
        READ(1,1004) NODAL                
        CLOSE(1)
   
        !SHOW ORIGINAL DATA READING FROM INPUT FILE
        WRITE(*,1004) NODAL

        !=====================================================
        !  MINIMIZE NORM |AX-b| TO GET THE WEIGHTS ADDING TO 1
        !=====================================================
        
        ! CREATE GLOBAL OPTIMIZATION PROBLEM OBJECT
        OPT = 0
        CALL NLO_CREATE( OPT,NLOPT_LN_COBYLA,NR )
        IF (OPT.LE.0) THEN
           WRITE(*,*) "NLOPT FAILED AT CREATING GLOBAL OBJECT!"
        ELSE
           WRITE(*,*) "CREATING GLOBAL OBJECT SUCCESSFULLY!"
        ENDIF

        ! SET GLOBAL TOLERANCE
        CALL NLO_SET_XTOL_REL( IRES,OPT,1.D-7 )
        CALL NLO_SET_FTOL_REL( IRES,OPT,1.D-15)
        IF (IRES.LT.0) THEN
           WRITE(*,*) "NLOPT FAILED AT SETTING LOCAL TOLERANCE!"
        ELSE
           WRITE(*,*) "SETTING GLOBAL TOLERANCE SUCCESSFULLY!"
        ENDIF
        
        ! SET GLOBAL BOUND CONSTRAINS
        CALL NLO_SET_LOWER_BOUNDS1( IRES,OPT,0.0 )
        CALL NLO_SET_UPPER_BOUNDS1( IRES,OPT,1.0 )
        IF (IRES.LT.0) THEN
           WRITE(*,*) "NLOPT FAILED AT SETTING LOWER & UPPER BOUNDS!"
        ELSE
           WRITE(*,*) "SETTING GLOBAL LOWER&UPPER BOUNDS SUCCESSFULLY!"
        ENDIF
        
        ! SET GLOBAL CONSTRAINT FUNCTION
        CALL NLO_ADD_EQUALITY_CONSTRAINT(IRES,OPT,MYCONSTRAINT,0,1.D-15)
        IF (IRES.LT.0) THEN       
           WRITE(*,*) "NLOPT FAILED AT SETTING UP CONSTRAINT FUNCTION!"
        ELSE                      
           WRITE(*,*) "SETTING GLOBAL CONSTRAINT FUNCTION SUCCESSFULLY!"
        ENDIF                     
        
        ! SET GLOBAL OBJECTIVE FUNCTION
        CALL NLO_SET_MIN_OBJECTIVE( IRES,OPT,MYFUNC,0 )!we need both
        IF (IRES.LT.0) THEN       
           WRITE(*,*) "NLOPT FAILED AT SETTING UP OBJECTIVE FUNCTION!"
        ELSE                      
           WRITE(*,*) "SETTING GLOBAL OBJECTIVE FUNCTION SUCCESSFULLY!"
        ENDIF                     
        
        ! GIVEN INITIAL POINT FOR X(N)
        CALL INIT_RANDOM_SEED(PP)
        CALL RANDOM_NUMBER(RANDOM)
        
        X(1:NR) = RANDOM(1:NR)        
        
        WRITE(*,*) "-----------------------------------------"
        WRITE(*,*) "RANDOM INITIAL WEIGHTS(=100) ARE LISTED AS FOLLOWS:"
        WRITE(*,1004) RANDOM 
        WRITE(*,*) "RANDOM INITIAL WEIGHTS(<5) ARE LISTED AS FOLLOWS:"
        WRITE(*,1004) X

1004    FORMAT(F18.8,F18.8,F18.8,F18.8) 
        
        ! WORKING ON OPTIMIZATION PROBLEM FOR ME!
        CALL NLO_OPTIMIZE( IRES,OPT,X,MINF)
1005    FORMAT(A,I3,A)
        IF (IRES.EQ.-1) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," GENERIC FAILURE CODE!"
           WRITE(*,*) "FOUND GLOBAL WEIHGT OPTIMIZER AT: "
           WRITE(*,1004) X
           WRITE(*,*) "MINIMUM NORM OF OBJECTIVE FUNCVAL = ",MINF
        ELSEIF (IRES.EQ.-2) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," INVALID ARGUMENTS:LIKE LOWERBO&
                        UND ARE BIGGER THAN UPPER BOUND, OR UNKNOWN A&
                        LGORITHM WAS SPECIFIED!"
           WRITE(*,*) "FOUND GLOBAL WEIHGT OPTIMIZER AT: "
           WRITE(*,1004) X
           WRITE(*,*) "MINIMUM NORM OF OBJECTIVE FUNCVAL = ",MINF
        ELSEIF (IRES.EQ.-3) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," RAN OUT OF MEMORY!"
           WRITE(*,*) "FOUND GLOBAL WEIHGT OPTIMIZER AT: "
           WRITE(*,1004) X
           WRITE(*,*) "MINIMUM NORM OF OBJECTIVE FUNCVAL = ",MINF
        ELSEIF (IRES.EQ.-4) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," HALTED BECAUSE ROUNDOFF ERRORS&
                        LIMITED PROCESS(RESULT IS STILL TRUSTABLE)!"
           WRITE(*,*) "FOUND GLOBAL WEIHGT OPTIMIZER AT: "
           WRITE(*,1004) X
           WRITE(*,*) "MINIMUM NORM OF OBJECTIVE FUNCVAL = ",MINF
        ELSEIF (IRES.EQ.-5) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," HALTED BECAUSE OF A FORCED TER&
                        MINATION: USER CALLED NLO_FORCE_STOP(OPT)!"
           WRITE(*,*) "FOUND GLOBAL WEIHGT OPTIMIZER AT: "
           WRITE(*,1004) X
           WRITE(*,*) "MINIMUM NORM OF OBJECTIVE FUNCVAL = ",MINF
        ELSEIF (IRES.EQ.1) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," GENERIC SUCCESS RETURN VALUE!"
           WRITE(*,*) "FOUND GLOBAL WEIHGT OPTIMIZER AT: "
           WRITE(*,1004) X
           WRITE(*,*) "MINIMUM NORM OF OBJECTIVE FUNCVAL = ",MINF
        ELSEIF (IRES.EQ.2) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," OPTIMIZATION STOPPED BECAUSE STO&
                        PVAL WAS REACHED!"
           WRITE(*,*) "FOUND GLOBAL WEIHGT OPTIMIZER AT: "
           WRITE(*,1004) X
           WRITE(*,*) "MINIMUM NORM OF OBJECTIVE FUNCVAL = ",MINF
        ELSEIF (IRES.EQ.3) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," OPTIMIZATION STOPPED BECAUSE FTO&
                        L WAS REACHED!"
           WRITE(*,*) "FOUND GLOBAL WEIHGT OPTIMIZER AT: "
           WRITE(*,1004) X
           WRITE(*,*) "MINIMUM NORM OF OBJECTIVE FUNCVAL = ",MINF
        ELSEIF (IRES.EQ.4) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," OPTIMIZATION STOPPED BECAUSE XTO&
                        L WAS REACHED!"
           WRITE(*,*) "FOUND GLOBAL WEIHGT OPTIMIZER AT: "
           WRITE(*,1004) X
           WRITE(*,*) "MINIMUM NORM OF OBJECTIVE FUNCVAL = ",MINF
        ELSEIF (IRES.EQ.5) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," OPTIMIZATION STOPPED BECAUSE MAX&
                        EVAL WAS REACHED!"
           WRITE(*,*) "FOUND GLOBAL WEIHGT OPTIMIZER AT: "
           WRITE(*,1004) X
           WRITE(*,*) "MINIMUM NORM OF OBJECTIVE FUNCVAL = ",MINF
        ELSEIF (IRES.EQ.6) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," OPTIMIZATION STOPPED BECAUSE MAX&
                        TIME WAS REACHED!"
           WRITE(*,*) "FOUND GLOBAL WEIHGT OPTIMIZER AT: "
           WRITE(*,1004) X
           WRITE(*,*) "MINIMUM NORM OF OBJECTIVE FUNCVAL = ",MINF
        ELSE
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," OPTIMIZATION DOESN'T WORK ANY MO&
                        RE FOR YOU,BAD BOY!"
        ENDIF
        
        !ADD PROBABILIRY WEIGHTS INTO FILE
1003    FORMAT(F18.8,F18.8,F18.8,F18.8,F18.8) 
        OPEN( UNIT=2,FILE="DATA_WEIGHTS.CSV" )
        DO QQ=1,NR
           WRITE(2,1003) NODAL(QQ),NODAL(NR+QQ),NODAL(2*NR+QQ),&
                                NODAL(3*NR+QQ),X(QQ)
        ENDDO
        CLOSE(2)         
        
        DEALLOCATE(X)
        DEALLOCATE(NODAL)
        END PROGRAM CALWEIGHT


        SUBROUTINE READRC(NR)
        
        IMPLICIT NONE
        INTEGER,INTENT(OUT)  :: NR 
        INTEGER              :: I,STAT
        CHARACTER(80)        :: COLUMN

        !read the row of the sample
        NR = 0
        OPEN(UNIT=1,FILE="DATA_ORIGINAL.CSV")
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
        close(1)
       
        END SUBROUTINE  READRC


        SUBROUTINE INIT_RANDOM_SEED(N)
   
        IMPLICIT NONE
        INTEGER :: I, CLOCK
        INTEGER :: N
        INTEGER :: SEED(N)
   
        CALL SYSTEM_CLOCK(COUNT=CLOCK)
        SEED = CLOCK + 37 * (/ (I - 1, I = 1, N) /)
        CALL RANDOM_SEED(PUT = SEED)
   
        END SUBROUTINE INIT_RANDOM_SEED


        SUBROUTINE MYFUNC( FUNCNORM,L,X,GRAD,NEED_GRADIENT,F_DATA )

        IMPLICIT NONE

        INTEGER               :: I,J,K                    !FOR LOOP
        INTEGER               :: L                        !FOR OPTIM
        INTEGER               :: NEED_GRADIENT            !FOR WEIGHTS
        DOUBLE PRECISION      :: TEMP2,TEMP3              !FOR WEIGHTS
        DOUBLE PRECISION      :: FUNCNORM,VECDOT,F_DATA   !FOR WEIGHTS
        DOUBLE PRECISION      :: CORELATION(4)            !FOR WEIGHTS
        DOUBLE PRECISION      :: NODAL(4*L)               !FOR WEIGHTS
        DOUBLE PRECISION      :: UD(L),GRAD(L)            !FOR WEIGHTS
        DOUBLE PRECISION      :: X(L),VECTOR(L)           !FOR OPTIM
        DOUBLE PRECISION      :: KORIGIN(L,L)             !FOR WEIGHTS
        DOUBLE PRECISION,PARAMETER   :: PI = 3.1415926    !FOR WEIGHTS
        
        !=====================================================
        ! CALCULATING THE WEIGHTS FOR EXCEEDENCE PROBABILITY 
        !=====================================================
        
        !CALCULATE CORELATION
        CORELATION(1:4) =[5.0,4.0,2.5,6.0]/DSQRT(PI)! HEADING,PRESSURE,RADIUS,VFSPEED

        !CALCULATE UD=TEMP2
        DO I = 0,L-1 
        TEMP2= 1.0
        DO J = 1,4
        TEMP2 = TEMP2 * ( 1.0+2.0/CORELATION(J)**2 )**(-1.0/2.0) &
                      * EXP( -NODAL(4*I+J)**2/(CORELATION(J)**2+2.0) ) 
        ENDDO 
        UD(I+1) = TEMP2 
        ENDDO 

1004    FORMAT(F18.8,F18.8,F18.8,F18.8) !SHOULD BE CHANGED L=5

        !CALCULATE KORIGIN
        DO I = 0,L-1
        DO J = 0,L-1
        TEMP3= 1.0
        DO K = 1,4
          TEMP3 = TEMP3 * EXP(-((NODAL(4*I+K)-NODAL(4*J+K))/&
                                CORELATION(K))**2)
        ENDDO
        KORIGIN(I+1,J+1) = TEMP3
        ENDDO
        ENDDO

        VECTOR   = MATMUL( KORIGIN,X  ) - UD    
        VECDOT   = DOT_PRODUCT(VECTOR,VECTOR)
        FUNCNORM = DSQRT(VECDOT)
        

        WRITE(*,*) "-----------------------------------------"
        WRITE(*,*)  "THE WEIGHT DURING THIS ITERATION IS:"
        WRITE(*,1004) X
        WRITE(*,*)  "THE NORM OF OBJECTIVE FUNCTION |AX-b| IS THIS:"&
                        ,FUNCNORM
        

        END SUBROUTINE MYFUNC


        SUBROUTINE MYCONSTRAINT(CON,L,X,GRAD,NEED_GRADIENT,H_DATA)

        INTEGER          :: NEED_GRADIENT
        DOUBLE PRECISION :: CON,H_DATA
        DOUBLE PRECISION :: X(L),GRAD(L)
        
        CON = ABS( SUM(X(1:L)) - 1 )

        WRITE(*,*)  "THE CONSTRAINT ERROR DURING THIS ITERATION IS:",CON
        
          END SUBROUTINE MYCONSTRAINT
                                                            
                       
