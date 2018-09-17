        PROGRAM CALWEIGHT        
        

        !=====================================================
        ! IN THIS PROGRAM, L SHOULD MATCH SIZE WITH NODAL VALUE
        ! L SHOULD BE CHANGED;
        ! NODAL VALUE SHOULD BE CHANGED;
        !=====================================================



        IMPLICIT NONE        
        
        INTEGER,PARAMETER            :: L=7           !FOR OPTIM
        INTEGER,PARAMETER            :: PP=100        !FOR 5 RANDOM SEED
        INTEGER*8                    :: OPT           !FOR OPTIM
        INTEGER                      :: IRES          !FOR OPTIM
        DOUBLE PRECISION             :: MINF          !FOR OPTIM
        DOUBLE PRECISION             :: X(L)          !WEIGHT(L),FOR OPTIM
        DOUBLE PRECISION             :: RANDOM(PP)    !FOR 5 RANDOM SEED 

        EXTERNAL MYFUNC, MYCONSTRAINT
        INCLUDE 'nlopt.f'        

                
        !=====================================================
        !  MINIMIZE NORM |AX-b| TO GET THE WEIGHTS ADDING TO 1
        !=====================================================
        
        ! CREATE GLOBAL OPTIMIZATION PROBLEM OBJECT
        OPT = 0
        CALL NLO_CREATE( OPT,NLOPT_LN_COBYLA,L )
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
        CALL INIT_RANDOM_SEED(PP)!! the most tricky thing (work)
        ! CALL RANDOM_SEED() !! the most tricky thing (work)
        CALL RANDOM_NUMBER(RANDOM)
        
        X(1:L) = RANDOM(1:L)        
        
        WRITE(*,*) "-----------------------------------------"
        WRITE(*,*) "RANDOM INITIAL WEIGHTS(=100) ARE LISTED AS FOLLOWS:"
        WRITE(*,1004) RANDOM 
        WRITE(*,*) "RANDOM INITIAL WEIGHTS(<5) ARE LISTED AS FOLLOWS:"
        WRITE(*,1004) X

1004    FORMAT(F18.8,F18.8,F18.8,F18.8) !SHOULD BE CHANGED L=5
        
        ! WORKING ON OPTIMIZATION PROBLEM FOR ME!
        ! ires = 1.0 ! THIS IS SUPER IMPORTANT,otherwise,directly fail!!
        CALL NLO_OPTIMIZE( IRES,OPT,X,MINF)!change to local_opti still worked!!!!
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
        

        END PROGRAM CALWEIGHT


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
        
!        !TYPING IN OPTIMIZER,THE NODAL VALUES,STORM PARAMETER,CATEGORY1
!        NODAL=(/-0.01885646,  0.00000013, -1.30708724, -0.00003627,&
!                -0.01885637, -0.00000001,  1.30708725, -0.00003644,&
!                -0.56355349, -1.14705095, -0.00000005, -0.00011303,&
!                 1.59271948,  0.00000071, -0.00000018,  0.00042243,&
!                -0.56355473,  1.14705079, -0.00000006, -0.00011325,&
!                 2.17487358,  2.11974243,  2.73568275,  1.31883086,&
!                 0.29007506, -1.14608709,  0.91401634,  2.62814870,&
!                 2.97868675, -1.82072309, -1.89681372, -1.91967468,&
!                -0.73816532,  0.51502721,  1.02513569,  0.18806785,&
!                -2.85173522,  0.37982843, -1.52078679,  2.56061970,&
!                 0.58099215, -2.75823670,  0.17532410, -2.91280402,&
!                 0.37944737,  0.54133879,  0.36309603, -1.73980694,&
!                -2.45177028,  1.25562154,  2.38947096, -0.12557591/)


!        !TYPING IN OPTIMIZER,THE NODAL VALUES,STORM PARAMETER,CATEGORY2
!        NODAL=(/ 0.57639573,  1.14386851,  0.00000003,  0.00000494,&
!                 0.55077536, -1.15011738, -0.00000016,  0.00030604,&
!                -1.59254266,  0.01889689,  0.00000013, -0.00022354,&
!                 0.01884644, -0.00032124, -1.30708731, -0.00006575,&
!                 0.01884653, -0.00032118,  1.30708738, -0.00006584/)

        
!        !TYPING IN OPTIMIZER,THE NODAL VALUES,STORM PARAMETER,CATEGORY3
!        NODAL=(/ 0.56259111,  1.14728639,  0.00000021, -0.00000021,&
!                 0.01885606,  0.00002404,  1.30708719, -0.00000017,&
!                 0.01885610,  0.00002425, -1.30708713,  0.00000007,&
!                -1.59271698, -0.00142117,  0.00000005,  0.00000013,&
!                 0.56451997, -1.14681528, -0.00000003,  0.00000005,&
!                -1.61426997,  1.49155776,  1.37909018, -1.67559054,&
!                 1.51926185,  2.87686159, -0.52495761, -0.58271178/)


        !TYPING IN OPTIMIZER,THE NODAL VALUES,STORM PARAMETER,CATEGORY4
        NODAL=(/ 0.54715187, -1.15096919, -0.00000038,  0.00000016,&
                 0.58006318,  1.14294194, -0.00000001, -0.00000006,&
                 0.01883979, -0.00041217, -1.30708742, -0.00000019,&
                 0.01884004, -0.00041254,  1.30708732, -0.00000018,&
                -1.59242548,  0.02427223,  0.00000001,  0.00000024,&
                -2.97203553,  1.15282058,  0.59348701,  1.90770162,&
                 0.47498376, -2.11400852,  2.99167022, -0.18781675/)
        
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
!        WRITE(*,*)  "GIVE ME UD:"
!        WRITE(*,1004) UD  

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

!       WRITE(*,*)  "GIVE ME KORIGIN:"
!       WRITE(*,1004) KORIGIN(:,1),KORIGIN(:,2),KORIGIN(:,3),&
!                      KORIGIN(:,4),KORIGIN(:,5)   

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
                                                            
                       
