
      PROGRAM  main
        
        IMPLICIT NONE       
        !=====================================================
        !=================  DATA DECLARATION  ================
        !=====================================================
        INTEGER,PARAMETER            :: PP=100            !FOR RANDOM SEED
        INTEGER                      :: N,NR,QQ            !FOR #OPTIMPARA & #STORMS 
        CHARACTER(80)                :: OUTPUTNAME        !FOR THE NAME OF OUTPUT FILE
        DOUBLE PRECISION             :: RANDOM(PP)        !FOR RANDOM SEED    

        INTEGER(4)                   :: F_DATA            !FOR OPTIM : LIBRARY FIXED
        INTEGER(8)                   :: IRES              !FOR OPTIM : LIBRARY FIXED
        INTEGER*8                    :: OPT               !FOR OPTIM : LIBRARY FIXED
        INTEGER*8                    :: LOCAL_OPT         !FOR OPTIM : LIBRARY FIXED
        DOUBLE PRECISION             :: MINF              !FOR OPTIM : LIBRARY FIXED

        DOUBLE PRECISION,PARAMETER   :: PI = 3.1415926    !FOR MYFUNC    
        DOUBLE PRECISION,ALLOCATABLE :: X(:)              !FOR MYFUNC   
     
        EXTERNAL MYFUNC
        INCLUDE 'nlopt.f'        

        !=====================================================
        !======== INPUT NUMBER OF STORMS YOU WANT GENERATE ===
        !=====================================================
        WRITE(*,*) "ENTER THE NUMBER OF SYNTHETIC STORMS YOU WANT&
            TO GENERATE(PLEASE TYPE IN AN INTEGER FROM(1:13)):"
        READ(*,*)  NR ! #GOALNODES,13,5,7,7
!        WRITE(*,*) "ENTER THE NAME OF OUTPUT FILE(PLEASE ENTER &
!                    A STRING.csv):"
!        READ(*,*)  OUTPUTNAME

        !=====================================================
        !======== NODAL VALUE OPTIMIZATION ===================
        !=====================================================
        
        ! SETTING OPTIMPARAS N FOR THE USE OF NLOPT LIBRARY
        N = 4*NR   
        ALLOCATE(X(N))

        ! CREATE GLOBAL OPTIMIZATION PROBLEM OBJECT
        OPT = 0
        CALL NLO_CREATE( OPT,NLOPT_AUGLAG,N )
        IF (OPT.LE.0) THEN
           WRITE(*,*) "NLOPT FAILED AT CREATING GLOBAL OBJECT!"
        ELSE
           WRITE(*,*) "CREATING GLOBAL OBJECT SUCCESSFULLY!"
        ENDIF

        ! CREATE LOCAL OPTIMIZATION PROBLEM OBJECT
        LOCAL_OPT = 1        
        CALL NLO_CREATE( LOCAL_OPT,NLOPT_LN_COBYLA,N )
        IF (LOCAL_OPT.LE.0) THEN
           WRITE(*,*) "NLOPT FAILED AT CREATING LOCAL OBJECT!"
        ELSE
           WRITE(*,*) "CREATING LOCAL OBJECT SUCCESSFULLY!"
        ENDIF

        ! SET GLOBAL TOLERANCE
        CALL NLO_SET_XTOL_REL( IRES,OPT,1.D-4 )
        CALL NLO_SET_FTOL_REL( IRES,OPT,1.D-8 )
        CALL NLO_SET_MAXTIME ( IRES,OPT,3600 )
        IF (IRES.LT.0) THEN
           WRITE(*,*) "NLOPT FAILED AT SETTING LOCAL TOLERANCE!"
        ELSE
           WRITE(*,*) "SETTING GLOBAL TOLERANCE SUCCESSFULLY!"
        ENDIF

        ! SET LOCAL TOLERANCE
        CALL NLO_SET_XTOL_REL( IRES,LOCAL_OPT,1.D-7 )
        CALL NLO_SET_FTOL_REL( IRES,LOCAL_OPT,1.D-15)
        IF (IRES.LT.0) THEN
           WRITE(*,*) "NLOPT FAILED AT SETTING LOCAL TOLERANCE!"
        ELSE
           WRITE(*,*) "SETTING LOCAL TOLERANCE SUCCESSFULLY!"
        ENDIF
        
        ! SET LOCAL OBJECTIVE FUNCTION
        F_DATA = 0
        CALL NLO_SET_MIN_OBJECTIVE( IRES,LOCAL_OPT,MYFUNC,F_DATA )!we need both
        IF (IRES.LT.0) THEN       
           WRITE(*,*) "NLOPT FAILED AT SETTING UP OBJECTIVE FUNCTION!"
        ELSE                      
           WRITE(*,*) "SETTING LOCAL OBJECTIVE FUNCTION SUCCESSFULLY!"
        ENDIF                     
        
        ! SET GLOBAL OBJECTIVE FUNCTION
        CALL NLO_SET_MIN_OBJECTIVE( IRES,OPT,MYFUNC,F_DATA )!we need both
        IF (IRES.LT.0) THEN       
           WRITE(*,*) "NLOPT FAILED AT SETTING UP OBJECTIVE FUNCTION!"
        ELSE                      
           WRITE(*,*) "SETTING GLOBAL OBJECTIVE FUNCTION SUCCESSFULLY!"
        ENDIF                     
        
        ! CONNECT LOCAL ALGORITHM WITH GLOBAL ALGORITHM 
        CALL NLO_SET_LOCAL_OPTIMIZER( IRES,OPT,LOCAL_OPT )
        IF (IRES.LT.0) THEN       
           WRITE(*,*) "NLOPT FAILED AT SETTING LOCAL OPTIMIZER!"
        ELSE                      
           WRITE(*,*) "CONNECTING LOCAL AND GLOBAL SUCCESSFULLY!"
        ENDIF                     
        
        ! GIVEN INITIAL POINT FOR X(N)
        CALL INIT_RANDOM_SEED(PP) 
        CALL RANDOM_NUMBER(RANDOM)        
        
        X(1:N) = 6.0*RANDOM(1:N)-3.0 

1003    FORMAT(F18.8,F18.8,F18.8,F18.8)
        WRITE(*,*) "-----------------------------------------"
        WRITE(*,*) "RANDOM INITIAL POINTS ARE LISTED AS FOLLOWS:"
        WRITE(*,1003) X

        ! WORKING ON OPTIMIZATION PROBLEM FOR ME!
        IRES = 1.0 ! THIS INITIALIZATION SUPER IMPORTANT,otherwise,directly fail!!
        CALL NLO_OPTIMIZE( IRES,LOCAL_OPT,X,MINF)!change to opt also worked!!!!
1005    FORMAT(A,I3,A)
        IF (IRES.EQ.-1) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," GENERIC FAILURE CODE!"
           WRITE(*,*) "FOUND GLOBAL MINIMUM AT: "
           WRITE(*,1003) X
           WRITE(*,*) "MIN OBJECTIVE FUNCVAL = ",MINF
           WRITE(*,*) "NLOPT HAS DONE",F_DATA,"ITERATION TIMES FOR YOU!" 
        ELSEIF (IRES.EQ.-2) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," INVALID ARGUMENTS:LIKE LOWERBO&
                        UND ARE BIGGER THAN UPPER BOUND, OR UNKNOWN A&
                        LGORITHM WAS SPECIFIED!"
           WRITE(*,*) "FOUND GLOBAL MINIMUM AT: "
           WRITE(*,1003) X
           WRITE(*,*) "MIN OBJECTIVE FUNCVAL = ",MINF
           WRITE(*,*) "NLOPT HAS DONE",F_DATA,"ITERATION TIMES FOR YOU!" 
        ELSEIF (IRES.EQ.-3) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," RAN OUT OF MEMORY!"
           WRITE(*,*) "FOUND GLOBAL MINIMUM AT: "
           WRITE(*,1003) X
           WRITE(*,*) "MIN OBJECTIVE FUNCVAL = ",MINF
           WRITE(*,*) "NLOPT HAS DONE",F_DATA,"ITERATION TIMES FOR YOU!" 
        ELSEIF (IRES.EQ.-4) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," HALTED BECAUSE ROUNDOFF ERRORS&
                        LIMITED PROCESS(RESULT IS STILL TRUSTABLE)!"
           WRITE(*,*) "FOUND GLOBAL MINIMUM AT: "
           WRITE(*,1003) X
           WRITE(*,*) "MIN OBJECTIVE FUNCVAL = ",MINF
           WRITE(*,*) "NLOPT HAS DONE",F_DATA,"ITERATION TIMES FOR YOU!" 
        ELSEIF (IRES.EQ.-5) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," HALTED BECAUSE OF A FORCED TER&
                        MINATION: USER CALLED NLO_FORCE_STOP(OPT)!"
           WRITE(*,*) "FOUND GLOBAL MINIMUM AT: "
           WRITE(*,1003) X
           WRITE(*,*) "MIN OBJECTIVE FUNCVAL = ",MINF
           WRITE(*,*) "NLOPT HAS DONE",F_DATA,"ITERATION TIMES FOR YOU!" 
        ELSEIF (IRES.EQ.1) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," GENERIC SUCCESS RETURN VALUE!"
           WRITE(*,*) "FOUND GLOBAL MINIMUM AT: "
           WRITE(*,1003) X
           WRITE(*,*) "MIN OBJECTIVE FUNCVAL = ",MINF
           WRITE(*,*) "NLOPT HAS DONE",F_DATA,"ITERATION TIMES FOR YOU!" 
        ELSEIF (IRES.EQ.2) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," OPTIMIZATION STOPPED BECAUSE STO&
                        PVAL WAS REACHED!"
           WRITE(*,*) "FOUND GLOBAL MINIMUM AT: "
           WRITE(*,1003) X
           WRITE(*,*) "MIN OBJECTIVE FUNCVAL = ",MINF
           WRITE(*,*) "NLOPT HAS DONE",F_DATA,"ITERATION TIMES FOR YOU!" 
        ELSEIF (IRES.EQ.3) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," OPTIMIZATION STOPPED BECAUSE FTO&
                       L WAS REACHED!"
           WRITE(*,*) "FOUND GLOBAL MINIMUM AT: "
           WRITE(*,1003) X
           WRITE(*,*) "MIN OBJECTIVE FUNCVAL = ",MINF
           WRITE(*,*) "NLOPT HAS DONE",F_DATA,"ITERATION TIMES FOR YOU!" 
        ELSEIF (IRES.EQ.4) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," OPTIMIZATION STOPPED BECAUSE XTO&
                        L WAS REACHED!"
           WRITE(*,*) "FOUND GLOBAL MINIMUM AT: "
           WRITE(*,1003) X
           WRITE(*,*) "MIN OBJECTIVE FUNCVAL = ",MINF
           WRITE(*,*) "NLOPT HAS DONE",F_DATA,"ITERATION TIMES FOR YOU!" 
        ELSEIF (IRES.EQ.5) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," OPTIMIZATION STOPPED BECAUSE MAX&
                        EVAL WAS REACHED!"
           WRITE(*,*) "FOUND GLOBAL MINIMUM AT: "
           WRITE(*,1003) X
           WRITE(*,*) "MIN OBJECTIVE FUNCVAL = ",MINF
           WRITE(*,*) "NLOPT HAS DONE",F_DATA,"ITERATION TIMES FOR YOU!" 
        ELSEIF (IRES.EQ.6) THEN
           WRITE(*,*) "----------------------------------------------"
           WRITE(*,1005) "IRES = ",IRES," OPTIMIZATION STOPPED BECAUSE MAX&
                        TIME WAS REACHED!"
           WRITE(*,*) "FOUND GLOBAL MINIMUM AT: "
           WRITE(*,1003) X
           WRITE(*,*) "MIN OBJECTIVE FUNCVAL = ",MINF
           WRITE(*,*) "NLOPT HAS DONE",F_DATA,"ITERATION TIMES FOR YOU!" 
        ELSE
           WRITE(*,*) IRES
           WRITE(*,*) "NLO_OPTIMIZE DOESN'T WORK ANY MORE FOR YOU,&
                        BAD BOY!"
        ENDIF
        
        !WRITE THE ORIGINAL RESULT INTO FILE 
        OPEN( UNIT=2,FILE="DATA_ORIGINAL.CSV" )
        WRITE(2,1003) X
        CLOSE(2)         

        !SORTING RESULT X(N) INTO ACSENDING ORDER
        CALL SORT(N,X)
 
        !WRITE THE RESULT INTO FILE 
        OPEN( UNIT=1,FILE="DATA_ORDERING.CSV" )
        DO QQ=1,NR
           WRITE(1,1003) X(QQ),X(NR+QQ),X(2*NR+QQ),X(3*NR+QQ)
        ENDDO
        CLOSE(1)         
        
        END PROGRAM main
        

        !================================================
        ! GENERATING INITIAL RANDOM SEED FOR OPTIMIZATION 
        !================================================
        SUBROUTINE INIT_RANDOM_SEED(N)
          
        IMPLICIT NONE
        INTEGER :: I, CLOCK
        INTEGER :: N
        INTEGER :: SEED(N)
        
        CALL SYSTEM_CLOCK(COUNT=CLOCK)
        SEED = CLOCK + 37 * (/ (I - 1, I = 1, N) /)
        CALL RANDOM_SEED(PUT = SEED)

        END SUBROUTINE INIT_RANDOM_SEED

         
        !================================================
        ! COST MYFUNC FOR GLOBAL OPTIMIZATION 
        !================================================
        SUBROUTINE MYFUNC(VAL,N,X,GRAD,NEED_GRADIENT,F_DATA) 
          
        INTEGER                       :: N  ! #OPTINPARA,52,20,28,28    
        INTEGER                       :: L  ! #GOALNODES,13,5,7,7     
        INTEGER                       :: F_DATA ! #ITERATION TIMES     
        INTEGER                       :: I,J,K,NEED_GRADIENT
        DOUBLE PRECISION              :: TEMP1
        DOUBLE PRECISION              :: TEMP2
        DOUBLE PRECISION              :: TEMP3
        DOUBLE PRECISION              :: VAL
        DOUBLE PRECISION              :: X(N),GRAD(N) 
        DOUBLE PRECISION              :: CORELATION(4),DOTPRO
        DOUBLE PRECISION,ALLOCATABLE  :: KORIGIN(:,:)
        DOUBLE PRECISION,ALLOCATABLE  :: KINVERSE(:,:) 
        DOUBLE PRECISION,ALLOCATABLE  :: UD(:),MAXVEC(:)
        DOUBLE PRECISION,PARAMETER    :: PI = 3.1415926   

        L = N / 4      
        
        ALLOCATE( KORIGIN(L,L) )        
        ALLOCATE( KINVERSE(L,L) )        
        ALLOCATE( UD(L) )        
        ALLOCATE( MAXVEC(L) )        
        
        !INITIALIZING GRAD(N)
        IF (NEED_GRADIENT.NE.0) THEN
           GRAD(1:N)=0.0
        ENDIF

        !CALCULATE CORELATION
        CORELATION(1:4) =[5.0,4.0,2.5,6.0]/DSQRT(PI)! HEADING,PRESSURE,RADIUS,VFSPEED

        !CALCULATE U=TEMP1
        TEMP1 = 1.0
        DO I = 1,4
        TEMP1 = TEMP1 / DSQRT(1.0+4.0/CORELATION(I)**2)
        ENDDO

        !CALCULATE UD=TEMP2
        DO I = 0,L-1 
        TEMP2= 1.0 
        DO J = 1,4
        TEMP2 = TEMP2 * ( 1.0+2.0/CORELATION(J)**2 )**(-1.0/2.0) &
                      * EXP( -X(4*I+J)**2/(CORELATION(J)**2+2.0) ) 
        ENDDO 
        UD(I+1) = TEMP2 
        ENDDO 

        !CALCULATE KORIGIN
        DO I = 0,L-1
        DO J = 0,L-1
        TEMP3= 1.0
        DO K = 1,4
        TEMP3 = TEMP3 * EXP(-((X(4*I+K)-X(4*J+K))/CORELATION(K))**2)
        ENDDO
        KORIGIN(I+1,J+1) = TEMP3
        ENDDO
        ENDDO

        !CALCULATE KINVERSE
        CALL INVERSE(KINVERSE,KORIGIN,L)
        
        MAXVEC =  MATMUL(KINVERSE,UD)
        DOTPRO =  DOT_PRODUCT(UD,MAXVEC)

        VAL = TEMP1 - DOTPRO
        F_DATA = F_DATA + 1

        WRITE(*,*) "FUNCTION VALUE AT OPTIMAL POINT IS:",VAL     

        DEALLOCATE( KORIGIN )        
        DEALLOCATE( KINVERSE )        
        DEALLOCATE( UD )        
        DEALLOCATE( MAXVEC )        

        END SUBROUTINE MYFUNC 


        !================================================
        ! INVERSE OF A MATRIX FOR GLOBAL OPTIMIZATION 
        !================================================
        SUBROUTINE INVERSE(KINVERSE,KORIGIN,M)

        !======= THE ORIGINAL K(M,M) WILL BE DESTROYED
        IMPLICIT NONE
        INTEGER(4),INTENT(IN) :: M
        DOUBLE PRECISION,INTENT(INOUT)  :: KORIGIN(M,M)
        DOUBLE PRECISION,INTENT(OUT) :: KINVERSE(M,M)
        DOUBLE PRECISION :: L(M,M),U(M,M),B(M),D(M),X(M)
        DOUBLE PRECISION :: COEFF
        INTEGER :: II,JJ,KK

        !======= INITIALIZATION FOR MATRICE U & L & B
        L = 0.0
        U = 0.0
        B = 0.0
        
        !======= FORWARD ELIMINATION
        DO KK = 1,M-1
        DO II = KK+1,M
        COEFF = KORIGIN(II,KK)/KORIGIN(KK,KK)
        L(II,KK) = COEFF
        DO JJ = KK+1,M
        KORIGIN(II,JJ) = KORIGIN(II,JJ) - COEFF*KORIGIN(KK,JJ)
        ENDDO
        ENDDO
        ENDDO
        DO II = 1,M
        L(II,II) = 1.0
        ENDDO

        !====== U IS UPPER TRIANGULAR PART OF K(M,M)
        DO JJ = 1,M
        DO II = 1,JJ
        U(II,JJ) = KORIGIN(II,JJ)
        ENDDO
        ENDDO

        !====== COMPUTE COLUMN OF KINVERSE(M,M)
        DO KK = 1,M
         B(KK) = 1.0
         D(1)  = B(1)
         !===== (A)SOLVE LD=B USING FORWARD SUBSTITUTION
         DO II = 2,M
           D(II) = B(II)
           DO JJ = 1,II-1
              D(II) = D(II) - L(II,JJ)*D(JJ)
           ENDDO
         ENDDO
         !===== (B)SOLVE UX=D USING BACK SUBTITUTION
         X(M) = D(M)/U(M,M)
         DO II = M-1,1,-1
            X(II) = D(II)
            DO JJ = M,II+1,-1
               X(II) = X(II) - U(II,JJ)*X(JJ)
            ENDDO
            X(II) = X(II)/U(II,II)
         ENDDO
         !===== (C)FILL THE SOLUTIONS X(M) INTO COLUMN KK OF KINVERSE
         DO II = 1,M
            KINVERSE(II,KK) = X(II)
         ENDDO
         B(KK) = 0.0
        ENDDO

        END SUBROUTINE INVERSE
       
        
        !================================================
        ! REARRANGE X(N) IN ASCENDING ORDER  
        !================================================
        FUNCTION  FindMinimum(x, Start, ActualSize)
           IMPLICIT  NONE
           INTEGER, INTENT(IN)                :: Start, ActualSize
           DOUBLE PRECISION, INTENT(IN)       :: x(ActualSize)
           DOUBLE PRECISION                   :: Minimum
           INTEGER                            :: Location
           INTEGER                            :: i
           INTEGER                            :: FindMinimum      
           Minimum  = x(Start)               ! assume the first is the min
           Location = Start                  ! record its position
           DO i = Start+1, ActualSize               ! start with next elements
              IF (x(i) < Minimum) THEN       !   if x(i) less than the min?
                 Minimum  = x(i)             !      Yes, a new minimum found
                 Location = i                !      record its position
              END IF
           END DO
           FindMinimum = Location            ! return the position
        END FUNCTION  FindMinimum
        
        SUBROUTINE  Swap(a, b)
           IMPLICIT  NONE
           DOUBLE PRECISION, INTENT(INOUT) :: a, b
           DOUBLE PRECISION                :: Temp
           Temp = a
           a    = b
           b    = Temp
        END SUBROUTINE  Swap
        
        SUBROUTINE  Sort(ActualSize,x)
           IMPLICIT  NONE
           INTEGER, INTENT(IN)                   :: ActualSize
           DOUBLE PRECISION, INTENT(INOUT)       :: x(ActualSize)
           INTEGER                               :: i
           INTEGER                               :: Location
           INTEGER,EXTERNAL                      :: FindMinimum     
           DO i = 1, ActualSize-1                  ! except for the last
              Location = FindMinimum(x, i, ActualSize)     ! find min from this to last
              CALL  Swap(x(i), x(Location))  ! swap this and the minimum
           END DO
        END SUBROUTINE  Sort


