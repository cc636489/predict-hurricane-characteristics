      PROGRAM  main
        
        IMPLICIT NONE       
        !=====================================================
        !=================  DATA DECLARATION  ================
        !=====================================================
        
        INTEGER(2)                   :: I,J,K,PP          !FOR LOOP
        INTEGER(2)                   :: NR,NC,OPENSTATUS  !FOR READ FILE
        INTEGER(2),DIMENSION(4)      :: P(1:4)=0,Q(1:4)=1 !FOR READ FILE
        CHARACTER(80)                :: FILENAME          !FOR READ FILE
        DOUBLE PRECISION,ALLOCATABLE :: SAMPLE(:)         !FOR READ FILE

        INTEGER(4)                   :: N,L               !FOR OPTIM : LIBRARY FIXED
        INTEGER(4)                   :: F_DATA            !FOR OPTIM : LIBRARY FIXED
        INTEGER(8)                   :: IRES              !FOR OPTIM : LIBRARY FIXED
        INTEGER*8                    :: OPT               !FOR OPTIM : LIBRARY FIXED
        INTEGER*8                    :: LOCAL_OPT         !FOR OPTIM : LIBRARY FIXED
        DOUBLE PRECISION             :: MINF              !FOR OPTIM : LIBRARY FIXED
        DOUBLE PRECISION             :: TEMP2,TEMP3       !FOR OPTIM : LIBRARY FIXED

        DOUBLE PRECISION             :: CORELATION(4)     !FOR WEIGHTS
        DOUBLE PRECISION             :: MIU               !FOR WEIGHTS
        DOUBLE PRECISION,PARAMETER   :: PI = 3.1415926    !FOR WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: X(:)              !FOR WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: UD(:)             !FOR WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: WEIGHT(:)         !FOR WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: KORIGIN(:,:)      !FOR WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: KINVERSE(:,:)     !FOR WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: UD_AUG(:)         !FOR WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: WEIGHT_AUG(:)     !FOR WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: KORIGIN_AUG(:,:)  !FOR WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: KINVERSE_AUG(:,:) !FOR WEIGHTS
      
        TYPE STORM_CATEGORY
                CHARACTER(80)        :: FULLNAME
                INTEGER(2)           :: CATENNODE 
                INTEGER(2)           :: GOALNNODE 
                REAL(8)              :: CATEPROBA
                REAL(8),ALLOCATABLE  :: GOALPROBA(:)  !WAIT TO CALCULATE
                REAL(8),ALLOCATABLE  :: CATEPARA(:,:)
                REAL(8),ALLOCATABLE  :: GOALPARA(:,:) !WAIT TO OPTIM
        END TYPE STORM_CATEGORY
        TYPE(STORM_CATEGORY) :: CATEGORY(4)
       
        EXTERNAL MYFUNC
        INCLUDE 'nlopt.f'        

        !=====================================================
        !=========  IMPORT AND REARRANGE DATA ================
        !=====================================================
       
        !OPEN FILE
        WRITE(*,*) "ENTER AN FILE NAME WITH .CSV EXTENSION: "
        READ(*,*) FILENAME
        OPEN(UNIT=1,FILE=FILENAME,IOSTAT=OPENSTATUS)
        
        !READ ROWS TO NR && COLUMNS TO NC
        CALL READRC(NR,NC)   
        ALLOCATE(SAMPLE(NC))
        
        DO I = 1,NR
        READ(1,*) SAMPLE(1:NC)
        IF ( (SAMPLE(4).GE.0) .AND. (SAMPLE(4).LT.5) ) THEN
           WRITE(*,*) 'WARNING:THIS PARTICULAR STORM PRESSURE IS TOO &
                       LOW TO BE CONSIDERED! WE DISCARD IT!' 
           CONTINUE

        ELSEIF ( (SAMPLE(4).GE.5) .AND. (SAMPLE(4).LT.45) ) THEN
           P(1) = P(1) + 1
        ELSEIF ( (SAMPLE(4).GE.45) .AND. (SAMPLE(4).LT.70) ) THEN
           P(2) = P(2) + 1
        ELSEIF ( (SAMPLE(4).GE.70) .AND. (SAMPLE(4).LT.95) ) THEN
           P(3) = P(3) + 1
        ELSEIF ( (SAMPLE(4).GE.95) ) THEN
           P(4) = P(4) + 1
        ELSE
           WRITE(*,*) 'ERROR: ORIGINAL DATA CONTAINS NEGATIVE PRESSURE!'
           CALL EXIT 
        ENDIF
        ENDDO
        REWIND(1)
        
        CATEGORY(1)%FULLNAME  = 'LESSER-STORM'
        CATEGORY(2)%FULLNAME  = 'GREATER-STORM-SIMPHSON 2'
        CATEGORY(3)%FULLNAME  = 'GREATER-STORM-SIMPHSON 3'
        CATEGORY(4)%FULLNAME  = 'GREATER-STORM-SIMPHSON 4'
        
        CATEGORY(1)%CATENNODE = P(1)
        CATEGORY(2)%CATENNODE = P(2)
        CATEGORY(3)%CATENNODE = P(3)
        CATEGORY(4)%CATENNODE = P(4)
        
        CATEGORY(1)%GOALNNODE = 13
        CATEGORY(2)%GOALNNODE = 5
        CATEGORY(3)%GOALNNODE = 7
        CATEGORY(4)%GOALNNODE = 7
        
        CATEGORY(1)%CATEPROBA = REAL(P(1))/P(1)
        CATEGORY(2)%CATEPROBA = REAL(P(2))/SUM(P(2:4))
        CATEGORY(3)%CATEPROBA = REAL(P(3))/SUM(P(2:4))
        CATEGORY(4)%CATEPROBA = REAL(P(4))/SUM(P(2:4))
        
        ALLOCATE( CATEGORY(1)%GOALPROBA(CATEGORY(1)%GOALNNODE) )
        ALLOCATE( CATEGORY(2)%GOALPROBA(CATEGORY(2)%GOALNNODE) )
        ALLOCATE( CATEGORY(3)%GOALPROBA(CATEGORY(3)%GOALNNODE) )
        ALLOCATE( CATEGORY(4)%GOALPROBA(CATEGORY(4)%GOALNNODE) )
       
        ALLOCATE( CATEGORY(1)%CATEPARA(P(1),NC) )
        ALLOCATE( CATEGORY(2)%CATEPARA(P(2),NC) )
        ALLOCATE( CATEGORY(3)%CATEPARA(P(3),NC) )
        ALLOCATE( CATEGORY(4)%CATEPARA(P(4),NC) )
        
        ALLOCATE( CATEGORY(1)%GOALPARA(CATEGORY(1)%GOALNNODE,NC) )
        ALLOCATE( CATEGORY(2)%GOALPARA(CATEGORY(2)%GOALNNODE,NC) )
        ALLOCATE( CATEGORY(3)%GOALPARA(CATEGORY(3)%GOALNNODE,NC) )
        ALLOCATE( CATEGORY(4)%GOALPARA(CATEGORY(4)%GOALNNODE,NC) )
        
        DO I = 1,NR
        READ(1,*) SAMPLE(1:NC)
        IF ( (SAMPLE(4).GE.5) .AND. (SAMPLE(4).LT.45) ) THEN
           CATEGORY(1)%CATEPARA(Q(1),:) = SAMPLE(:)
           Q(1) = Q(1) + 1
        ELSEIF ( (SAMPLE(4).GE.45) .AND. (SAMPLE(4).LT.70) ) THEN
           CATEGORY(2)%CATEPARA(Q(2),:) = SAMPLE(:)
           Q(2) = Q(2) + 1
        ELSEIF  ( (SAMPLE(4).GE.70) .AND. (SAMPLE(4).LT.95) ) THEN
           CATEGORY(3)%CATEPARA(Q(3),:) = SAMPLE(:)
           Q(3) = Q(3) + 1
        ELSEIF  ( (SAMPLE(4).GE.95) ) THEN
           CATEGORY(4)%CATEPARA(Q(4),:) = SAMPLE(:)
           Q(4) = Q(4) + 1
        ELSE
           CONTINUE    
        ENDIF
        ENDDO
        REWIND(1)

        !CHECK DATA IMPORT
        IF ((SUM(P(:)).EQ.SUM(Q(:))-4) .AND. (SUM(P(:)).LE.NR)) THEN
           WRITE(*,*) "CONGRATS: IMPORT & PARTITION DATA SUCCESSFULLY!"
        ELSE
           WRITE(*,*) "ERROR: #IMPORT DATA DON'T MATCH ORIGNIAL ONE!"
        ENDIF
        CLOSE(1)       

1001    FORMAT (A,I3,A,F11.3,F11.3,F11.3,F11.3,F11.3,F11.3)
1002    FORMAT (A,I3,A,A)
        DO I = 1,4
          WRITE(*,*) "-----------------------------------------"
          WRITE(*,1002) "FOR CATAGORY",I," : ",CATEGORY(I)%FULLNAME
          WRITE(*,*) "                YEAR.","   DISTANCE","    HEAD&
                ING","   PRESSURE","     RADIUS","    VFSPEED"
          DO J = 1,P(I)
          WRITE(*,1001) "LINE ",J," : ",CATEGORY(I)%CATEPARA(J,1:NC)
          ENDDO
        ENDDO
       
        DEALLOCATE(SAMPLE)
        
        !=====================================================
        !======== NODAL VALUE OPTIMIZATION (CASE PP=2) =======
        !=====================================================
        
        PP = 1  
        
        ! SETTING OPTIMPARAS N FOR THE USE OF NLOPT LIBRARY
        L = CATEGORY(PP)%GOALNNODE ! #GOALNODES,13,5,7,7    
        N = (NC-2)*L               ! #OPTINPARA,52,20,28,28
        ALLOCATE(X(N))
        
        ! CREATE GLOBAL OPTIMIZATION PROBLEM OBJECT
        OPT = 0
        CALL NLO_CREATE( OPT,NLOPT_AUGLAG,20 )
        IF (OPT.LE.0) THEN
           WRITE(*,*) "NLOPT FAILED AT CREATING GLOBAL OBJECT!"
        ELSE
           WRITE(*,*) "CREATING GLOBAL OBJECT SUCCESSFULLY!"
        ENDIF

        ! CREATE LOCAL OPTIMIZATION PROBLEM OBJECT
        LOCAL_OPT = 1        
        CALL NLO_CREATE( LOCAL_OPT,NLOPT_LN_COBYLA,20 )
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
        
!        ! SET GLOBAL BOUND CONSTRAINS
!        CALL NLO_SET_LOWER_BOUNDS1( IRES,OPT,-3.0 )
!        CALL NLO_SET_UPPER_BOUNDS1( IRES,OPT,3.0 )
!        IF (IRES.LT.0) THEN
!           WRITE(*,*) "NLOPT FAILED AT SETTING LOWER & UPPER BOUNDS!"
!        ELSE
!           WRITE(*,*) "SETTING GLOBAL LOWER&UPPER BOUNDS SUCCESSFULLY!"
!        ENDIF
!        
!        ! SET LOCAL BOUND CONSTRAINS
!        CALL NLO_SET_LOWER_BOUNDS1( IRES,LOCAL_OPT,-3.0 )
!        CALL NLO_SET_UPPER_BOUNDS1( IRES,LOCAL_OPT,3.0 )
!        IF (IRES.LT.0) THEN
!           WRITE(*,*) "NLOPT FAILED AT SETTING LOWER & UPPER BOUNDS!"
!        ELSE
!           WRITE(*,*) "SETTING LOCAL LOWER&UPPER BOUNDS SUCCESSFULLY!"
!        ENDIF
        
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
        CALL INIT_RANDOM_SEED(N) !! the most tricky thing (work)
        ! CALL RANDOM_SEED()!! random seed is valid but don't work don't
        ! kow why !!!!!!!
        CALL RANDOM_NUMBER(X)        
        
        X(:) = 6.0*X(:)-3.0 

1003    FORMAT(F18.8,F18.8,F18.8,F18.8)
        WRITE(*,*) "-----------------------------------------"
        WRITE(*,*) "RANDOM INITIAL POINTS ARE LISTED AS FOLLOWS:"
        WRITE(*,1003) X

        ! WORKING ON OPTIMIZATION PROBLEM FOR ME!
        ires = 1.0 ! THIS IS SUPER IMPORTANT,otherwise,directly fail!!
        CALL NLO_OPTIMIZE( IRES,LOCAL_OPT,X,MINF)!change to local_opt worked!!!!
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
        



      END PROGRAM main 



      SUBROUTINE INIT_RANDOM_SEED(N)
        
      IMPLICIT NONE
      INTEGER :: I, CLOCK
      INTEGER :: N
      INTEGER :: SEED(N)
      
      CALL SYSTEM_CLOCK(COUNT=CLOCK)
      SEED = CLOCK + 37 * (/ (I - 1, I = 1, N) /)
      CALL RANDOM_SEED(PUT = SEED)

        
      END SUBROUTINE INIT_RANDOM_SEED

       
 
      SUBROUTINE MYFUNC(VAL,N,X,GRAD,NEED_GRADIENT,f_data) 
        
        INTEGER :: N  ! #OPTINPARA,52,20,28,28    
        INTEGER :: L  ! #GOALNODES,13,5,7,7     
        INTEGER :: F_DATA ! #ITERATION TIMES     
        INTEGER :: I,J,K,NEED_GRADIENT
        DOUBLE PRECISION :: TEMP1
        DOUBLE PRECISION :: TEMP2
        DOUBLE PRECISION :: TEMP3
        DOUBLE PRECISION :: CORELATION(4),DOTPRO
        DOUBLE PRECISION,ALLOCATABLE  :: KORIGIN(:,:)
        DOUBLE PRECISION,ALLOCATABLE  :: KINVERSE(:,:) 
        DOUBLE PRECISION,ALLOCATABLE  :: UD(:),MAXVEC(:)
        DOUBLE PRECISION :: PI = 3.141592656
        DOUBLE PRECISION :: VAL
        DOUBLE PRECISION :: X(20),GRAD(20) 
        L = N / 4      
        ALLOCATE( KORIGIN(L,L) )        
        ALLOCATE( KINVERSE(L,L) )        
        ALLOCATE( UD(L) )        
        ALLOCATE( MAXVEC(L) )        
        
        IF (NEED_GRADIENT.NE.0) THEN
                GRAD(1:20)=0.0
        ENDIF

1003    FORMAT(F18.8,F18.8,F18.8,F18.8)
        !CALCULATE CORELATION
        CORELATION(1:4) =[5.0,4.0,2.5,6.0]/DSQRT(PI)! HEADING,PRESSURE,RADIUS,VFSPEED

        !CALCULATE U=TEMP1
        TEMP1 = 1.0
        DO I = 1,4
        TEMP1 = TEMP1 / DSQRT(1.0+4.0/CORELATION(I)**2)
        ENDDO
!        WRITE(*,*)  "GIVE ME U:",TEMP1

        !CALCULATE UD=TEMP2
        DO I = 0,L-1 
        TEMP2= 1.0 
        DO J = 1,4
        TEMP2 = TEMP2 * ( 1.0+2.0/CORELATION(J)**2 )**(-1.0/2.0) &
                      * EXP( -X(4*I+J)**2/(CORELATION(J)**2+2.0) ) 
        ENDDO 
        UD(I+1) = TEMP2 
        ENDDO 

1004    FORMAT(F11.5,F11.5,F11.5,F11.5,F11.5) !SHOULD BE CHANGED L=5
!        WRITE(*,*)  "GIVE ME UD:"
!        WRITE(*,1004) UD  

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

!        WRITE(*,*)  "GIVE ME KORIGIN:"
!        WRITE(*,1004) KORIGIN
!KORIGIN(:,1),KORIGIN(:,2),KORIGIN(:,3),&
!KORIGIN(:,4), KORIGIN(:,5)       
        
        !CALCULATE KINVERSE
        CALL INVERSE(KINVERSE,KORIGIN,L)
!        WRITE(*,*)  "GIVE ME KINVERSE:"
!        WRITE(*,1004) KINVERSE
!(:,1),KINVERSE(:,2),KINVERSE(:,3),&
!KINVERSE(:,4), KINVERSE(:,5)       

        !CALCULATE FVAL = TEMP1 - UD(T)*KINVERSE(-1)*UD
        MAXVEC =  MATMUL(KINVERSE,UD)
        DOTPRO =  DOT_PRODUCT(UD,MAXVEC)
!        WRITE(*,*) "GIVE ME MATMUL OF KINVERSE AND UD:"
!        WRITE(*,1004) MAXVEC         
!        WRITE(*,*) "GIVE ME DOT PRODUCT OF UD AND MATMUL:",DOTPRO     

        VAL = TEMP1 - DOTPRO
        F_DATA = F_DATA + 1

!        WRITE(*,*) "-----------------------------------------"
!        WRITE(*,*) "OPTIMAL POINT AT X(N) IS LISTED AS:"
!        WRITE(*,1003) X 
        WRITE(*,*) "FUNCTION VALUE AT OPTIMAL POINT IS:",VAL     
!        WRITE(*,*) "-----------------------------------------"

        END SUBROUTINE MYFUNC 


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

        !======= PREPARE L & U MATRICES
        !L IS A MATRIX OF ELIMINATION COEFFICIENT WITH THE DIAGONAL 1
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
     


        SUBROUTINE READRC(NR,NC)

        INTEGER(2),INTENT(INOUT)  :: NR,NC 
        INTEGER(2)    :: I,STAT
        CHARACTER(80) :: COLUMN

        !read the column of the sample 
        read(1,'(A)')  column
        NC = 1
        do I = 1,len(column)
               if (column(I:I)==',') then
                 NC = NC + 1
               endif
        enddo
        rewind(1)
 
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
     
        !check the row and column of the sample
1000    FORMAT(A,I3,A,I3) 
        write(*,1000) " THE NUMBER OF ROW IS",NR,&
                  "; THE NUMBER OF COLUMN IS",NC

      END SUBROUTINE  READRC


