      PROGRAM  main
        
        IMPLICIT NONE       
        !=====================================================
        !=================  DATA DECLARATION  ================
        !=====================================================
        
        INTEGER(4)                   :: N,L               !FOR #OPTIMPARA & #STORMS 
        INTEGER(4)                   :: PP                !FOR DIFFERNET CATEGORY LOOP
        INTEGER(4)                   :: RR                !FOR PLOTING STORM PARAMTER(END)

        INTEGER(4)                   :: F_DATA            !FOR OPTIM : LIBRARY FIXED
        INTEGER(8)                   :: IRES              !FOR OPTIM : LIBRARY FIXED
        INTEGER*8                    :: OPT               !FOR OPTIM : LIBRARY FIXED
        INTEGER*8                    :: LOCAL_OPT         !FOR OPTIM : LIBRARY FIXED
        DOUBLE PRECISION             :: MINF              !FOR OPTIM : LIBRARY FIXED

        DOUBLE PRECISION,PARAMETER   :: PI = 3.1415926    !FOR OPTIM AND WEIGHTS
        DOUBLE PRECISION             :: CORELATION(4)     !FOR OPTIM AND WEIGHTS
        DOUBLE PRECISION             :: TEMP2,TEMP3       !FOR OPTIM AND WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: X(:)              !FOR OPTIM AND WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: UD(:)             !FOR OPTIM AND WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: WEIGHT(:)         !FOR OPTIM AND WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: KORIGIN(:,:)      !FOR OPTIM AND WEIGHTS
        DOUBLE PRECISION,ALLOCATABLE :: KINVERSE(:,:)     !FOR OPTIM AND WEIGHTS
     
        DOUBLE PRECISION             :: ALPHA,BETA              !FOR ROSENBLATT
        DOUBLE PRECISION             :: MIU,K                   !FOR ROSENBLATT
        DOUBLE PRECISION             :: RPMEAN,RPSIGMA,RPCOE    !FOR ROSENBLATT
        DOUBLE PRECISION             :: VFMEAN,VFSIGMA          !FOR ROSENBLATT

        DOUBLE PRECISION,EXTERNAL    :: BETACDF,BETAPDF
        DOUBLE PRECISION,EXTERNAL    :: WEIBULLCDF,WEIBULLPDF
        DOUBLE PRECISION,EXTERNAL    :: LOGNORMALCDF,LOGNORMALPDF
        
        DOUBLE PRECISION,ALLOCATABLE :: HEADING(:)
        DOUBLE PRECISION,ALLOCATABLE :: PRESSURE(:)
        DOUBLE PRECISION,ALLOCATABLE :: RADIUS(:)
        DOUBLE PRECISION,ALLOCATABLE :: VFSPEED(:)

        EXTERNAL MYFUNC
        INCLUDE 'nlopt.f'        

        ALPHA= 10.229
        BETA = 11.747
        MIU = 48.6
        K   = 1.8               
        RPMEAN = 406.2
        RPSIGMA= 187.7
        RPCOE = -0.711
        VFMEAN = 6.6 
        VFSIGMA= 2.8

        !=====================================================
        !======== NODAL VALUE OPTIMIZATION (CASE PP=2) =======
        !=====================================================
        
        PP = 1  
        
        ! SETTING OPTIMPARAS N FOR THE USE OF NLOPT LIBRARY
        L = 13    ! #GOALNODES,13,5,7,7    
        N = 4*L   ! #OPTINPARA,52,20,28,28

        ALLOCATE(X(N))
        ALLOCATE(UD(L))
        ALLOCATE(WEIGHT(L))        
        ALLOCATE(KORIGIN(L,L))        
        ALLOCATE(KINVERSE(L,L))        
        ALLOCATE(HEADING(L))
        ALLOCATE(PRESSURE(L))
        ALLOCATE(RADIUS(L))
        ALLOCATE(VFSPEED(L))


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
        IRES = 1.0 ! THIS IS SUPER IMPORTANT,otherwise,directly fail!!
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
        
1002    FORMAT(A,I3,A,F18.8,F18.8,F18.8,F18.8)   
     
        !===============================================
        !== PRESENTING NODAL VALUES IN NORMAL  SPACE ===
        !===============================================
        CALL PRE_ROSENBLATT(L,N,X,HEADING,PRESSURE,RADIUS,VFSPEED)

        WRITE(*,*) "                         HEADING",&
                   "           PRESSURE",&
                   "           RADIUS",&
                   "            VFSPEED"

        DO RR = 1,L
        WRITE(*,1002) "SYNTHETIC STORM",RR,":",HEADING(RR),&
                       PRESSURE(RR),RADIUS(RR),VFSPEED(RR)
        ENDDO        
        
        !===============================================
        !== PRESENTING NODAL VALUES IN PHYSICAL SPACE  =
        !===============================================
        CALL ROSENBLATT(L,ALPHA,BETA,MIU,K,&                       
                        RPMEAN,RPSIGMA,RPCOE,VFMEAN,VFSIGMA,&
                        HEADING,PRESSURE,RADIUS,VFSPEED)

        WRITE(*,*) "                         HEADING",&
                   "           PRESSURE",&
                   "           RADIUS",&
                   "            VFSPEED"

        DO RR = 1,L
        WRITE(*,1002) "SYNTHETIC STORM",RR,":",HEADING(RR),&
                       PRESSURE(RR),RADIUS(RR),VFSPEED(RR)
        ENDDO        



        END PROGRAM main
        
        ! --------------------------------------------------------------------
        ! INTEGER FUNCTION  FindMinimum():
        !    This function returns the location of the minimum in the section
        ! between Start and End.
        ! --------------------------------------------------------------------
        
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
        
        ! --------------------------------------------------------------------
        ! SUBROUTINE  Swap():
        !    This subroutine swaps the values of its two formal arguments.
        ! --------------------------------------------------------------------
        
           SUBROUTINE  Swap(a, b)

              IMPLICIT  NONE

              DOUBLE PRECISION, INTENT(INOUT) :: a, b
              DOUBLE PRECISION                :: Temp
        
              Temp = a
              a    = b
              b    = Temp

           END SUBROUTINE  Swap
        
        ! --------------------------------------------------------------------
        ! SUBROUTINE  Sort():
        !    This subroutine receives an array x() and sorts it into ascending
        ! order.
        ! --------------------------------------------------------------------
        
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

        !=======================================
        ! REARRANGE TERMS IN X(4L) 
        !=======================================
        SUBROUTINE PRE_ROSENBLATT(L,N,X,HEADING,PRESSURE,RADIUS,VFSPEED)

        IMPLICIT NONE
        
        INTEGER                        :: TT 
        INTEGER,INTENT(IN)             :: L
        INTEGER,INTENT(IN)             :: N
        DOUBLE PRECISION,INTENT(INOUT) :: X(N)
        DOUBLE PRECISION,INTENT(INOUT) :: HEADING(L)
        DOUBLE PRECISION,INTENT(INOUT) :: PRESSURE(L)
        DOUBLE PRECISION,INTENT(INOUT) :: RADIUS(L)
        DOUBLE PRECISION,INTENT(INOUT) :: VFSPEED(L)

        !BEGIN ARRANGE THE NODAL VALUES IN NORMAL SPACE GIVEN ME X(4L)
        !TO REORDER X(N) INTO HEADING,PRESSURE,RADIUS,VFSPEED
        CALL Sort(N,X) 
        
        DO TT = 1,L
        
        HEADING(TT)  = X(TT)
        PRESSURE(TT) = X(TT+L)        
        RADIUS(TT)   = X(TT+L*2) 
        VFSPEED(TT)  = X(TT+L*3) 
        
        ENDDO

        !END ARRANGE THE NODAL VALUES IN NOMAL SPACE RETURN YOU NODAL(4L)
        
        END SUBROUTINE PRE_ROSENBLATT


        !==============================================
        ! ROSENBLATT TRANSFORMATION AFTER ODERING TERMS
        !==============================================
        SUBROUTINE ROSENBLATT(L,ALPHA,BETA,MIU,K,&
                                RPMEAN,RPSIGMA,RPCOE,VFMEAN,VFSIGMA,&
                                HEADING,PRESSURE,RADIUS,VFSPEED)

        IMPLICIT NONE

        INTEGER :: PP,QQ,MM,NN !DO ITERATION
        INTEGER,INTENT(IN) :: L !NUMBER OF STORMS TO BE GENERATED
        DOUBLE PRECISION,PARAMETER  :: PI = 3.141592656        

        DOUBLE PRECISION,INTENT(IN) :: ALPHA,BETA     
        DOUBLE PRECISION,INTENT(IN) :: MIU,K
        DOUBLE PRECISION,INTENT(IN) :: RPMEAN,RPSIGMA,RPCOE
        DOUBLE PRECISION,INTENT(IN) :: VFMEAN,VFSIGMA

        DOUBLE PRECISION :: X0,CDF,PDF,ERROR
        DOUBLE PRECISION :: RPLOGMEAN,RPLOGSIGMA
        DOUBLE PRECISION :: VFLOGMEAN,VFLOGSIGMA
        DOUBLE PRECISION :: FVAL !STANDARD NORMALCDF FOR 4 STORM PARA

        DOUBLE PRECISION,INTENT(INOUT) :: HEADING(L)
        DOUBLE PRECISION,INTENT(INOUT) :: PRESSURE(L)
        DOUBLE PRECISION,INTENT(INOUT) :: RADIUS(L)
        DOUBLE PRECISION,INTENT(INOUT) :: VFSPEED(L)

        DOUBLE PRECISION :: BETACDF,BETAPDF
        DOUBLE PRECISION :: WEIBULLCDF,WEIBULLPDF
        DOUBLE PRECISION :: LOGNORMALCDF,LOGNORMALPDF
        
1000    FORMAT(A,I3,A)
1001    FORMAT(A,F18.8) 
        

        OUTERLOOP: DO QQ = 1,L
        !=======================================
        !== SOLVING EQUATION FOR HEADING ANGLE =
        !=======================================
        X0 = 0.5!THE RANGE IS (0.0,1.0) FOR STANDARD BETA DISTRIBUTION
        FVAL = 1.0/2.0*(1.0+ERF( HEADING(QQ)/SQRT(2.0) ))
        
        write(*,*) "=========================================="
        write(*,1001) "fvalheading is:",fval

        CDF = BETACDF(X0,ALPHA,BETA)
        PDF = BETAPDF(X0,ALPHA,BETA)
        ERROR = ABS(CDF-FVAL)
        
        write(*,1001) "initial betacdf is :",cdf 
        write(*,1001) "initial betapdf is:",pdf 

        PP = 1
        DO WHILE( ERROR>1.D-8  )
           
           X0=X0-(CDF-FVAL)/PDF

           CDF = BETACDF(X0,ALPHA,BETA)
           PDF = BETAPDF(X0,ALPHA,BETA)
           ERROR = ABS(CDF-FVAL)   
        
           write(*,1000) "-------------- loop",pp,"-------------------"
           write(*,1001) "x0 in this loop is:",X0
           write(*,1001) "betacdf in this loop is :",cdf 
           write(*,1001) "betapdf in this loop is :",pdf 
           write(*,*) "betaerror in this loop is :",ERROR   
           
           pp = pp + 1
        
        ENDDO

         
        HEADING(QQ) = X0*360.0 - 180.0
        
        write(*,*) "=================================================="
        write(*,1001) "HEADING(DEG) AT THE END IS:" ,HEADING(QQ)
        write(*,*) "=================================================="

        !=======================================
        !=== SOLVING EQUATION FOR PRESSURE =====
        !=======================================
        X0  = 78.0!THE RANGE IS [48,111]
        FVAL = 1.0/2.0*(1.0+ERF( PRESSURE(QQ)/SQRT(2.0)) )
 
        write(*,1001) "fvalpressure is:",fval

        CDF = WEIBULLCDF(X0,MIU,K)
        PDF = WEIBULLPDF(X0,MIU,K)
        ERROR = ABS(CDF-FVAL)

        write(*,*) "initial weibullcdf is :",cdf 
        write(*,*) "initial weibullpdf is :",pdf 

        PP = 1
        DO WHILE( ERROR>1.D-8  )
           
           X0=X0-(CDF-FVAL)/PDF

           CDF = WEIBULLCDF(X0,MIU,K)
           PDF = WEIBULLPDF(X0,MIU,K)           
           ERROR = ABS(CDF-FVAL)

           write(*,1000) "-------------- loop",pp,"-------------------"
           write(*,1001) "x0 in this loop is:",X0
           write(*,1001) "weibullcdf in this loop is :",cdf 
           write(*,1001) "weibullpdf in this loop is :",pdf 
           write(*,*) "weibull error in this loop is :",ERROR   
           
           pp = pp + 1
        
        ENDDO
        
        PRESSURE(QQ) = X0
        
        write(*,*) "=================================================="
        write(*,1001) "PRESSURE(MB) AT THE END IS:" ,PRESSURE(QQ)
        write(*,*) "=================================================="

        !=======================================
        !=== SOLVING EQUATION FOR RADIUS =======
        !=======================================
        X0 = 20.0 !ARBITRARY POSITIVE VALUE, WILL NOT CHANGE RESULT  
        MM = RPMEAN*PRESSURE(QQ)**(RPCOE) !-0.711 NEED TO PARAMETRIZE
        NN = RPSIGMA*PRESSURE(QQ)**(RPCOE) !-0.711 NEED TO PARAMETRIZE
        RPLOGMEAN  = LOG( MM/SQRT(1.0+NN/MM**2.0)  )
        RPLOGSIGMA = SQRT( LOG(1.0+NN/MM**2.0)  )
        FVAL = 1.0/2.0*(1.0+ERF( RADIUS(QQ)/SQRT(2.0) ))

        write(*,1001) "fvalradius is:",fval
        
        CDF = LOGNORMALCDF(X0,RPLOGMEAN,RPLOGSIGMA) !CDF IN PHYSICAL SPACE
        PDF = LOGNORMALPDF(X0,RPLOGMEAN,RPLOGSIGMA)
        ERROR = ABS(CDF-FVAL)        

        write(*,1001) "initial lognormalcdf for radius is :",cdf
        write(*,1001) "initial lognormalpdf for radius is :",pdf
        
        !START SOLVE EQUATION: FVAL=CDF using newton method 
        PP = 1
        DO WHILE( ERROR>1.D-8  )
           
           X0=X0-(CDF-FVAL)/PDF
           
           CDF = LOGNORMALCDF(X0,RPLOGMEAN,RPLOGSIGMA)
           PDF = LOGNORMALPDF(X0,RPLOGMEAN,RPLOGSIGMA)
           ERROR = ABS(CDF-FVAL)        
        
           write(*,1000) "-------------- loop",pp,"-------------------"
           write(*,1001) "x0 in this loop is:",X0
           write(*,1001) "lognormalcdf for radius in this loop is :",cdf
           write(*,1001) "lognormalpdf for radius in this loop is :",pdf
           write(*,*) "radius error in this loop is :",ERROR   
        
           pp = pp + 1

        ENDDO
        
        RADIUS(QQ) = X0

        write(*,*) "=================================================="
        write(*,1001) "RADIUS(NMI) AT THE END IS:" ,RADIUS(QQ)
        write(*,*) "=================================================="
 
        !=======================================
        !== SOLVING EQUATION FOR FORWARD SPEED =
        !=======================================
        X0  = 10.0
        VFLOGMEAN  = LOG( VFMEAN/SQRT(1.0+VFSIGMA/VFMEAN**2.0)  )
        VFLOGSIGMA = SQRT( LOG(1.0+VFSIGMA/VFMEAN**2.0)  )
        FVAL = 1.0/2.0*(1.0+ERF( VFSPEED(QQ)/SQRT(2.0) )) 

        write(*,1001) "fvalvfspeed is:",fval
        
        CDF = LOGNORMALCDF(X0,VFLOGMEAN,VFLOGSIGMA)
        PDF = LOGNORMALPDF(X0,VFLOGMEAN,VFLOGSIGMA)
        ERROR = ABS(CDF-FVAL)        

        write(*,1001) "initial lognormalcdf for vf is :",cdf
        write(*,1001) "initial lognormalpdf for vf is :",pdf
        
        PP = 1
        DO WHILE( ERROR>1.D-8  )

           X0=X0-(CDF-FVAL)/PDF
 
           CDF = LOGNORMALCDF(X0,VFLOGMEAN,VFLOGSIGMA)
           PDF = LOGNORMALPDF(X0,VFLOGMEAN,VFLOGSIGMA)
           ERROR = ABS(CDF-FVAL)        
        
           write(*,1000) "-------------- loop",pp,"-------------------"
           write(*,1001) "x0 in this loop is:",X0
           write(*,1001) "lognormalcdf for vf in this loop is :",cdf 
           write(*,1001) "lognormalpdf for vf in this loop is :",pdf
           write(*,*) "vf error in this loop is :",ERROR   

           pp = pp + 1

        ENDDO
        
        VFSPEED(QQ) = X0

        write(*,*) "=================================================="
        write(*,1001) "VFSPEED(M/S) AT THE END IS:" ,VFSPEED(QQ)
        write(*,*) "=================================================="

        ENDDO OUTERLOOP


        END SUBROUTINE ROSENBLATT


      !=======================================
      !== CALCULATE CDF OF BETA DISTRIBUTION =
      !=======================================      
      FUNCTION BETACDF( X,ALPHA,BETA  )
         
      IMPLICIT NONE
      

      DOUBLE PRECISION                  :: BETACDF
      DOUBLE PRECISION                  :: X
      DOUBLE PRECISION                  :: ALPHA,BETA
      DOUBLE PRECISION                  :: A
      DOUBLE PRECISION                  :: ABSERR,RESULTS
      DOUBLE PRECISION                  :: EPSABS = 0.0E+00     
      DOUBLE PRECISION                  :: EPSREL = 0.001E+00   
      INTEGER(4)                        :: LIMIT = 1000
      INTEGER(4)                        :: LENW  = 4000
      INTEGER(4)                        :: IER
      INTEGER(4)                        :: LAST,NEVAL       
      INTEGER(4),ALLOCATABLE            :: IWORK(:),WORK(:)
      INTEGER(4)                        :: KEY = 6           

      DOUBLE PRECISION,EXTERNAL :: F02
      
      ALLOCATE  ( IWORK(LIMIT)  ) 
      ALLOCATE  ( WORK(LENW)  ) 

      A = 0.0
      
!      CALL DQAG(F02,A,X,EPSABS,EPSREL,KEY,RESULTS,ABSERR,NEVAL,IER,&
!                 LIMIT,LENW,LAST,IWORK,WORK)

      !THIS IS TRAPEDOIZAL METHODS FOR CALCULATING INTEGRATION
      RESULTS = 0.5*X*X**(ALPHA-1)*(1-X)**(BETA-1)

      BETACDF = RESULTS * GAMMA(ALPHA+BETA)/GAMMA(ALPHA)/GAMMA(BETA)
      
      RETURN           

      END FUNCTION BETACDF
   
 
!      FUNCTION F02(X)
!      ! THERE IS ONE DOUBT WHETHER I COULD USE ALPHA AND BETA AS
!      ! PASSING ARGUMENT IN THE FUNCTION F02(X)
!      IMPLICIT NONE
!      DOUBLE PRECISION :: X
!      DOUBLE PRECISION :: F02
!      DOUBLE PRECISION :: ALPHA = 10.774
!      DOUBLE PRECISION :: BETA = 11.565
!
!      F02 =  X**(ALPHA-1.0) * (1.0-X)**(BETA-1.0)
!      
!      RETURN
!
!      END FUNCTION F02
      

      !=======================================
      !== CALCULATE PDF OF BETA DISTRIBUTION =
      !=======================================
      FUNCTION BETAPDF(X,ALPHA,BETA)

      IMPLICIT NONE
      DOUBLE PRECISION :: BETAPDF
      DOUBLE PRECISION :: X
      DOUBLE PRECISION :: ALPHA,BETA

      BETAPDF = GAMMA(ALPHA+BETA)/GAMMA(ALPHA)/GAMMA(BETA)&
                   * (1.0-X)**(BETA-1.0) * X**(ALPHA-1.0)
      
      RETURN          

      END FUNCTION BETAPDF        


      !=======================================
      ! CALCULATE CDF OF WEIBULL DISTRIBUTION 
      !=======================================
      FUNCTION WEIBULLCDF(X,MIU,K)

      IMPLICIT NONE
      DOUBLE PRECISION :: WEIBULLCDF
      DOUBLE PRECISION :: X
      DOUBLE PRECISION :: MIU,K
      
      WEIBULLCDF = 1.0 - (  EXP(-(X/MIU)**K)-EXP(-(70.0/MIU)**K) )/&
                         (  EXP(-(48.0/MIU)**K)-EXP(-(70.0/MIU)**K) )  


      RETURN
      
      END FUNCTION WEIBULLCDF


      !=======================================
      ! CALCULATE PDF OF WEIBULL DISTRIBUTION 
      !=======================================
      FUNCTION WEIBULLPDF(X,MIU,K)

      IMPLICIT NONE
      DOUBLE PRECISION :: WEIBULLPDF
      DOUBLE PRECISION :: X
      DOUBLE PRECISION :: MIU,K


      WEIBULLPDF = (K/X)*(X/MIU)**K * EXP(-(X/MIU)**K)/&
                   (  EXP(-(48.0/MIU)**K)-EXP(-(70.0/MIU)**K) )
      RETURN
      
      END FUNCTION WEIBULLPDF

      
      !=======================================
      !CALCULATE CDF OF LOGNORMAL DISTRIBUTION 
      !=======================================
      FUNCTION LOGNORMALCDF(X,LOGMEAN,LOGSIGMA)
      
      IMPLICIT NONE
      DOUBLE PRECISION :: LOGNORMALCDF
      DOUBLE PRECISION :: X
      DOUBLE PRECISION :: LOGMEAN,LOGSIGMA

      LOGNORMALCDF=1.0/2.0+1.0/2.0*ERF((LOG(X)-LOGMEAN)/&
                           SQRT(2.0)/LOGSIGMA)

      RETURN

      END FUNCTION LOGNORMALCDF

      
      !=======================================
      !CALCULATE PDF OF LOGNORMAL DISTRIBUTION 
      !=======================================
      FUNCTION LOGNORMALPDF(X,LOGMEAN,LOGSIGMA)
      
      IMPLICIT NONE
      DOUBLE PRECISION :: LOGNORMALPDF
      DOUBLE PRECISION :: X
      DOUBLE PRECISION :: LOGMEAN,LOGSIGMA
      DOUBLE PRECISION,PARAMETER :: PI = 3.141592656
      
      LOGNORMALPDF = 1.0/X/LOGSIGMA/SQRT(2.0*PI)&
                      *EXP(-(LOG(X)-LOGMEAN)**2.0/2.0/LOGSIGMA**2.0)

      RETURN

      END FUNCTION LOGNORMALPDF
      
 
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
      ! COST FUNCTION FOR GLOBAL OPTIMIZATION 
      !================================================
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

1003  FORMAT(F18.8,F18.8,F18.8,F18.8)
      !CALCULATE CORELATION
      CORELATION(1:4) =[5.0,4.0,2.5,6.0]/DSQRT(PI)! HEADING,PRESSURE,RADIUS,VFSPEED

      !CALCULATE U=TEMP1
      TEMP1 = 1.0
      DO I = 1,4
      TEMP1 = TEMP1 / DSQRT(1.0+4.0/CORELATION(I)**2)
      ENDDO
!      WRITE(*,*)  "GIVE ME U:",TEMP1

      !CALCULATE UD=TEMP2
      DO I = 0,L-1 
      TEMP2= 1.0 
      DO J = 1,4
      TEMP2 = TEMP2 * ( 1.0+2.0/CORELATION(J)**2 )**(-1.0/2.0) &
                    * EXP( -X(4*I+J)**2/(CORELATION(J)**2+2.0) ) 
      ENDDO 
      UD(I+1) = TEMP2 
      ENDDO 

1004  FORMAT(F11.5,F11.5,F11.5,F11.5,F11.5) !SHOULD BE CHANGED L=5
!      WRITE(*,*)  "GIVE ME UD:"
!      WRITE(*,1004) UD  

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

!      WRITE(*,*)  "GIVE ME KORIGIN:"
!      WRITE(*,1004) KORIGIN
!KORIGIN(:,1),KORIGIN(:,2),KORIGIN(:,3),&
!KORIGIN(:,4), KORIGIN(:,5)       
        
      !CALCULATE KINVERSE
      CALL INVERSE(KINVERSE,KORIGIN,L)
!      WRITE(*,*)  "GIVE ME KINVERSE:"
!      WRITE(*,1004) KINVERSE
!(:,1),KINVERSE(:,2),KINVERSE(:,3),&
!KINVERSE(:,4), KINVERSE(:,5)       

      !CALCULATE FVAL = TEMP1 - UD(T)*KINVERSE(-1)*UD
      MAXVEC =  MATMUL(KINVERSE,UD)
      DOTPRO =  DOT_PRODUCT(UD,MAXVEC)
!      WRITE(*,*) "GIVE ME MATMUL OF KINVERSE AND UD:"
!      WRITE(*,1004) MAXVEC         
!      WRITE(*,*) "GIVE ME DOT PRODUCT OF UD AND MATMUL:",DOTPRO     

      VAL = TEMP1 - DOTPRO
      F_DATA = F_DATA + 1

!      WRITE(*,*) "-----------------------------------------"
!      WRITE(*,*) "OPTIMAL POINT AT X(N) IS LISTED AS:"
!      WRITE(*,1003) X 
      WRITE(*,*) "FUNCTION VALUE AT OPTIMAL POINT IS:",VAL     
!      WRITE(*,*) "-----------------------------------------"

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
      ! READ EXCEL.CSV FILE FOR GLOBAL OPTIMIZATION 
      !================================================
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


