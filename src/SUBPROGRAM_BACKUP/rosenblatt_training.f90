        PROGRAM ROSENBLATT

        IMPLICIT NONE

        INTEGER :: PP,QQ,RR !DO ITERATION
        DOUBLE PRECISION,PARAMETER :: PI = 3.141592656        

        DOUBLE PRECISION :: NODAL(20)
        DOUBLE PRECISION :: X0,CDF,PDF,ERROR
        DOUBLE PRECISION :: ALPHA,BETA     
        DOUBLE PRECISION :: MIU,K
        DOUBLE PRECISION :: RPMEAN,RPSIGMA,RPLOGMEAN,RPLOGSIGMA
        DOUBLE PRECISION :: VFMEAN,VFSIGMA,VFLOGMEAN,VFLOGSIGMA
        DOUBLE PRECISION :: FVAL !STANDARD NORMALCDF FOR 4 STORM PARA
        DOUBLE PRECISION :: HEADING(5),PRESSURE(5),RADIUS(5),VFSPEED(5)

        DOUBLE PRECISION,EXTERNAL :: BETACDF,BETAPDF
        DOUBLE PRECISION,EXTERNAL :: WEIBULLCDF,WEIBULLPDF
        DOUBLE PRECISION,EXTERNAL :: LOGNORMALCDF,LOGNORMALPDF


!IRES =   4 OPTIMIZATION STOPPED BECAUSE XTOL WAS REACHED!
! FOUND GLOBAL MINIMUM AT: 
!        0.01885636       -0.00000032       -1.30708725        0.00004492
!        0.56355330       -1.14705116        0.00000015       -0.00003945
!        0.56355491        1.14705054       -0.00000014       -0.00003957
!       -1.59271938        0.00000110       -0.00000001       -0.00004311
!        0.01885637        0.00000010        1.30708729        0.00013410
! MIN OBJECTIVE FUNCVAL =    7.5660829833268939E-003
! NLOPT HAS DONE     2377836 ITERATION TIMES FOR YOU!



!IRES =   4 OPTIMIZATION STOPPED BECAUSE XTOL WAS REACHED!
! FOUND GLOBAL MINIMUM AT: 
!       -0.55394364        1.14936579       -0.00000538        0.00000025
!       -0.01890471        0.00024167        1.30710344       -0.00000014
!       -0.57320007       -1.14467047       -0.00000560       -0.00000003
!       -0.01879674        0.00024117       -1.30707114        0.00000006
!        1.59261940       -0.01420361        0.00007895       -0.00000018
! MIN OBJECTIVE FUNCVAL =    7.5660840826334286E-003
! NLOPT HAS DONE     1039290 ITERATION TIMES FOR YOU!



!IRES =   4 OPTIMIZATION STOPPED BECAUSE XTOL WAS REACHED!
! FOUND GLOBAL MINIMUM AT: 
!       -1.24578235        0.78227758        0.00000870        0.00000025
!        1.20801172        0.81734448       -0.00000806       -0.00000053
!        0.00091388       -0.00995582        1.30727724        0.00000012
!        0.01573051       -1.20766178       -0.00000018        0.00000015
!        0.00089023       -0.00995543       -1.30727709       -0.00000003
! MIN OBJECTIVE FUNCVAL =    7.5680109458772793E-003
! NLOPT HAS DONE       77880 ITERATION TIMES FOR YOU!



!IRES =   4 OPTIMIZATION STOPPED BECAUSE XTOL WAS REACHED!
! FOUND GLOBAL MINIMUM AT: 
!       -0.56355388        1.14705111       -0.00001944        0.00000012
!       -0.01882937        0.00006925        1.30710881        0.00000001
!       -0.01888326       -0.00006923       -1.30706545       -0.00000008
!        1.59271926       -0.00000061       -0.00000933        0.00000010
!       -0.56355459       -1.14705081        0.00008722       -0.00000022
! MIN OBJECTIVE FUNCVAL =    7.5660829649334960E-003
! NLOPT HAS DONE     1167015 ITERATION TIMES FOR YOU!
       


!IRESIRES =   4 OPTIMIZATION STOPPED BECAUSE XTOL WAS REACHED!
! FOUND OBAL MINIMUM AT:
!        0.56355278        1.14705143       -0.00000011       -0.00009934
!       -1.59271924       -0.00000216        0.00000012        0.00037086
!        0.01885629       -0.00000016       -1.30708712       -0.00003155
!        0.01885646        0.00000017        1.30708723       -0.00003154
!        0.56355578       -1.14705036        0.00000012       -0.00010027
! MIN OBJECTIVE FUNCVAL =    7.5660830233109655E-003
! NLOPT HAS DONE     4072166 ITERATION TIMES FOR YOU! 
        
        
        !TYPING IN OPTIMIZER,THE NODAL VALUES,STORM PARAMETER(I AM USING
        !THE LAST SET OF NODAL VALUES from above)
        NODAL(1:20)=(/-1.59271924, 0.00037086, 0.00000012, 1.30708723,&
                      -1.14705036, 0.56355278,-0.00000216, 1.14705143,&
                      -1.30708712, 0.01885629,-0.00000016,-0.00003155,& 
                      -0.00009934, 0.01885646, 0.00000017,-0.00003154,& 
                      -0.00010027, 0.56355578, 0.00000012,-0.00000011/)

1000    FORMAT(A,I3,A)
1001    FORMAT(A,F18.8) 
1002    FORMAT(A,I3,A,F18.8,F18.8,F18.8,F18.8)   
        

        DO QQ = 0,4


        !=======================================
        !== SOLVING EQUATION FOR HEADING ANGLE =
        !=======================================
        
        X0 = 0.5!THE RANGE IS (0.0,1.0) FOR STANDARD BETA DISTRIBUTION
        ALPHA= 10.229
        BETA = 11.747
        FVAL = 1.0/2.0*(1.0+ERF( NODAL(1+4*QQ)/SQRT(2.0) ))
        
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

         
        HEADING(QQ+1) = X0*360.0 - 180.0
        
        write(*,*) "=================================================="
        write(*,1001) "HEADING(DEG) AT THE END IS:" ,HEADING(QQ+1)
        write(*,*) "=================================================="

        

        !=======================================
        !=== SOLVING EQUATION FOR PRESSURE =====
        !=======================================
        X0  = 78.0!THE RANGE IS [48,111]
        MIU = 48.6
        K   = 1.8               
        FVAL = 1.0/2.0*(1.0+ERF( NODAL(2+4*QQ)/SQRT(2.0)) )
 
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
        
        PRESSURE(QQ+1) = X0
        
        write(*,*) "=================================================="
        write(*,1001) "PRESSURE(MB) AT THE END IS:" ,PRESSURE(QQ+1)
        write(*,*) "=================================================="


        !=======================================
        !=== SOLVING EQUATION FOR RADIUS =======
        !=======================================
        
!        LOGSIGMA = SQRT( LOG(1.0+SIGMA**2.0/MEAN**2.0)  ) TORO'S VERSION
!        LOGMEAN = LOG(MEAN) * EXP(-LOGSIGMA**2.0/2.0)

!        LOGMEAN =  LOG( MEAN/SQRT(1.0+SIGMA/MEAN**2.0)  ) WIKI VERSION
!        LOGSIGMA = SQRT( LOG(1.0+SIGMA/MEAN**2.0)  )

        X0 = 20.0 !ARBITRARY POSITIVE VALUE, WILL NOT CHANGE RESULT  
        RPMEAN = 406.2*PRESSURE(QQ+1)**(-0.711)
        RPSIGMA= 187.7*PRESSURE(QQ+1)**(-0.711)
        RPLOGMEAN  = LOG( RPMEAN/SQRT(1.0+RPSIGMA/RPMEAN**2.0)  )
        RPLOGSIGMA = SQRT( LOG(1.0+RPSIGMA/RPMEAN**2.0)  )
        FVAL = 1.0/2.0*(1.0+ERF(NODAL(3+4*QQ)/SQRT(2.0) ))

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
        
        RADIUS(QQ+1) = X0

        write(*,*) "=================================================="
        write(*,1001) "RADIUS(NMI) AT THE END IS:" ,RADIUS(QQ+1)
        write(*,*) "=================================================="
      

 
        !=======================================
        !== SOLVING EQUATION FOR FORWARD SPEED =
        !=======================================
       
        X0  = 10.0
        VFMEAN = 6.6 
        VFSIGMA= 2.8
        VFLOGMEAN  = LOG( VFMEAN/SQRT(1.0+VFSIGMA/VFMEAN**2.0)  )
        VFLOGSIGMA = SQRT( LOG(1.0+VFSIGMA/VFMEAN**2.0)  )
        FVAL = 1.0/2.0*(1.0+ERF(NODAL(4+4*QQ)/SQRT(2.0) )) 

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
        
        VFSPEED(QQ+1) = X0

        write(*,*) "=================================================="
        write(*,1001) "VFSPEED(M/S) AT THE END IS:" ,VFSPEED(QQ+1)
        write(*,*) "=================================================="

        ENDDO


        WRITE(*,*) "                         HEADING",&
                   "           PRESSURE",&
                   "           RADIUS",&
                   "            VFSPEED"

        DO RR = 1,5
        WRITE(*,1002) "SYNTHETIC STORM",RR,":",HEADING(RR),&
                       PRESSURE(RR),RADIUS(RR),VFSPEED(RR)
        ENDDO        

        END PROGRAM ROSENBLATT


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
           
           CALL DQAG(F02,A,X,EPSABS,EPSREL,KEY,RESULTS,ABSERR,NEVAL,IER,&
                      LIMIT,LENW,LAST,IWORK,WORK)

           !THIS IS TRAPEDOIZAL METHODS FOR CALCULATING INTEGRATION
!           RESULTS = 0.5*X*X**(ALPHA-1)*(1-X)**(BETA-1)

           BETACDF = RESULTS * GAMMA(ALPHA+BETA)/GAMMA(ALPHA)/GAMMA(BETA)
           
           RETURN           

        END FUNCTION BETACDF

        
        FUNCTION F02(X)

          ! THERE IS ONE DOUBT WHETHER I COULD USE ALPHA AND BETA AS
          ! PASSING ARGUMENT IN THE FUNCTION F02(X)
          IMPLICIT NONE
          DOUBLE PRECISION :: X
          DOUBLE PRECISION :: F02
          DOUBLE PRECISION :: ALPHA = 10.229
          DOUBLE PRECISION :: BETA  = 11.747

          F02 =  X**(ALPHA-1.0) * (1.0-X)**(BETA-1.0)
          
          RETURN

        END FUNCTION F02
        

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
       

