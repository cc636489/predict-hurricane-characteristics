        
        PROGRAM MAIN
        
        IMPLICIT NONE
        
        INTEGER                      :: I        !FOR LOOP
        INTEGER                      :: NR       !FOR READ FILE
        
        DOUBLE PRECISION,ALLOCATABLE :: HEADING(:)
        DOUBLE PRECISION,ALLOCATABLE :: PRESSURE(:)
        DOUBLE PRECISION,ALLOCATABLE :: RADIUS(:)
        DOUBLE PRECISION,ALLOCATABLE :: VFSPEED(:)
        DOUBLE PRECISION,ALLOCATABLE :: WEIGHTS(:)

        DOUBLE PRECISION             :: ALPHA,BETA
        DOUBLE PRECISION             :: MIU,K 
        DOUBLE PRECISION             :: RPMEAN,RPSIGMA,RPCOE
        DOUBLE PRECISION             :: VFMEAN,VFSIGMA
                

        !=====================================================
        !=========  IMPORT AND REARRANGE DATA ================
        !=====================================================
       
        !INPUT DATA
        WRITE(*,*) "INPUT DATA FOR ALPHA,BETA,MIU,K(USE COMMA TO&
                     SEPERATE THESE FOUR DATA):"
        READ(*,*) ALPHA,BETA,MIU,K   
        WRITE(*,*) "INPUT DATA FOR RPMEAN,RPSIGMA,RPCOE(USE COMMA TO&
                     SEPERATE THESE THREE DATA):"
        READ(*,*) RPMEAN,RPSIGMA,RPCOE   
        WRITE(*,*) "INPUT DATA FOR VFMEAN,VFSIGMA(USE COMMA TO&
                     SEPERATE THESE TWO DATA):"
        READ(*,*) VFMEAN,VFSIGMA   

        !READ ROWS TO NR
        CALL READRC(NR)
        
        !ALLOCATE STORM PARAMETERS
        ALLOCATE(HEADING(NR))
        ALLOCATE(PRESSURE(NR))
        ALLOCATE(RADIUS(NR))
        ALLOCATE(VFSPEED(NR))
        ALLOCATE(WEIGHTS(NR))
        
1000    FORMAT(F18.8,F18.8,F18.8,F18.8,F18.8)
1001    FORMAT(A,I3,A,F16.8,F16.8,F16.8,F16.8,F16.8)

        !READ DATA FOR EACH PARAMETERS
        OPEN(UNIT=1,FILE="DATA_WEIGHTS.CSV")
        DO I = 1,NR
           READ(1,1000) HEADING(I),PRESSURE(I),RADIUS(I),VFSPEED(I),&
                        WEIGHTS(I)
        ENDDO
        REWIND(1)
        CLOSE(1)

        !SHOW ORIGINAL DATA READING FROM INPUT FILE
        WRITE(*,*) "                         HEADING",&
                   "         PRESSURE",&
                   "          RADIUS",&
                   "        VFSPEED",&
                   "         WEIGHTS"

        DO I = 1,NR
        WRITE(*,1001) "SYNTHETIC STORM",I,":",HEADING(I),&
                       PRESSURE(I),RADIUS(I),VFSPEED(I),WEIGHTS(I)
        ENDDO        

        CALL ROSENBLATT(NR,ALPHA,BETA,MIU,K,&                       
                        RPMEAN,RPSIGMA,RPCOE,VFMEAN,VFSIGMA,&
                        HEADING,PRESSURE,RADIUS,VFSPEED)

        !SHOW PHYSICAL DATA READING AFTER TRANSFORMATION
        WRITE(*,*) "                         HEADING",&
                   "         PRESSURE",&
                   "          RADIUS",&
                   "        VFSPEED",&
                   "         WEIGHTS"

        DO I = 1,NR
        WRITE(*,1001) "SYNTHETIC STORM",I,":",HEADING(I),&
                       PRESSURE(I),RADIUS(I),VFSPEED(I),WEIGHTS(I)
        ENDDO        

        !READ DATA FOR EACH PARAMETERS
        OPEN(UNIT=2,FILE="DATA_ROSENBLATT.CSV")
        DO I = 1,NR
           WRITE(2,1000) HEADING(I),PRESSURE(I),RADIUS(I),VFSPEED(I),&
                        WEIGHTS(I)
        ENDDO
        CLOSE(2)
        
        DEALLOCATE(HEADING)
        DEALLOCATE(PRESSURE)
        DEALLOCATE(RADIUS)
        DEALLOCATE(VFSPEED)
        DEALLOCATE(WEIGHTS)
        
        END PROGRAM MAIN
        
        !================================================
        ! READ EXCEL.CSV FILE FOR GLOBAL OPTIMIZATION 
        !================================================
        SUBROUTINE READRC(NR)

        INTEGER,INTENT(OUT)  :: NR 
        INTEGER              :: I,STAT

        !read the row of the sample
        NR = 0
        OPEN(UNIT=1,FILE="DATA_WEIGHTS.CSV")
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
       

        !==============================================
        ! ROSENBLATT TRANSFORMATION AFTER ODERING TERMS
        !==============================================
        SUBROUTINE ROSENBLATT(L,ALPHA,BETA,MIU,K,&
                                RPMEAN,RPSIGMA,RPCOE,VFMEAN,VFSIGMA,&
                                HEADING,PRESSURE,RADIUS,VFSPEED)

        IMPLICIT NONE

        INTEGER                     :: PP,QQ        !DO ITERATION
        INTEGER,INTENT(IN)          :: L            !NUMBER OF STORMS TO BE GENERATED
        DOUBLE PRECISION            :: TEMP1,TEMP2  !DO ITERATION
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

1000    FORMAT(A,I3,A)
1001    FORMAT(A,F18.8) 
        

        OUTERLOOP: DO QQ = 1,L

        !SOLVING EQUATION FOR HEADING ANGLE 
        X0 = 0.5!THE RANGE IS (0.0,1.0) FOR STANDARD BETA DISTRIBUTION
        FVAL = 1.0/2.0*(1.0+ERF( HEADING(QQ)/SQRT(2.0) ))
        
        CALL BETACDF(X0,ALPHA,BETA,CDF)
        CALL BETAPDF(X0,ALPHA,BETA,PDF)
        ERROR = ABS(CDF-FVAL)
        
        PP = 1
        DO WHILE( ERROR>1.D-8  )
           
           X0=X0-(CDF-FVAL)/PDF

           CALL BETACDF(X0,ALPHA,BETA,CDF)
           CALL BETAPDF(X0,ALPHA,BETA,PDF)
           ERROR = ABS(CDF-FVAL)   
        
           write(*,*) "betaerror in this loop is :",ERROR   
           
           pp = pp + 1
        
        ENDDO

        HEADING(QQ) = X0*360.0 - 180.0
        
        write(*,1001) "HEADING(DEG) AT THE END IS:" ,HEADING(QQ)

        !SOLVING EQUATION FOR PRESSURE
        FVAL = 1.0/2.0*(1.0+ERF( PRESSURE(QQ)/SQRT(2.0)) )
        TEMP1  = EXP(-(48.0/MIU)**K)-EXP(-(115.0/MIU)**K)
        TEMP2  = EXP(-(115.0/MIU)**K)
        PRESSURE(QQ) = MIU*(-LOG( FVAL*TEMP1+TEMP2 ))**(1.0/K)
        write(*,1001) "PRESSURE(MB) AT THE END IS:" ,PRESSURE(QQ)

        !SOLVING EQUATION FOR RADIUS
        TEMP1         = RPMEAN*PRESSURE(QQ)**(RPCOE) 
        TEMP2         = RPSIGMA*PRESSURE(QQ)**(RPCOE)
        RPLOGMEAN  = LOG( TEMP1/SQRT(1.0+TEMP2/TEMP1**2.0)  )
        RPLOGSIGMA = SQRT( LOG(1.0+TEMP2/TEMP1**2.0)  )
        RADIUS(QQ) = EXP( RPLOGSIGMA*RADIUS(QQ)+RPLOGMEAN  )     
        write(*,1001) "RADIUS(NMI) AT THE END IS:" ,RADIUS(QQ)
 
        !SOLVING EQUATION FOR FORWARD SPEED
        VFLOGMEAN  = LOG( VFMEAN/SQRT(1.0+VFSIGMA/VFMEAN**2.0)  )
        VFLOGSIGMA = SQRT( LOG(1.0+VFSIGMA/VFMEAN**2.0)  )
        VFSPEED(QQ) = EXP( VFLOGSIGMA*VFSPEED(QQ)+VFLOGMEAN  )     
        write(*,1001) "VFSPEED(M/S) AT THE END IS:" ,VFSPEED(QQ)
        write(*,*) "=================================================="

        ENDDO OUTERLOOP

        END SUBROUTINE ROSENBLATT


        !=======================================
        !== CALCULATE CDF OF BETA DISTRIBUTION =
        !=======================================      
        SUBROUTINE BETACDF( X,ALPHA,BETA,CDF  )
           
        IMPLICIT NONE
        

        DOUBLE PRECISION,INTENT(INOUT)    :: CDF
        DOUBLE PRECISION,INTENT(IN)       :: X
        DOUBLE PRECISION,INTENT(IN)       :: ALPHA,BETA
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

        ALLOCATE  ( IWORK(LIMIT)  ) 
        ALLOCATE  ( WORK(LENW)  ) 

        A = 0.0
        
        CALL DQAG(F02,A,X,EPSABS,EPSREL,KEY,RESULTS,ABSERR,NEVAL,IER,&
                   LIMIT,LENW,LAST,IWORK,WORK)
        
        CDF = RESULTS * GAMMA(ALPHA+BETA)/GAMMA(ALPHA)/GAMMA(BETA)
        
        DEALLOCATE  ( IWORK ) 
        DEALLOCATE  ( WORK ) 
        

        CONTAINS        

           FUNCTION F02(X)
           ! THERE IS ONE DOUBT WHETHER I COULD USE ALPHA AND BETA AS
           ! PASSING ARGUMENT IN THE FUNCTION F02(X)
           IMPLICIT NONE
           DOUBLE PRECISION :: X
           DOUBLE PRECISION :: F02
   
           F02 =  X**(ALPHA-1.0) * (1.0-X)**(BETA-1.0)
           
           RETURN
   
           END FUNCTION F02

        END SUBROUTINE BETACDF
   
 
!        FUNCTION F02(X)
!        ! THERE IS ONE DOUBT WHETHER I COULD USE ALPHA AND BETA AS
!        ! PASSING ARGUMENT IN THE FUNCTION F02(X)
!        IMPLICIT NONE
!        DOUBLE PRECISION :: X
!        DOUBLE PRECISION :: F02
!
!        F02 =  X**(ALPHA-1.0) * (1.0-X)**(BETA-1.0)
!        
!        RETURN
!
!        END FUNCTION F02
        

        !=======================================
        !== CALCULATE PDF OF BETA DISTRIBUTION =
        !=======================================
        SUBROUTINE BETAPDF(X,ALPHA,BETA,PDF)

        IMPLICIT NONE
        DOUBLE PRECISION,INTENT(INOUT) :: PDF
        DOUBLE PRECISION,INTENT(IN)  :: X
        DOUBLE PRECISION,INTENT(IN)  :: ALPHA,BETA

        PDF = GAMMA(ALPHA+BETA)/GAMMA(ALPHA)/GAMMA(BETA)&
              * (1.0-X)**(BETA-1.0) * X**(ALPHA-1.0)
        
        END SUBROUTINE BETAPDF        
