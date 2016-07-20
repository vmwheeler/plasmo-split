#define BHFIELD_VERSION 'Oct 5, 2012'
#define cubareal real*8    

      PROGRAM full_spectrum

      IMPLICIT NONE
      REAL*8 UFTOL
      PARAMETER (UFTOL=1.0D-6)
      
      INTEGER I,NSTOPF,NSHELLS
      
      REAL*8 WLFAC(3)
      
      REAL*8 RADCOR,RADCOT,WAVEL,PI,CC,EPSVAC,MU,OMEGA
      REAL*8 REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      REAL*8 X,Y,Y1,Y2,Y3,Y4,YMAX,EXTMAX,RMAX
      REAL*8 QEXT,QSCA,QBACK,QABS,ANS
      REAL*8 INTENSITY,TEMP,IB
      
      COMPLEX*16 RFREL1,RFREL2
      
      REAL*8 XPMIN(3),XPMAX(3)
      REAL*8 XCMIN(3),XCMAX(3)
      
      CHARACTER*80 DIRNAM,FILNAM(3),FNAME(3),IFNAME,IFILNAM
      CHARACTER*50 FNLOGF
      
#ifdef CHECK_UNDERFLOW
      REAL*8 RDMIN
      CHARACTER*80 RDSTR
      common /RDMIN/RDMIN,RDSTR
      RDMIN=UFTOL
      RDSTR='Default'
#endif
      
#ifdef DATA_DIR
      DIRNAM=DATA_DIR
      PRINT *, DATA_DIR
#else
      DIRNAM='./'
#endif

C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Declare all relevant physical constants
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      WAVEL = 2.9D0
      RADCOR = 0.03D0
      RADCOT = 0.06D0
      PI=ACOS(-1.0D0)
      CC=2.99792458D8 ! light speed [m s-1]
      EPSVAC=1.0D7/(4.0D0*PI*CC*CC) ! eps0[F m-1]
      MU=4.0D0*PI*1.0D-7 ! assume non-magnetic (MU=MU0=const) [N A-2] 
      OMEGA=2.0D0*PI*CC/(WAVEL*1.0D-6) ! angular frequency [s-1]
      
      
C     log file
      FNLOGF='bhfield.log'
      OPEN(51,FILE=FNLOGF,STATUS='UNKNOWN')
      
C     The optical constants:
C     reference (background)
      FILNAM(1)='vac.nk'
      WLFAC(1)=1.0D-4 !factor so all wavelengths in micron
C     core
C      FILNAM(2)='test_material.nk'
C      FILNAM(2)='Ag_palik.nk'
      FILNAM(2)='Au_babar.nk'
      WLFAC(2)=1.0D0
C     shell
C      FILNAM(3)='test_material.nk'
C      FILNAM(3)='SiO2_palik.nk'
      FILNAM(3)='CeO2_patsalas.nk'
C       FILNAM(3)='vac.nk'
C      FILNAM(3)='Ag_palik.nk'
      WLFAC(3)=1.0D-3
      DO 901 I=1,3
       FNAME(I)=DIRNAM(1:MAX(INDEX(DIRNAM,' ')-1,1))//
     1          FILNAM(I)(1:MAX(INDEX(FILNAM(I),' ')-1,1))
  901 CONTINUE
      
      CALL OPTCON(1,WAVEL,FNAME,WLFAC, 
     1             REFMED,REFRE1,REFIM1,REFRE2,REFIM2)
      PRINT *, "Wavelength is: ", WAVEL
      PRINT *, "using this nk data for reference medium: ", FNAME(1)
      PRINT *, "got: ",REFMED," + i*",0.00," for ",WAVEL," micron"
      PRINT *, "using this nk data for core medium: ", FNAME(2)
      PRINT *, "got: ",REFRE1," + i*",REFIM1," for ",WAVEL," micron"
      PRINT *, "using this nk data for shell medium: ", FNAME(3)
      PRINT *, "got: ",REFRE2," + i*",REFIM2," for ",WAVEL," micron"
      
C     get relative refractive indices
      RFREL1=DCMPLX(REFRE1,REFIM1)/REFMED
      RFREL2=DCMPLX(REFRE2,REFIM2)/REFMED
 
C     Try re-calling OPTCON() to see if data is stored after loaded the first time
C      WAVEL = WAVEL + 1.0
C      CALL OPTCON(0,WAVEL,FNAME,WLFAC, 
C     1             REFMED,REFRE1,REFIM1,REFRE2,REFIM2)
C      PRINT *, ""
C      PRINT *, "AFTER CALLING FIRST TIME"
C      PRINT *, "Wavelength is: ", WAVEL
C      PRINT *, "using this nk data for reference medium: ", FNAME(1)
C      PRINT *, "got: ",REFMED," + i*",0.00," for ",WAVEL," micron"
C      PRINT *, "using this nk data for core medium: ", FNAME(2)
C      PRINT *, "got: ",REFRE1," + i*",REFIM1," for ",WAVEL," micron"
C      PRINT *, "using this nk data for shell medium: ", FNAME(3)
C      PRINT *, "got: ",REFRE2," + i*",REFIM2," for ",WAVEL," micron"      
 
 
C     define the particle size parameters for core and shell
      X=2.0D0*PI*RADCOR*REFMED/WAVEL
      Y=2.0D0*PI*RADCOT*REFMED/WAVEL
      
      
 
C
C max order of COEFF (A-W) necessary for field calc
C
      Y1=2.0D0*PI*RADCOR*ABS(DCMPLX(REFRE1,REFIM1))/WAVEL
      Y2=2.0D0*PI*RADCOT*ABS(DCMPLX(REFRE1,REFIM1))/WAVEL
      Y3=2.0D0*PI*RADCOT*ABS(DCMPLX(REFRE2,REFIM2))/WAVEL
      
C VMW fixed to only use radial coord... used to use cartesian and check each
C now uses worst case between radial and cart
      
      EXTMAX=MAX(ABS(XPMIN(1)),ABS(XPMAX(1)),
     1            ABS(XCMIN(1)),ABS(XCMIN(2)),ABS(XCMIN(3)),
     1            ABS(XCMAX(1)),ABS(XCMAX(2)),ABS(XCMAX(3)))


C     RMAX=EXTMAX*SQRT(2.0D0)
      RMAX=EXTMAX*SQRT(3.0D0)
      Y4=2.0D0*PI*RMAX*REFMED/WAVEL
      YMAX=MAX(Y1,Y2,Y3,Y4)
      NSTOPF=INT(YMAX+4.05D0*YMAX**0.3333D0+2.0D0)
C      NSTOPF = NSTOPF -4
C
   14 FORMAT ("CORE SIZE PARAM = ",F8.3,", COAT SIZE",
     1" PARAM = ",F8.3,", NSTOP(estim) = ",I3)
      WRITE(51,14) X,Y,NSTOPF
      
C
C calculate coefficients
C
      CALL BHCOAT(X,Y,RFREL1,RFREL2,NSTOPF,QEXT,QSCA,QBACK)
      QABS=QEXT-QSCA
      
      CALL IntegrateTestSemiInf(ANS)
      
      PRINT *, "Here it is: ", ANS
      
      IFILNAM="ASTMG173.csv"
      IFNAME=DIRNAM(1:MAX(INDEX(DIRNAM,' ')-1,1))//
     1          IFILNAM(1:MAX(INDEX(IFILNAM,' ')-1,1))
      CALL SOLARINTENSITY(1,WAVEL,IFNAME,INTENSITY,2)
      
      PRINT *, "And here is the intensity at ", WAVEL, "um"
      PRINT *, "  " , INTENSITY
      
      TEMP = 1000.0D0
      CALL BBINTENSITY(WAVEL,TEMP,IB)
      PRINT *, "AND now the blackbody intensity there at T=1000"
      PRINT *, "  ", IB
      
      UPPER = 
C     integrate over blackbody spectrum
      CALL INTEGRATEBB(WAVEL,TEMP,UPPER)
      
      STOP
      END

      
      SUBROUTINE IntegrateTestSemiInf(ANS)

      implicit none

      integer ndim, ncomp, nvec, last, seed, mineval, maxeval
      cubareal epsrel, epsabs
      
      parameter (ndim = 3)
      parameter (ncomp = 1)
      parameter (nvec = 1)
      parameter (epsrel = 1D-4)
      parameter (epsabs = 1D-12)
      parameter (last = 4)
      parameter (seed = 0)
      parameter (mineval = 0)
      parameter (maxeval = 50000)

      integer nstart, nincrease, nbatch, gridno
      integer*8 spin
      character*(*) statefile
      parameter (nstart = 1000)
      parameter (nincrease = 500)
      parameter (nbatch = 1000)
      parameter (gridno = 0)
      parameter (statefile = "")
      parameter (spin = -1)

      integer nnew, nmin
      cubareal flatness
      parameter (nnew = 1000)
      parameter (nmin = 2)
      parameter (flatness = 25D0)

      integer key1, key2, key3, maxpass
      cubareal border, maxchisq, mindeviation
      integer ngiven, ldxgiven, nextra
      parameter (key1 = 47)
      parameter (key2 = 1)
      parameter (key3 = 1)
      parameter (maxpass = 5)
      parameter (border = 0D0)
      parameter (maxchisq = 10D0)
      parameter (mindeviation = .25D0)
      parameter (ngiven = 0)
      parameter (ldxgiven = ndim)
      parameter (nextra = 0)

      integer key
      parameter (key = 0)
      
      REAL*8 ANS, userdata

      external sinexp

      cubareal integral(ncomp), error(ncomp), prob(ncomp)
      integer verbose, neval, fail, nregions
      character*16 env

      integer c

      
      call getenv("CUBAVERBOSE", env)
      verbose = 2
      read(env, *, iostat=fail, end=999, err=999) verbose
  999 continue

      print *, "----- integrating full particle volume ------"

      call Vegas(ndim, ncomp, sinexp, userdata, nvec, 
     &   epsrel, epsabs, verbose, seed, 
     &   mineval, maxeval, nstart, nincrease, nbatch, 
     &   gridno, statefile, spin, 
     &   neval, fail, integral, error, prob)

      print *, "neval    =", neval
      print *, "fail     =", fail
      print '(F20.12," +- ",F20.12,"   p = ",F8.3)', 
     &   (integral(c), error(c), prob(c), c = 1, ncomp)
      
      ANS = integral(1)

      RETURN
      END 
      
C************************************************************************

        integer function sinexp(ndim, xx, ncomp, ff)
        implicit none
        integer ndim, ncomp
        cubareal xx(*), ff(*),xp
        
#define x xx(1)
#define f ff(1)

C        f = SIN(x)*EXP(-x)
        xp = x/(1.0D0-x)
        f = SIN(xp)*EXP(-xp)/(1.0D0-x)/(1.0D0-x)

        sinexp = 0
        end
        
C ********************************************************************
C This thing should return the ASTM solar spectrum 
C ********************************************************************
C   WHICH = 1 --> Extraterrestrial radiation
C   WHICH = 2 --> Global tilt
C   WHICH = 3 --> Direct + circumsolar

      SUBROUTINE SOLARINTENSITY(INIT,WAVEL,FNAME,INTENSITY,WHICH)
C
C get n,k values at specified wavelength
C
      IMPLICIT NONE
      INTEGER NDMAX, WHICH
      PARAMETER (NDMAX=5000)
      INTEGER INIT,J
      REAL*8 WAVEL
      REAL*8 DATWL,ETR,GTILT,DIRCIRC,RAT,INTENSITY
      REAL*8 DOPT(1,4,NDMAX),WLFAC
      INTEGER ND
      CHARACTER*80 FNAME
      CHARACTER*255 BUF
      SAVE DOPT,ND
C
      IF(INIT.EQ.1) THEN
        OPEN(1,FILE=FNAME,STATUS='OLD')
        J=0
  562   CONTINUE
        READ(1,'(A)',END=561) BUF
        IF(BUF(1:1).EQ.'#'.OR.BUF(1:1).EQ.';') GO TO 562
        IF(INDEX(BUF,'.').EQ.0) GO TO 562
        READ(BUF,*,ERR=562) DATWL,ETR,GTILT,DIRCIRC
        IF(DATWL.EQ.0.0D0) GO TO 562
        J=J+1
        IF(J.GE.NDMAX) THEN
         WRITE(*,*) 'ndata > max for ',FNAME
         STOP
        END IF
        WLFAC = 1.0D-3 !hardcoded!
        DOPT(1,1,J)=DATWL*WLFAC
        DOPT(1,2,J)=ETR
        DOPT(1,3,J)=GTILT
        DOPT(1,4,J)=DIRCIRC
        GO TO 562
  561   CONTINUE
        CLOSE(1)
        ND=J
  570   FORMAT('Reading n,k(',I1,'): ndata = ',I4,' wavelength = ',
     *  E13.6,' - ',E13.6,' file = [',A,']')
        WRITE(*,570) 1,ND,DOPT(1,1,1),DOPT(1,1,ND),TRIM(FNAME)
        WRITE(51,570) 1,ND,DOPT(1,1,1),DOPT(1,1,ND),TRIM(FNAME)

      END IF
C
       IF(WAVEL.LT.DOPT(1,1,1).OR.WAVEL.GT.DOPT(1,1,ND)) THEN
        WRITE(*,*) 'wavelength out of range! ',WAVEL,' SI'
        STOP
       END IF
       DO 510 J=2,ND
C wavelength is in the increasing order
        IF(WAVEL.LE.DOPT(1,1,J)) GO TO 511
  510  CONTINUE
  511  CONTINUE
C simple interpolation
       RAT=(WAVEL-DOPT(1,1,J-1))/(DOPT(1,1,J)-DOPT(1,1,J-1))
       
       INTENSITY=DOPT(1,WHICH+1,J-1)+(DOPT(1,WHICH+1,J)
     &  -DOPT(1,WHICH+1,J-1))*RAT

C
C
      RETURN
      END
      
      
C ***********************************************************************
C Function to calculate blackbody emissive power at a given temperature
C  as a function of wavelength (in micron)
C ***********************************************************************
      
      SUBROUTINE BBINTENSITY(lambda, t, eb) 

      implicit none


      real*8 lambda, t, lambda_m
      real*8 eb
      REAL*8 c_first_radiation
      REAL*8 c_second_radiation
      ! -- Local declarations --

      real*8 beta
      
      lambda_m = lambda*1.0D-6
      
      PRINT *, "WAVEL in BB: ", lambda_m
      
      c_first_radiation = 3.74177118e-16
      c_second_radiation = 14387.75225e-06
      beta  = c_second_radiation / (lambda_m * t)
      eb = c_first_radiation  / (lambda_m ** 5 * (exp(beta) - 1.0D0))

      

      return
      end
      
      
C ***********************************************************************
C Function to calculate temperature derivative of the 
C  blackbody emissive power at a given temperature
C  as a function of wavelength (in micron)
C ***********************************************************************
      
      SUBROUTINE BBINTENSITYDT(lambda, t, ebdt) 

      implicit none


      real*8 lambda, t, lambda_m
      real*8 eb
      REAL*8 c_first_radiation
      REAL*8 c_second_radiation
      ! -- Local declarations --

      real*8 beta

      lambda_m = lambda*1.0D-6
      
      PRINT *, "WAVEL in BB: ", lambda_m
      
      c1 = 3.74177118e-16
      c2 = 14387.75225e-06

      dbdt = (c1*c2*E**(c2/(t*lambda_m)))
    &  /((-1.0D0 + E**(c2/(t*lambda_m)))**2*T**2*lambda_m**6)
    
      return
      end
