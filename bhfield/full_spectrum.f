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
      
      integer npar
      real*8 params(1)
      external diffIsolPabs, bbfunc, bbint
      
      
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
      WAVEL = 0.28D0
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
      
      PRINT *,"initializing optical data"
      CALL OPTCON(1,WAVEL,FNAME,WLFAC, 
     1             REFMED,REFRE1,REFIM1,REFRE2,REFIM2)
      PRINT *, "Wavelength is: ", WAVEL
      PRINT *, "using this nk data for reference medium: ", FNAME(1)
      PRINT *, "got: ",REFMED," + i*",0.00," for ",WAVEL," micron"
      PRINT *, "using this nk data for core medium: ", FNAME(2)
      PRINT *, "got: ",REFRE1," + i*",REFIM1," for ",WAVEL," micron"
      PRINT *, "using this nk data for shell medium: ", FNAME(3)
      PRINT *, "got: ",REFRE2," + i*",REFIM2," for ",WAVEL," micron"
      
      
      PRINT *, "initializing solar intensity data"
      IFILNAM="ASTMG173.csv"
      IFNAME=DIRNAM(1:MAX(INDEX(DIRNAM,' ')-1,1))//
     1          IFILNAM(1:MAX(INDEX(IFILNAM,' ')-1,1))
      CALL SOLARINTENSITY(1,WAVEL,IFNAME,INTENSITY,2)
      
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
      
      PRINT *, "Here it is: ", ANS
      
      
      CALL SOLARINTENSITY(0,WAVEL,IFNAME,INTENSITY,2)
      PRINT *, "And here is the intensity at ", WAVEL, "um"
      PRINT *, "  " , INTENSITY
      
      TEMP = 2911.0D0
      CALL BBEPOW(WAVEL,TEMP,IB)
      PRINT *, "AND now the blackbody emmisive power there at T=1000"
      PRINT *, "  ", IB
      
      CALL BBEPOWDT(WAVEL,TEMP,IB)
      PRINT *, "AND now the blackbody emmisive power DT there at T=1000"
      PRINT *, "  ", IB
      
      
      !call cubacores(0,1000)
      !call initfun(diffIsolPabs,32768)
      !CALL IntegrateBB(ANS,TEMP,1)
      PRINT *,"integral over whole spectrum should give sigT^4"
      PRINT *, "MINE: ", ANS
      PRINT *, "Should be: ", 5.670E-8*TEMP**4.0D0
      PRINT *, "check: ", ANS/(5.670E-8*TEMP**4.0D0)
      
      npar = 1
      params(1) = TEMP
      
      
      CALL IntegrateGeneric(ANS,diffIsolPabs,params,npar,1)
      PRINT *, ""
      PRINT *, "Integrating over solar spec. should get ~1000W/m^2"
      PRINT *, "How did I do?"
      PRINT *, ANS
      
      !CALL IntegrateGeneric(ANS,bbfunc,params,npar,4)
      PRINT *, ""
      PRINT *, "This thingy should should give me the same result as "
      PRINT *, "IntegrateBB: ", ANS
      
      !params(1) = 5777.0D0
      !CALL IntegrateGeneric(ANS,bbint,params,npar,1)
      
      !print *, "check: ", ANS/(5.670E-8*5777.0D0**4.0D0/pi)
      
      stop
      end
      

  
  
C ********************************************************************
C integral over blackbody at given temp using 4 choices for integration scheme
C ********************************************************************
C   choice = 1 --> Vegas
C   choice = 2 --> Suave  
C   choice = 3 --> Divonne  
C   choice = 4 --> Cuhre
  
      SUBROUTINE IntegrateGeneric(ANS,extfunc,userdata,npar,choice)

      implicit none

      integer ndim, ncomp, nvec, last, seed, mineval, maxeval
      cubareal epsrel, epsabs
      
      parameter (ndim = 3)
      parameter (ncomp = 1)
      parameter (nvec = 1)
      parameter (epsrel = 1D-5)
      parameter (epsabs = 1D-12)
      parameter (last = 4)
      parameter (seed = 0)
      parameter (mineval = 0)
      parameter (maxeval = 100000)

      integer nstart, nincrease, nbatch, gridno,choice
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
      
      integer npar
      real*8 userdata(npar)
      
      REAL*8 ANS

      external extfunc

      cubareal integral(ncomp), error(ncomp), prob(ncomp)
      integer verbose, neval, fail, nregions
      character*16 env

      integer c
     
      !this is where I impose a serial eval... hopefully not too slow...
      
      
      call getenv("CUBAVERBOSE", env)
      verbose = 2
      read(env, *, iostat=fail, end=999, err=999) verbose
  999 continue


      if (choice == 1) then
        print *, "----- using Vegas -----"
        call Vegas(ndim, ncomp, extfunc, userdata, nvec, 
     &    epsrel, epsabs, verbose, seed, 
     &    mineval, maxeval, nstart, nincrease, nbatch, 
     &    gridno, statefile, spin, 
     &    neval, fail, integral, error, prob)

        print *, "neval    =", neval
        print *, "fail     =", fail
        print '(F20.12," +- ",F20.12,"   p = ",F8.3)', 
     &     (integral(c), error(c), prob(c), c = 1, ncomp)
      
      else if (choice == 2) then
        print *, "----- using Suave -----"
        call suave(ndim, ncomp, extfunc, userdata, nvec,
     &    epsrel, epsabs, verbose + last, seed,
     &    mineval, maxeval, nnew, nmin, flatness,
     &    statefile, spin,
     &    nregions, neval, fail, integral, error, prob)

        print *, "nregions =", nregions
        print *, "neval    =", neval
        print *, "fail     =", fail
        print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)
      
      else if (choice == 3) then
        print *, "----- using Divonne -----"
        call divonne(ndim, ncomp, extfunc, userdata, nvec,
     &    epsrel, epsabs, verbose, seed,
     &    mineval, maxeval, key1, key2, key3, maxpass,
     &    border, maxchisq, mindeviation,
     &    ngiven, ldxgiven, 0, nextra, 0,
     &    statefile, spin,
     &    nregions, neval, fail, integral, error, prob)

        print *, "nregions =", nregions
        print *, "neval    =", neval
        print *, "fail     =", fail
        print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)
      
      else if (choice == 4) then
        print *, "----- using Cuhre -----"
        call cuhre(ndim, ncomp, extfunc, userdata, nvec,
     &    epsrel, epsabs, verbose + last,
     &    mineval, maxeval, key,
     &    statefile, spin,
     &    nregions, neval, fail, integral, error, prob)

        print *, "nregions =", nregions
        print *, "neval    =", neval
        print *, "fail     =", fail
        print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)
      
      end if
      
      
      
      ANS = integral(1)

      RETURN
      END 
        
C ********************************************************************
C This is the integrand I_sol*p_abs to be integrated over angle and wavelength
C ********************************************************************

      integer function diffIsolPabs(ndim,x,ncomp,f,userdata,nvec,core)
      implicit none
      real*8 pi, solint
      real*8 WLFAC(3),REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      CHARACTER*80 DIRNAM,FILNAM(3),FNAME(3),IFNAME,IFILNAM
      integer ndim, ncomp, I, core, nvec
      double precision x(ndim), f(ncomp), th,ph,lam
      real*8 userdata


      DIRNAM=DATA_DIR
      
      pi=acos(-1.0D0)
      
      th = x(1)*pi
      ph = x(2)*2.0D0*pi
      
      lam = 0.28D0 + (4.0D0-0.28D0)*x(3) !shift to integral over solar spectrum

      
      
      ! note first chunk is necessary to change integration interval
      ! second 10^-6 is a Jacobian since we changed the integral from meters to micron
      ! f = (2.0D0*pi*pi* (4.0D0-0.28D0) ) * th*ph*lam*1.0D-6*1.0D-6
      ! checks out with analytical one given by mathematica
      
C This is ugly.  May want to make it so the file is set only once
C  in the entire program
      
C     The optical constants:
C     reference (background)
      FILNAM(1)='vac.nk'
      WLFAC(1)=1.0D-4 !factor so all wavelengths in micron
C     core
      FILNAM(2)='Au_babar.nk'
      WLFAC(2)=1.0D0
C     shell
      FILNAM(3)='CeO2_patsalas.nk'
      WLFAC(3)=1.0D-3
      DO 901 I=1,3
       FNAME(I)=DIRNAM(1:MAX(INDEX(DIRNAM,' ')-1,1))//
     1          FILNAM(I)(1:MAX(INDEX(FILNAM(I),' ')-1,1))
  901 CONTINUE
  
      CALL OPTCON(0,lam,FNAME,WLFAC, 
     1             REFMED,REFRE1,REFIM1,REFRE2,REFIM2)

      
      IFILNAM="ASTMG173.csv"
      IFNAME=DIRNAM(1:MAX(INDEX(DIRNAM,' ')-1,1))//
     1          IFILNAM(1:MAX(INDEX(IFILNAM,' ')-1,1))
      CALL SOLARINTENSITY(0,lam,IFNAME,solint,2)

  229 FORMAT('**pid ' I1, 
     & '**  |wavelength:', E12.5, 'um',  
     & '| |core n+ik:', E12.5, '+i',E12.5, 
     & '| |shell n+ik:', E12.5, '+i',E12.5,
     & '| |intensity:',E12.5,'|')
      !WRITE(*,229) core,lam,REFRE1,REFIM1,REFRE2,REFIM2,solint

      
      !test... set wavelength, but have bhcoat called here to see
      ! if the final result is fuckered.
      
      ! power of three since we are integrating over meter (10^-6) and
      ! intensity is given in per nm (10^9)
      f = (4.0D0-0.28D0 )* solint*1.0D3
      !f = 1
      
      diffIsolPabs = 0
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
      
      SUBROUTINE BBEPOW(lambda, t, eb) 

      implicit none


      real*8 lambda, t, lambda_m
      real*8 eb
      REAL*8 c_first_radiation
      REAL*8 c_second_radiation
      ! -- Local declarations --

      real*8 beta
      
      print *, "did I make it?"
      
      lambda_m = lambda*1.0D-6 !converstion from micron to meter
      !lambda_m = lambda!converstion from micron to meter
      
      c_first_radiation = 3.74177118e-16
      c_second_radiation = 14387.75225e-06
      beta  = c_second_radiation / (lambda_m * t)
      eb = c_first_radiation  / (lambda_m ** 5 * (exp(beta) - 1.0D0))

      

      return
      end
      
      
      SUBROUTINE BBEPOWM(lambda, t, eb) 

      implicit none


      real*8 lambda, t, lambda_m
      real*8 eb
      REAL*8 c_first_radiation
      REAL*8 c_second_radiation
      ! -- Local declarations --

      real*8 beta
      
      !lambda_m = lambda*1.0D-6 !converstion from micron to meter
      lambda_m = lambda*1.0D0!NO converstion from micron to meter
      
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
      
      SUBROUTINE BBEPOWDT(lambda, t, ebdt) 

      implicit none


      real*8 lambda, t, lambda_m
      real*8 ebdt
      REAL*8 c1
      REAL*8 c2
      ! -- Local declarations --

      
      lambda_m = lambda*1.0D-6 !converstion from micron to meter
      
      c1 = 3.74177118e-16
      c2 = 14387.75225e-06

      ebdt = (c1*c2*EXP(c2/(t*lambda_m)))
      ebdt = ebdt/((-1.0D0+EXP(c2/(t*lambda_m)))**2*t**2*lambda_m**6)
    
      return
      end
      
      SUBROUTINE BBEPOWDTM(lambda, t, ebdt) 

      implicit none


      real*8 lambda, t, lambda_m
      real*8 ebdt
      REAL*8 c1
      REAL*8 c2
      ! -- Local declarations --

      
      lambda_m = lambda*1.0D0 !NO converstion from micron to meter
      
      c1 = 3.74177118e-16
      c2 = 14387.75225e-06

      ebdt = (c1*c2*EXP(c2/(t*lambda_m)))
      ebdt = ebdt/((-1.0D0+EXP(c2/(t*lambda_m)))**2*t**2*lambda_m**6)
    
      return
      end

      
C ********************************************************************
C integral over blackbody at given temp using 4 choices for integration scheme
C ********************************************************************
C   choice = 1 --> Vegas
C   choice = 2 --> Suave  
C   choice = 3 --> Divonne  
C   choice = 4 --> Cuhre

      SUBROUTINE IntegrateBB(ANS,TEMP,choice)

      implicit none

      integer ndim, ncomp, nvec, last, seed, mineval, maxeval
      cubareal epsrel, epsabs
      
      parameter (ndim = 3)
      parameter (ncomp = 1)
      parameter (nvec = 1)
      parameter (epsrel = 1D-5)
      parameter (epsabs = 1D-12)
      parameter (last = 4)
      parameter (seed = 0)
      parameter (mineval = 0)
      parameter (maxeval = 100000)

      integer nstart, nincrease, nbatch, gridno,choice
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
      
      REAL*8 ANS, userdata,TEMP

      external bbfunc

      cubareal integral(ncomp), error(ncomp), prob(ncomp)
      integer verbose, neval, fail, nregions
      character*16 env

      integer c

      
      userdata = TEMP
      
      call getenv("CUBAVERBOSE", env)
      verbose = 2
      read(env, *, iostat=fail, end=999, err=999) verbose
  999 continue


      if (choice == 1) then
        print *, "----- using Vegas -----"
        call Vegas(ndim, ncomp, bbfunc, userdata, nvec, 
     &    epsrel, epsabs, verbose, seed, 
     &    mineval, maxeval, nstart, nincrease, nbatch, 
     &    gridno, statefile, spin, 
     &    neval, fail, integral, error, prob)

        print *, "neval    =", neval
        print *, "fail     =", fail
        print '(F20.12," +- ",F20.12,"   p = ",F8.3)', 
     &     (integral(c), error(c), prob(c), c = 1, ncomp)
      
      else if (choice == 2) then
        print *, "----- using Suave -----"
        call suave(ndim, ncomp, bbfunc, userdata, nvec,
     &    epsrel, epsabs, verbose + last, seed,
     &    mineval, maxeval, nnew, nmin, flatness,
     &    statefile, spin,
     &    nregions, neval, fail, integral, error, prob)

        print *, "nregions =", nregions
        print *, "neval    =", neval
        print *, "fail     =", fail
        print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)
      
      else if (choice == 3) then
        print *, "----- using Divonne -----"
        call divonne(ndim, ncomp, bbfunc, userdata, nvec,
     &    epsrel, epsabs, verbose, seed,
     &    mineval, maxeval, key1, key2, key3, maxpass,
     &    border, maxchisq, mindeviation,
     &    ngiven, ldxgiven, 0, nextra, 0,
     &    statefile, spin,
     &    nregions, neval, fail, integral, error, prob)

        print *, "nregions =", nregions
        print *, "neval    =", neval
        print *, "fail     =", fail
        print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)
      
      else if (choice == 4) then
        print *, "----- using Cuhre -----"
        call cuhre(ndim, ncomp, bbfunc, userdata, nvec,
     &    epsrel, epsabs, verbose + last,
     &    mineval, maxeval, key,
     &    statefile, spin,
     &    nregions, neval, fail, integral, error, prob)

        print *, "nregions =", nregions
        print *, "neval    =", neval
        print *, "fail     =", fail
        print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &    (integral(c), error(c), prob(c), c = 1, ncomp)
      
      end if
      
      
      
      ANS = integral(1)

      RETURN
      END 
      

C**************
C blackbody emmisive power function in format for cuba
C ******************
        integer function bbfunc(ndim, x, ncomp, f,userdata)
        implicit none
        integer ndim, ncomp
        cubareal x(ndim), f(ncomp),xi
        real*8 t
        REAL*8 c1
        REAL*8 c2
        ! -- Local declarations --

        real*8 userdata
        

        !converstion from micron to meter and transformation to 0-1 interval
        !lambda_m = x/(1.0D0-x) 
        
        t = userdata
C        CALL BBEPOWM(lambda_m,t, eb)
C        f = eb/(1.0D0-x)/(1.0D0-x)

        c1 = 3.74177118e-16
        c2 = 14387.75225e-06
        
        xi = x(1)/(1.0D0-x(1)) 
        
        !this gives sig T4
        f = xi**3.0D0/(EXP(xi)-1.0D0)/(1.0D0-x(1))/(1.0D0-x(1))
        f = f*c1/c2**4.0D0 * t**4.0D0

        !test to see if math works (of course it does, agrees with mathematica)
        !f = xi**3.0D0/(EXP(xi)-1.0D0)/(1.0D0-x)/(1.0D0-x)
        !f = f*c1/c2**4.0D0 *4.0D0 * t**3.0D0

C        f = SIN(x)*EXP(-x)
C        xp = x/(1.0D0-x)
C        f = SIN(xp)*EXP(-xp)/(1.0D0-x)/(1.0D0-x)

        bbfunc = 0
        end


C**************
C blackbody intensity function in format for cuba
C ******************
        integer function bbint(ndim, x, ncomp, f,userdata)
        implicit none
        integer ndim, ncomp
        cubareal x(ndim), f(ncomp),xi

        real*8 t,pi
        REAL*8 c1
        REAL*8 c2
        ! -- Local declarations --

        real*8 userdata
        
        !converstion from micron to meter and transformation to 0-1 interval
        !lambda_m = x/(1.0D0-x) 
        
        pi=ACOS(-1.0D0)
        
        t = userdata

        c1 = 3.74177118e-16
        c2 = 14387.75225e-06
        
        xi = x(1)/(1.0D0-x(1)) 
        
        !this gives sig T4
        f = xi**3.0D0/(EXP(xi)-1.0D0)/(1.0D0-x(1))/(1.0D0-x(1))
        f = f*c1/c2**4.0D0 * t**4.0D0 / pi


        bbint = 0
        end
        

C************************************************************************

        integer function sinexp(ndim, x, ncomp, f)
        implicit none
        integer ndim, ncomp
        cubareal x(ndim), f(ncomp),xp
        
C        f = SIN(x)*EXP(-x)
        xp = x(1)/(1.0D0-x(1))
        f = SIN(xp)*EXP(-xp)/(1.0D0-x(1))/(1.0D0-x(1))

        sinexp = 0
        end
  
