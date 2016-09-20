#define BHFIELD_VERSION 'Oct 5, 2012'
#define cubareal real*8    

      PROGRAM create_files_for_FEMsim

      IMPLICIT NONE
      REAL*8 UFTOL
      PARAMETER (UFTOL=1.0D-6)
      
      INTEGER I,nshells,ntemps
      
      REAL*8 WLFAC(3)
      
      REAL*8 RADCOR,RADCOT,WAVEL,PI,CC,EPSVAC,MU
      REAL*8 REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      REAL*8 INTENSITY,TEMP,IB!,rad
            
      CHARACTER*80 DIRNAM,FILNAM(3),FNAME(3),IFNAME,IFILNAM
      CHARACTER*50 FNLOGF
      
      integer npar_s
      real*8 userdata_s(7)
      external diffIsolPabs, diffShellIsolPabs, bbfunc, bbint
      
      
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
      WAVEL = 0.507D0 ! dummy wavelength to initialize datafiles
      RADCOR = 0.03D0
      RADCOT = 0.06D0
      PI=ACOS(-1.0D0)
      CC=2.99792458D8 ! light speed [m s-1]
      EPSVAC=1.0D7/(4.0D0*PI*CC*CC) ! eps0[F m-1]
      MU=4.0D0*PI*1.0D-7 ! assume non-magnetic (MU=MU0=const) [N A-2] 
      !OMEGA=2.0D0*PI*CC/(WAVEL*1.0D-6) ! angular frequency [s-1]
      
      
C     log file
C      FNLOGF='bhfield.log'
C      OPEN(51,FILE=FNLOGF,STATUS='UNKNOWN')
      
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

      
      CALL SOLARINTENSITY(0,WAVEL,IFNAME,INTENSITY,2)
      PRINT *, "And here is the intensity at ", WAVEL, "um"
      PRINT *, "  " , INTENSITY
      
      TEMP = 1000.0D0
      CALL BBINT(WAVEL,TEMP,IB)
      PRINT *, "blackbody emmisive power there at T=",TEMP
      PRINT *, "  ", IB
      
      CALL BBINTDT(WAVEL,TEMP,IB)
      PRINT *, "blackbody emmisive power DT there at T=",TEMP
      PRINT *, "  ", IB

       
C      rad = 0.025
      userdata_s(1) = RADCOR
      userdata_s(2) = RADCOT
      userdata_s(3) = EPSVAC
      userdata_s(4) = CC
      userdata_s(5) = MU
C      userdata_s(6) = rad !set in subroutine
      npar_s = 7
      print *, "starting the business"
      nshells = 50
      ntemps = 50
      CALL IntegrateShellsFull(userdata_s,npar_s,nshells,ntemps)

   
      
      stop
      end
      
      
C **********************************************************************
C Should be final subroutine that produces files for MATLAB
C **********************************************************************
      
      SUBROUTINE IntegrateShellsFull(userdata,npar,nshells,ntemps)
      
      implicit none
      
      integer nshells,npar,i,j,k
      integer ntemps
      real*8 userdata(npar)
      real*8 rads(nshells),uabs(nshells),temps(ntemps)
      real*8 uem(nshells,ntemps),uemdt(nshells,ntemps)
      real*8 radcor,radcot,rstep,t,delt,tlow,thigh
      
      character(3) crad, srad
      character(5) csrad
      
      external diffShellIsolPabs,diffShellIbbPabs
      external diffShellIbbdTPabs
      
      radcor = userdata(1)
      radcot = userdata(2)
      rstep = radcot/dble(nshells-1)
      
      ! first do the absorption at each radial distance
      do 291 i=1,nshells
        rads(i) = dble(i-1)*rstep
        userdata(6) = rads(i)
        call IntegrateGeneric(uabs(i),diffShellIsolPabs,
     &              userdata,npar,3,1)
  291 continue
      
      
      tlow = 300
      thigh = 2000
      delt = (thigh-tlow)/dble(ntemps-1)
      t = tlow
      !then do the emission and temp deriv of emission
      do 292 j=1,ntemps
        !rads(i) = dble(i-1)*rstep
        userdata(7) = t
        temps(j)=t
        print *,delt
        do 293 i=1,nshells
          userdata(6) = rads(i)
          call IntegrateGeneric(uem(i,j),diffShellIbbPabs,
     &              userdata,npar,3,1)
          call IntegrateGeneric(uemdt(i,j),diffShellIbbdTPabs,
     &              userdata,npar,3,1)
  293   continue
        t = t + delt
  292 continue
      
      

  701 FORMAT(3E13.5,E13.5)
  702 FORMAT(3E13.5)
  
C     make some file names according to particle size
      write(crad,'(I2.0)') NINT(radcor*1.E3)
      write(srad,'(I2.0)') NINT(radcot*1.E3)
      csrad = trim(crad)//'-'//trim(srad)
      print *, '=================================== '
      print *, csrad
      print *, '=================================== '
  
      
      OPEN(44,FILE=csrad//'_'//'UabsVr.dat',STATUS='UNKNOWN')
      WRITE(44,*) 'Absorption per volume'
      WRITE(44,*) '  units [=] (W m^-3)'
      WRITE(44,*) '-----------'
C      
      OPEN(45,FILE=csrad//'_'//'UemVr.dat',STATUS='UNKNOWN')
      WRITE(45,*) 'Emission per volume'
      WRITE(45,*) '  units [=] (W m^-3)'
      WRITE(45,'(A13)',advance='no') 'temperatures:'
      do 123 k=1,ntemps
        WRITE(45,702,advance='no') temps(k)
  123 continue
      WRITE(45,*)
      WRITE(45,*) '-----------'
C      
      OPEN(46,FILE=csrad//'_'//'UemdTVr.dat',STATUS='UNKNOWN')
      WRITE(46,*) 'Emission per volume'
      WRITE(46,*) '  units [=] (W m^-3)'
      WRITE(46,'(A13)',advance='no') 'temperatures:'
      do 124 k=1,ntemps
        WRITE(46,702,advance='no') temps(k)
  124 continue
      WRITE(46,*)
      WRITE(46,*) '-----------'
      
      
      do 119 i=1,nshells
        WRITE(44,701) rads(i)*1.0E-6, uabs(i)
  119 continue   
  
      
      do 120 i=1,nshells
        WRITE(45,702,advance='no') rads(i)
        WRITE(46,702,advance='no') rads(i)
        do 121 j = 1,ntemps
          WRITE(45,702,advance='no') uem(i,j)
          WRITE(46,702,advance='no') uemdt(i,j)
  121   continue
        WRITE(45,*)
        WRITE(46,*)
  120 continue
      
  
      return
      end
  
  
C ********************************************************************
C This is the integrand I_bb*p_abs to be integrated over angle and 
C wavelength only at a given radius
C ********************************************************************

      integer function diffShellIbbPabs(ndim,x,
     &                        ncomp,f,userdata,nvec,pid)
      implicit none
      
      real*8 pi,ibb,t
      real*8 WLFAC(3),REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      real*8 QEXT,QSCA,QABS,QBACK,XP(3),UABS,MU,OMEGA,CC
      real*8 EFSQ,I0,EPSVAC,rad,lammin
      real*8 RADCOT,RADCOR,RMAX,EXTMAX,Xpara,Ypara,Y1,Y2,Y3,Y4,YMAX
      
      CHARACTER*80 DIRNAM,FILNAM(3),FNAME(3)
      
      integer ndim,ncomp,I,pid,nvec,NSTOPF,IWHERE
      
      double precision x(ndim), f(ncomp),lam
      real*8 userdata(7)

      COMPLEX*16 RFREL1,RFREL2,EC(3),HC(3)

      DIRNAM=DATA_DIR
      
      RADCOR = userdata(1)
      RADCOT = userdata(2)
      EPSVAC = userdata(3)
      CC = userdata(4)
      MU = userdata(5)
      rad = userdata(6)
      t = userdata(7)

      
      pi=acos(-1.0D0)
      lammin = 0.28D0
      lam = lammin + x(1)/(1.0D0-x(1))  !shift to integral over infinite interval

      
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

      RFREL1=DCMPLX(REFRE1,REFIM1)/REFMED
      RFREL2=DCMPLX(REFRE2,REFIM2)/REFMED
        
      Xpara=2.0D0*PI*RADCOR*REFMED/lam
      Ypara=2.0D0*PI*RADCOT*REFMED/lam
       
      Y1=2.0D0*PI*RADCOR*ABS(DCMPLX(REFRE1,REFIM1))/lam
      Y2=2.0D0*PI*RADCOT*ABS(DCMPLX(REFRE1,REFIM1))/lam
      Y3=2.0D0*PI*RADCOT*ABS(DCMPLX(REFRE2,REFIM2))/lam
    
      EXTMAX=RADCOT
      RMAX=EXTMAX*SQRT(3.0D0)
      Y4=2.0D0*PI*RMAX*REFMED/lam
      YMAX=MAX(Y1,Y2,Y3,Y4)
      NSTOPF=INT(YMAX+4.05D0*YMAX**0.3333D0+2.0D0)


      CALL BHCOAT(Xpara,Ypara,RFREL1,RFREL2,NSTOPF,QEXT,QSCA,QBACK)
      QABS=QEXT-QSCA
      
      XP(1) = rad
      XP(2) = pi*x(2) !theta
      XP(3) = 2.0D0*pi*x(3) !phi
     
      CALL FIELDVMW(lam,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1                 RADCOR,RADCOT,XP,IWHERE,EC,HC)
      I0 = 0.5*SQRT(EPSVAC/MU)*1.
      EFSQ=ABS(EC(1))**2.0D0+ABS(EC(2))**2.0D0+ABS(EC(3))**2.0D0
      OMEGA=2.0D0*PI*CC/(lam*1.0D-6) ! angular frequency [s-1]
      IF(IWHERE.EQ.1) THEN
        UABS=EPSVAC*OMEGA*REFRE1*REFIM1*EFSQ/I0
      ELSE IF(IWHERE.EQ.2) THEN
        UABS=EPSVAC*OMEGA*REFRE2*REFIM2*EFSQ/I0
      ELSE
        UABS=0.0D0
      END IF

      call BBINT(lam,t,ibb)
      
      !PRINT *, lam,t,ibb
      
C*****finalmente
      f(1) = 2.0D0*pi*pi*ibb*UABS*sin(pi*x(2))
      f(1) = f(1)/(1.0D0-x(1))/(1.0D0-x(1)) !infinite interval shift
C *****this confirms that integration over infinite interval works for ibb
C      f(1) = ibb/(1.0D0-x(1))/(1.0D0-x(1))
      
      diffShellIbbPabs = 0
      end
      
      
C ********************************************************************
C This is the integrand d(I_bb)/dT*p_abs to be integrated over angle and 
C wavelength only at a given radius
C ********************************************************************

      integer function diffShellIbbdTPabs(ndim,x,
     &                        ncomp,f,userdata,nvec,pid)
      implicit none
      
      real*8 pi,ibbdt,t
      real*8 WLFAC(3),REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      real*8 QEXT,QSCA,QABS,QBACK,XP(3),UABS,MU,OMEGA,CC
      real*8 EFSQ,I0,EPSVAC,rad,lammin
      real*8 RADCOT,RADCOR,RMAX,EXTMAX,Xpara,Ypara,Y1,Y2,Y3,Y4,YMAX
      
      CHARACTER*80 DIRNAM,FILNAM(3),FNAME(3)
      
      integer ndim,ncomp,I,pid,nvec,NSTOPF,IWHERE
      
      double precision x(ndim), f(ncomp),lam
      real*8 userdata(7)

      COMPLEX*16 RFREL1,RFREL2,EC(3),HC(3)

      DIRNAM=DATA_DIR
      
      RADCOR = userdata(1)
      RADCOT = userdata(2)
      EPSVAC = userdata(3)
      CC = userdata(4)
      MU = userdata(5)
      rad = userdata(6)
      t = userdata(7)

      
      pi=acos(-1.0D0)
      lammin = 0.28D0
      lam = lammin + x(1)/(1.0D0-x(1))  !shift to integral over infinite interval

      
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

      RFREL1=DCMPLX(REFRE1,REFIM1)/REFMED
      RFREL2=DCMPLX(REFRE2,REFIM2)/REFMED
        
      Xpara=2.0D0*PI*RADCOR*REFMED/lam
      Ypara=2.0D0*PI*RADCOT*REFMED/lam
       
      Y1=2.0D0*PI*RADCOR*ABS(DCMPLX(REFRE1,REFIM1))/lam
      Y2=2.0D0*PI*RADCOT*ABS(DCMPLX(REFRE1,REFIM1))/lam
      Y3=2.0D0*PI*RADCOT*ABS(DCMPLX(REFRE2,REFIM2))/lam
    
      EXTMAX=RADCOT
      RMAX=EXTMAX*SQRT(3.0D0)
      Y4=2.0D0*PI*RMAX*REFMED/lam
      YMAX=MAX(Y1,Y2,Y3,Y4)
      NSTOPF=INT(YMAX+4.05D0*YMAX**0.3333D0+2.0D0)


      CALL BHCOAT(Xpara,Ypara,RFREL1,RFREL2,NSTOPF,QEXT,QSCA,QBACK)
      QABS=QEXT-QSCA
      
      XP(1) = rad
      XP(2) = pi*x(2) !theta
      XP(3) = 2.0D0*pi*x(3) !phi
     
      CALL FIELDVMW(lam,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1                 RADCOR,RADCOT,XP,IWHERE,EC,HC)
      I0 = 0.5*SQRT(EPSVAC/MU)*1.
      EFSQ=ABS(EC(1))**2.0D0+ABS(EC(2))**2.0D0+ABS(EC(3))**2.0D0
      OMEGA=2.0D0*PI*CC/(lam*1.0D-6) ! angular frequency [s-1]
      IF(IWHERE.EQ.1) THEN
        UABS=EPSVAC*OMEGA*REFRE1*REFIM1*EFSQ/I0
      ELSE IF(IWHERE.EQ.2) THEN
        UABS=EPSVAC*OMEGA*REFRE2*REFIM2*EFSQ/I0
      ELSE
        UABS=0.0D0
      END IF

      call BBINTDT(lam,t,ibbdt)
      
      !PRINT *, lam,t,ibb
      
C*****finalmente
      f(1) = 2.0D0*pi*pi*ibbdt*UABS*sin(pi*x(2))
      f(1) = f(1)/(1.0D0-x(1))/(1.0D0-x(1)) !infinite interval shift
C *****this confirms that integration over infinite interval works for ibb
C      f(1) = ibb/(1.0D0-x(1))/(1.0D0-x(1))
      
      diffShellIbbdTPabs = 0
      end
  

        
C ********************************************************************
C This is the integrand I_sol*p_abs to be integrated over angle and 
C wavelength only at a given radius
C ********************************************************************

      integer function diffShellIsolPabs(ndim,x,
     &                        ncomp,f,userdata,nvec,pid)
      implicit none
      
      real*8 pi, solint
      real*8 WLFAC(3),REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      real*8 QEXT,QSCA,QABS,QBACK,XP(3),UABS,MU,OMEGA,CC
      real*8 EFSQ,I0,EPSVAC,lammin,dlam,rad
      real*8 RADCOT,RADCOR,RMAX,EXTMAX,Xpara,Ypara,Y1,Y2,Y3,Y4,YMAX
      
      CHARACTER*80 DIRNAM,FILNAM(3),FNAME(3),IFNAME,IFILNAM
      
      integer ndim,ncomp,I,pid,nvec,NSTOPF,IWHERE
      
      double precision x(ndim), f(ncomp),lam
      real*8 userdata(6)

      COMPLEX*16 RFREL1,RFREL2,EC(3),HC(3)

      DIRNAM=DATA_DIR
      
      RADCOR = userdata(1)
      RADCOT = userdata(2)
      EPSVAC = userdata(3)
      CC = userdata(4)
      MU = userdata(5)
      rad = userdata(6)

      
      pi=acos(-1.0D0)

      lammin = 0.28D0
      dlam = 4.0D0-0.28D0
      lam = lammin + dlam*x(1) !shift to integral over solar spectrum

      
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

      
      ! power of three since we are integrating over meter (10^-6) and
      ! intensity is given in per nm (10^9)
      ! uncomment this f to get 1000 W/m/m(as you should)
      !f = (4.0D0-0.28D0 )* solint*1.0D3
      !f = 1
      RFREL1=DCMPLX(REFRE1,REFIM1)/REFMED
      RFREL2=DCMPLX(REFRE2,REFIM2)/REFMED
        
      Xpara=2.0D0*PI*RADCOR*REFMED/lam
      Ypara=2.0D0*PI*RADCOT*REFMED/lam
       
      Y1=2.0D0*PI*RADCOR*ABS(DCMPLX(REFRE1,REFIM1))/lam
      Y2=2.0D0*PI*RADCOT*ABS(DCMPLX(REFRE1,REFIM1))/lam
      Y3=2.0D0*PI*RADCOT*ABS(DCMPLX(REFRE2,REFIM2))/lam
    
      EXTMAX=RADCOT
      RMAX=EXTMAX*SQRT(3.0D0)
      Y4=2.0D0*PI*RMAX*REFMED/lam
      YMAX=MAX(Y1,Y2,Y3,Y4)
      NSTOPF=INT(YMAX+4.05D0*YMAX**0.3333D0+2.0D0)


      CALL BHCOAT(Xpara,Ypara,RFREL1,RFREL2,NSTOPF,QEXT,QSCA,QBACK)
      QABS=QEXT-QSCA
      
      XP(1) = rad
      XP(2) = pi*x(2) !theta
      XP(3) = 2.0D0*pi*x(3) !phi
     
      CALL FIELDVMW(lam,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1                 RADCOR,RADCOT,XP,IWHERE,EC,HC)
      I0 = 0.5*SQRT(EPSVAC/MU)*1.
      EFSQ=ABS(EC(1))**2.0D0+ABS(EC(2))**2.0D0+ABS(EC(3))**2.0D0
      OMEGA=2.0D0*PI*CC/(lam*1.0D-6) ! angular frequency [s-1]
      IF(IWHERE.EQ.1) THEN
        UABS=EPSVAC*OMEGA*REFRE1*REFIM1*EFSQ/I0
      ELSE IF(IWHERE.EQ.2) THEN
        UABS=EPSVAC*OMEGA*REFRE2*REFIM2*EFSQ/I0
      ELSE
        UABS=0.0D0
      END IF

C*****finalmente
      f(1) = dlam*2.0D0*pi*pi*solint*UABS*sin(pi*x(2))

      diffShellIsolPabs = 0
      end
        
      
      

        
C ********************************************************************
C This thing should return the ASTM solar spectrum 
C ********************************************************************
C   WHICH = 1 --> Extraterrestrial radiation
C   WHICH = 2 --> Global tilt
C   WHICH = 3 --> Direct + circumsolar
C *NB: raw data in per nm, but returns in per micron

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
       INTENSITY=INTENSITY*1.0D3
C
C
      RETURN
      END
      
      
C ***********************************************************************
C Function to calculate blackbody emissive power at a given temperature
C  as a function of wavelength (in micron)
C ***********************************************************************
      
      SUBROUTINE BBINT(lambda, t, eb) 

      implicit none


      real*8 lambda, t, lambda_m
      real*8 eb,pi
      REAL*8 c_first_radiation
      REAL*8 c_second_radiation
      ! -- Local declarations --

      real*8 beta
      
      pi=ACOS(-1.0D0)
      
      lambda_m = lambda*1.0D-6 !converstion from micron to meter
      !lambda_m = lambda!converstion from micron to meter
      
      c_first_radiation = 3.74177118e-16
      c_second_radiation = 14387.75225e-06
      beta  = c_second_radiation / (lambda_m * t)
      eb = c_first_radiation  / (lambda_m ** 5 * (exp(beta) - 1.0D0))
      eb = eb/pi*1.0D-6 ! per meter to per micron
      

      return
      end
      

C ***********************************************************************
C Function to calculate temperature derivative of the 
C  blackbody emissive power at a given temperature
C  as a function of wavelength (in micron)
C ***********************************************************************
      
      SUBROUTINE BBINTDT(lambda, t, ebdt) 

      implicit none


      real*8 lambda, t, lambda_m
      real*8 ebdt,pi
      REAL*8 c1
      REAL*8 c2
      ! -- Local declarations --

      pi=ACOS(-1.0D0)
      
      lambda_m = lambda*1.0D-6 !converstion from micron to meter
      
      c1 = 3.74177118e-16
      c2 = 14387.75225e-06

      ebdt = (c1*c2*EXP(c2/(t*lambda_m)))
      ebdt = ebdt/((-1.0D0+EXP(c2/(t*lambda_m)))**2*t**2*lambda_m**6)
      ebdt = ebdt/pi*1.0D-6
    
      return
      end
      

C ********************************************************************
C integral over "extfunc" (input) using 4 choices for integration scheme
C ********************************************************************
C   choice = 1 --> Vegas
C   choice = 2 --> Suave  
C   choice = 3 --> Divonne  
C   choice = 4 --> Cuhre
  
      SUBROUTINE IntegrateGeneric(ANS,extfunc,userdata,npar,ndim,choice)

      implicit none

      integer ndim, ncomp, nvec, last, seed, mineval, maxeval
      cubareal epsrel, epsabs
      
      !parameter (ndim = 4)
      parameter (ncomp = 1)
      parameter (nvec = 1)
      parameter (epsrel = 1D-6)
      parameter (epsabs = 1D-20)
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
      !parameter (ldxgiven = ndim)
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
     
      ldxgiven = ndim
     
      
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
