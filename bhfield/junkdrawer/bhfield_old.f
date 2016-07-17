#define BHFIELD_VERSION 'Oct 5, 2012'
#define cubareal real*8    
C bhfield.f
C Calc internal & external intensity distributions of a nanoshell 
C by Mie theory; core & coat radii = const; wavelength = const
C for Ag/silica nanoshells, liposomes, etc.

C Author: Honoh Suzuki (based on Bohren & Huffman BHCOAT)
C Department of Chemistry, University of Toyama 2007-2012
C License: GPL
C See:
C H. Suzuki and I-Y. S. Lee: Calculation of the Mie Scattering Field 
C inside and outside a Coated Spherical Particle, 
C Int. J. Phys. Sci., 3, 38-41 (2008).
C Errata: Int. J. Phys. Sci. 4, 615 (2009).
C Eq. (9): The second term in the numerator should read 
C -m_2 x B_n XI_(n-1)(x_2) instead of -m_2 x B_n XI_n(x_2).  
C ("n-1" is the subscript for XI)

C Refs:
C Bohren & Huffman (BH)
C C. F. Bohren & D. R. Huffman, Absorption and Scattering of Light
C by Small Particles, Wiley Interscience (1983).
C p. 486; p. 92-95; p. 182
C Kreibig mean free path correction(Ag): BH p.337; Kreibig74
C U. Kreibig, J. Phys. F 4, 999 (1974).
C Barber & Hill
C P. W. Barber & S. C. Hill, Light Scattering by Particles: Computational
C Methods, World Scientific (1990).
C p. 189, 222, 228; Program S7, S8

C Assumptions:
C * non-absorbing medium (refmed = real)
C * nonmagnetic (mu = mu0); no net surface charge (BH p.118)
C * incident field: monochromatic, plane x-polarized (parallel) wave 
C   propagating in +z direction with unit amplitude (E0 = 1)

C Tweaks:
C At origin (R=0) in grid field computation, R=small is used to 
C avoid zero-division, but fake XP=YP=ZP=0 is retained for vtk 3D 
C plot.  

C Computational notes:
C Complex spherical Bessels (particularly y_n) and scattering 
C coefficients cast a serious round-off problem in the case of 
C highly absorbing material (e.g. silver).  It is mainly due to 
C the hyperbolic nature of complex trigonometry, and easily exceeds 
C the limit of standard Fortran precision (COMPLEX*16).  Such 
C problems can be detected by checking numerical equality of 
C algebraically equivalent expressions for c_n (and underflows 
C in numerator/denominator expressions).  If it happens, we should 
C use arbitrary precision package (LBNL ARPREC).  
C Get ARPREC here: http://crd.lbl.gov/~dhbailey/mpdist/

C Compiler preprocessing options:
C NMAX_COEFFS   max number of coefficients (mandatory)
C DATA_DIR      where complex refractive index data (n, k) reside
C USE_ARPREC    use arprec package to avoid round-offs (slow...)
C USE_ARPREC_LEGENDRE use arprec for Legendre functions too
C CHECK_UNDERFLOW    check underflow and show appropriate mpdigit
C CHECK_SURFACE_MODE check components of large electric fields 
C inside the coat field
C CHECK_NUMERICAL_POYNTING calc numerical derivative to estimate 
C -div<S> (1st-3rd) and DIVSR (very slow)
C CHECK_TANGENTIAL_CONTINUITY check tangential continuity 
C DEBUG_BESSEL  check spherical Bessels and roundoff errors
C OUTPUT_YBESSE output y_n(x) for roundoff error checking
C DEBUG_BESSEL_DERIV
C DEBUG_BH183


C arprec case (mingw/msys) example:
C gfortran -cpp -DNMAX_COEFFS=500 -DDATA_DIR="'./'" -DUSE_ARPREC 
C -DCHECK_SURFACE_MODE -O2 -Wall -c bhfield.f -I /usr/local/lib/arprec
C g++ -O2 -Wall -o bhfield.exe bhfield.o -static 
C -L/usr/local/lib -larprec_f_main -larprecmod -larprec 
C -Lc:/mingw/bin/../lib/gcc/mingw32/4.4.0 -Lc:/mingw/bin/../lib/gcc 
C -Lc:/mingw/bin/../lib/gcc/mingw32/4.4.0/../../../../mingw32/lib 
C -Lc:/mingw/bin/../lib/gcc/mingw32/4.4.0/../../.. -L/mingw/lib 
C -lgfortranbegin -lgfortran -lmingw32 -lmoldname -lmingwex -lmsvcrt 
C -luser32 -lkernel32 -ladvapi32 -lshell32

C arprec porting memo:
C (1) PROGRAM BHFIELD -> subroutine f_main
C (2) use mpmodule
C (3) REAL*8 -> type (mp_real); COMPLEX*16 -> type (mp_complex)
C (4) call mpinit(60) 
C     or integer mpdigit; parameter (mpdigit=20); call mpinit(mpdigit)
C (5) ARG (name reserved) -> SARG
C (6) function list: see arprec/fortran/mp_mod.f (interface)
C     COMPLEX(x, y), REAL/DBLE/REALPART(z), IMAG/IMAGPART(z), INT(x)
C     -> cmplx(x, y), mpreal(z), aimag(z), aint(x)
C     MAX() -> up to 3 args
C     CMPLX(z), CONJG(z): ok
C (7) CMPLX(x, 0.0D0) -> CMPLX(x, mpreal(0.0D0)) or CMPLX(x, ZERO)
C (8) READ(SARG,*) WAVEL -> use REAL*8 var for relay
C (9) (-1.0D0)**RN: MPLOG error -> mpreal((-1)**N)
C (10) WRITE() can't be used -> WRITE(*,*) DBLE(mpreal(z)) etc, or
C      substitute into ordinary variables just for write()
C (11) x=1.0D0 -> x='1.0D0'


      PROGRAM BHFIELD

      IMPLICIT NONE
      REAL*8 UFTOL
      PARAMETER (UFTOL=1.0D-6)
      
      INTEGER I,J,K,NSTOPF,IW,CT,NPTSPH,NPTSTH,IWHERE,N,NSHELLS
      
      REAL*8 WLFAC(3),UABS_T
      
      REAL*8 RADCOR,RADCOT,WAVEL,PI,CC,EPSVAC,MU,OMEGA
      REAL*8 REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      REAL*8 X,Y,Y1,Y2,Y3,Y4,YMAX,EXTMAX,RMAX
      REAL*8 QEXT,QSCA,QBACK,QABS,EFSQ,UABS,MYQABS
      REAL*8 DELTAS1,RAD, I0, AP
      
      COMPLEX*16 RFREL1,RFREL2
      COMPLEX*16 EC(3),HC(3)
      
      INTEGER RGRID,SGRID,NGRID(3)
      REAL*8 XPMIN(3),XPMAX(3),RSTEP,THSTEP,PHSTEP,XP(3)
      REAL*8 XCMIN(3),XCMAX(3),XCSTEP(3),XC(3)
      
      CHARACTER*80 DIRNAM,FILNAM(3),FNAME(3)
      CHARACTER*50 FNLOGF,FNESQ(4),FNUAB(4),FNUABC(4)
      CHARACTER*24 STRTIM
      CHARACTER*10 SFIELD(4)
      
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
      WAVEL = 0.50D0
      RADCOR = 0.03D0
      RADCOT = 0.06D0
      PI=ACOS(-1.0D0)
      CC=2.99792458D8 ! light speed [m s-1]
      EPSVAC=1.0D7/(4.0D0*PI*CC*CC) ! eps0[F m-1]
      MU=4.0D0*PI*1.0D-7 ! assume non-magnetic (MU=MU0=const) [N A-2] 
      OMEGA=2.0D0*PI*CC/(WAVEL*1.0D-6) ! angular frequency [s-1]
      
      
C     The optical constants:
C     reference (background)
      FILNAM(1)='vac.nk'
      WLFAC(1)=1.0D-4 !factor so all wavelengths in micron
C     core
C      FILNAM(2)='test_material.nk'
      FILNAM(2)='Ag_palik.nk'
      WLFAC(2)=1.0D-4
C     shell
C      FILNAM(3)='test_material.nk'
      FILNAM(3)='SiO2_palik.nk'
C      FILNAM(3)='Ag_palik.nk'
      WLFAC(3)=1.0D-4
      DO 901 I=1,3
       FNAME(I)=DIRNAM(1:MAX(INDEX(DIRNAM,' ')-1,1))//
     1          FILNAM(I)(1:MAX(INDEX(FILNAM(I),' ')-1,1))
  901 CONTINUE
      PRINT *, "using this nk data for reference medium: ", FNAME(1)
      PRINT *, "using this nk data for core medium: ", FNAME(2)
      PRINT *, "using this nk data for shell medium: ", FNAME(3)
      CALL OPTCON(1,WAVEL,FNAME,WLFAC, 
     1             REFMED,REFRE1,REFIM1,REFRE2,REFIM2)
      
C     get relative refractive indices
      RFREL1=DCMPLX(REFRE1,REFIM1)/REFMED
      RFREL2=DCMPLX(REFRE2,REFIM2)/REFMED
  
C     define the particle size parameters for core and shell
      X=2.0D0*PI*RADCOR*REFMED/WAVEL
      Y=2.0D0*PI*RADCOT*REFMED/WAVEL
      
      
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Define the grid
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RGRID = 10 !r
      SGRID = 5
      !THGRID = 5 !theta
      !PHGRID = 5 !phi
      
      XPMIN(1) = 0
      XPMAX(1) = RADCOT
      XPMIN(2) = 0
      XPMAX(2) = PI
      XPMIN(3) = 0
      XPMAX(3) = 2*PI
     
C     Now in cartesian     
      NGRID(1) = 4 
      NGRID(2) = 4
      NGRID(3) = 4
      
      XCMIN(1) = -RADCOT
      XCMAX(1) = RADCOT
      XCMIN(2) = -RADCOT
      XCMAX(2) = RADCOT
      XCMIN(3) = -RADCOT
      XCMAX(3) = RADCOT   
      

C
C output file names
C
      FNLOGF='bhfield.log'
      FNESQ(4)='E_0allf.dat'
      FNESQ(1)='E_1core.dat'
      FNESQ(2)='E_2coat.dat'
      FNESQ(3)='E_3exte.dat'
      FNUAB(4)='U_0allf_rad.dat'
      FNUAB(1)='U_1core_rad.dat'
      FNUAB(2)='U_2coat_rad.dat'
C non-absorbing medium: external UABS = 0
      FNUAB(3)='U_3exte_rad.dat'
      !Cartesian
      FNUABC(4)='U_0allf_cart.dat'
      FNUABC(1)='U_1core_cart.dat'
      FNUABC(2)='U_2coat_cart.dat'
C non-absorbing medium: external UABS = 0
      FNUABC(3)='U_3exte_cart.dat'
C

C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     Start log file and checking for underflow problems
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      OPEN(51,FILE=FNLOGF,STATUS='UNKNOWN') ! log file

#if defined(DEBUG_BESSEL) || defined(CHECK_UNDERFLOW)
      OPEN(88,FILE='bhdebug.log',STATUS='UNKNOWN')
      WRITE(88,*) 'standard version'
#endif /* DEBUG_BESSEL or CHECK_UNDERFLOW */
      CALL FDATE(STRTIM)
      WRITE(51,*) 'Started: ',STRTIM
      WRITE(51,*) ''
   29 FORMAT ("Optical properties:"/
     1"Ref index of ref medium: ",F7.4,3X/
     2"Ref index of core medium: ",F7.4,3X," + ",F7.4,3X,"*i"/
     3"Ref index of shell medium: ",F7.4,3X," + ",F7.4,3X,"*i")
      WRITE(51,29) REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      
   11 FORMAT (/"COATED SPHERE SCATTERING: bhfield (version: ",A,
     1        1X,A,")"/)
      WRITE(51,11) BHFIELD_VERSION,'standard'
   13 FORMAT ("Input parameters:"/
     1"WAVELENGTH [um] =",F7.4,3X,
     2"CORE RADIUS[um] =",F7.3,3X,"COAT RADIUS[um] =",F7.3/
     3"Rad-Span: ","Ngrid =",I4," min[um] =",F7.3," max[um] =",F7.3/
     4"Tht-Span: ","Ngrid =",I4," min[um] =",F7.3," max[um] =",F7.3/
     5"Phi-Span: ","Ngrid =",I4," min[um] =",F7.3," max[um] =",F7.3)
      WRITE(51,13) WAVEL,RADCOR,RADCOT
     
C Bohren: X*REFIM1, X*REFIM2, AND Y*REFIM2 SHOULD BE LESS THAN ABOUT 30
   15 FORMAT('!!! CAUTION !!! ',A8,' > 30: ',D10.4)
      IF(X*REFIM1.GE.30.0D0) THEN
       WRITE(51,15) 'X*REFIM1', X*REFIM1
      END IF
      IF(X*REFIM2.GE.30.0D0) THEN
       WRITE(51,15) 'X*REFIM2', X*REFIM2
      END IF
      IF(Y*REFIM2.GE.30.0D0) THEN
       WRITE(51,15) 'Y*REFIM2', Y*REFIM2
      END IF
     
     
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
      

   18 FORMAT('wavelength[nm],epsilon(ext),epsilon(abs),Qext,Qabs,',
     1'Qsca,Qback')
      WRITE(51,18)
   16 FORMAT('',F7.2,6E11.4)
      WRITE(51,16) WAVEL*1.0D3,QEXT,QABS,QSCA
C

  950 FORMAT('# ',A8,' field: ',A/'# r[um] theta[] phi[] ',A)
  951 FORMAT('# ',A8,' field: ',A/'# x[um] y[um] z[um] ',A)


      SFIELD(1)='Core    '
      SFIELD(2)='Coat    '
      SFIELD(3)='External'
      SFIELD(4)='All     '
      
      
      DO 80 IW=1,4
       OPEN(10+IW,FILE=FNESQ(IW),STATUS='UNKNOWN')
       OPEN(20+IW,FILE=FNUAB(IW),STATUS='UNKNOWN')
       OPEN(30+IW,FILE=FNUABC(IW),STATUS='UNKNOWN')
       WRITE(10+IW,950) SFIELD(IW),'EFSQ','EFSQ[dimensionless]'
       WRITE(20+IW,950) SFIELD(IW),'Uabs','UABS[F m-1 s-1]'
       WRITE(30+IW,951) SFIELD(IW),'Uabs','UABS[F m-1 s-1]'
   80 CONTINUE 
      
  701 FORMAT(3E13.5,E13.5)
C  706 FORMAT(3E13.5,6E13.5)
C  711 FORMAT(3E13.5,6E13.5,F6.3,3F6.1,I3)
C  712 FORMAT(3E13.5,6E13.5,F6.3,4F6.1,I3)
C  713 FORMAT(3E13.5,3E13.5,F7.2,6E13.5,E11.3)
      
      
C     create grid spacing for solution

C scan steps

      DO 9874 I=1,3
       IF(NGRID(I).GT.1) THEN
        XCSTEP(I)=(XCMAX(I)-XCMIN(I))/DBLE(NGRID(I)-1)
       ELSE
        XCSTEP(I)=0.0D0
       END IF
 9874 CONTINUE
 
C     this garbage makes a prescribed arc length for integration
C     choose length for smallest circle and all others are chosen from there
      RSTEP = (XPMAX(1)-XPMIN(1))/DBLE(RGRID-1)
      DELTAS1 = RSTEP*2*PI/DBLE(SGRID-1)
      CT = 0
      RAD = 0
      
      PRINT *, "RGRID = ", RGRID
      PRINT *, "RSTEP = ", RSTEP
      
      DO 100 I=1,RGRID
       RAD = DBLE(I-1)*RSTEP
       
       NPTSTH = NINT(PI*RAD/DELTAS1)
       THSTEP = PI/DBLE(NPTSTH-1)
       NPTSPH = NINT(2*PI*RAD/DELTAS1)
       PHSTEP = 2*PI/DBLE(NPTSPH-1)
       
       PRINT *, "***new outer loop***"
       PRINT *, "NPTSTH: ", NPTSTH
       PRINT *, "THSTEP: ", THSTEP
       PRINT *, "NPTSPH: ", NPTSPH
       PRINT *, "PHSTEP: ", PHSTEP
       PRINT *, ""

       XP(1) = RAD
C      Grab a point at the center so python can interpolate later:
       IF (RAD.EQ.0) THEN
        XP(2) = 0
        XP(3) = 0
        CALL FIELDVMW(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1              RADCOR,RADCOT,XP, IWHERE,EC,HC)
        EFSQ=ABS(EC(1))**2.0D0+ABS(EC(2))**2.0D0+ABS(EC(3))**2.0D0
        UABS=EPSVAC*OMEGA*REFRE1*REFIM1*EFSQ
        PRINT *, "at 0: ", UABS
        DO 129 IW=1,4
          IF(IW.EQ.4.OR.IW.EQ.IWHERE) THEN
           WRITE(10+IW,701) (XP(N),N=1,3),EFSQ
           WRITE(20+IW,701) (XP(N),N=1,3),UABS
          ELSE
           WRITE(10+IW,701) (XP(N),N=1,3),0.0D0
           WRITE(20+IW,701) (XP(N),N=1,3),0.0D0
          ENDIF
  129    CONTINUE
         CT = CT + 1
         PRINT *, CT, ": ", XP(1), ", ",XP(2), ", ",XP(3), IWHERE
       END IF
       
C      Now do the rest of the points      
       
       DO 110 J=1,NPTSTH
        XP(2)=XPMIN(2)+DBLE(J-1)*THSTEP
        DO 120 K=1,NPTSPH
         XP(3)=XPMIN(3)+DBLE(K-1)*PHSTEP
         
         
         CALL FIELDVMW(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1              RADCOR,RADCOT,XP, IWHERE,EC,HC)
         
         CT = CT + 1
         PRINT *, CT, ": ", XP(1), ", ",XP(2), ", ",XP(3), IWHERE
C
C vector electric field: snapshots [Re (t=0), Im (t=period/4)]
C electric field strength E*(E^*) [(V m-1)^2] = EFSQ * E0**2
C EFSQ: dimensionless
C
         EFSQ=ABS(EC(1))**2.0D0+ABS(EC(2))**2.0D0+ABS(EC(3))**2.0D0

C
C absorbed energy per unit volume and time Ua [W m-3] = UABS * E0**2
C eps0[F m-1], omega[s-1], E0[V m-1]; [F]=[C V-1], [W]=[C V s-1]
C UABS: [F m-1 s-1] = [(W m-3) (V m-1)^-2]
C Our formula: Uabs = EPSVAC*OMEGA*REFRE*REFIM*EFSQ (nonmagnetic)
C
         IF(IWHERE.EQ.1) THEN
C
          UABS=EPSVAC*OMEGA*REFRE1*REFIM1*EFSQ
C
         ELSE IF(IWHERE.EQ.2) THEN
C
          UABS=EPSVAC*OMEGA*REFRE2*REFIM2*EFSQ
C
         ELSE
C
          UABS=0.0D0
C
         END IF
         
C
C output data
C

C each domain (core, coat, or external) & all fields


         DO 130 IW=1,4
          IF(IW.EQ.4.OR.IW.EQ.IWHERE) THEN
           WRITE(10+IW,701) (XP(N),N=1,3),EFSQ
           WRITE(20+IW,701) (XP(N),N=1,3),UABS
          ELSE
           WRITE(10+IW,701) (XP(N),N=1,3),0.0D0
           WRITE(20+IW,701) (XP(N),N=1,3),0.0D0
          ENDIF
  130    CONTINUE


  120   CONTINUE
  110  CONTINUE
  100 CONTINUE   
      
      PRINT *, "Now for cartesian"
      CT = 0
      
      DO 101 I=1,NGRID(3)
       XC(3)=XCMIN(3)+DBLE(I-1)*XCSTEP(3)
       DO 111 J=1,NGRID(2)
        XC(2)=XCMIN(2)+DBLE(J-1)*XCSTEP(2)
        DO 121 K=1,NGRID(1)
         XC(1)=XCMIN(1)+DBLE(K-1)*XCSTEP(1)
         
C         PRINT *, XC(3),XC(2),XC(1)
C         PRINT *, SQRT(XP(1)*XP(1)+XP(2)*XP(2)+XP(3)*XP(3))
         CALL FIELD(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1              RADCOR,RADCOT,XC, IWHERE,EC,HC)
         
         CT = CT + 1
         PRINT *, CT, ": ", XC(1), ", ",XC(2), ", ",XC(3), IWHERE
C
C vector electric field: snapshots [Re (t=0), Im (t=period/4)]
C electric field strength E*(E^*) [(V m-1)^2] = EFSQ * E0**2
C EFSQ: dimensionless
C
         EFSQ=ABS(EC(1))**2.0D0+ABS(EC(2))**2.0D0+ABS(EC(3))**2.0D0

C
C absorbed energy per unit volume and time Ua [W m-3] = UABS * E0**2
C eps0[F m-1], omega[s-1], E0[V m-1]; [F]=[C V-1], [W]=[C V s-1]
C UABS: [F m-1 s-1] = [(W m-3) (V m-1)^-2]
C Our formula: Uabs = EPSVAC*OMEGA*REFRE*REFIM*EFSQ (nonmagnetic)
C
         IF(IWHERE.EQ.1) THEN
C
          UABS=EPSVAC*OMEGA*REFRE1*REFIM1*EFSQ
C
         ELSE IF(IWHERE.EQ.2) THEN
C
          UABS=EPSVAC*OMEGA*REFRE2*REFIM2*EFSQ
C
         ELSE
C
          UABS=0.0D0
C
         END IF
         
C
C output data
C

C each domain (core, coat, or external) & all fields


         DO 131 IW=1,4
          IF(IW.EQ.4.OR.IW.EQ.IWHERE) THEN
           !WRITE(10+IW,701) (XP(N),N=1,3),EFSQ
           WRITE(30+IW,701) (XC(N),N=1,3),UABS
          ELSE
           !WRITE(10+IW,701) (XP(N),N=1,3),0.0D0
           WRITE(30+IW,701) (XC(N),N=1,3),0.0D0
          ENDIF
  131    CONTINUE


  121   CONTINUE
  111  CONTINUE
  101 CONTINUE   
      
      PRINT *,"Qabs = ", QABS
      PRINT *,"Qsca = ", QSCA
      
      CALL IntegrateUABS(RADCOR,RADCOT,WAVEL,
     &                   REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     &                   EPSVAC,OMEGA,UABS_T)
      
      PRINT *, "UABS_T = ", UABS_T
      
      I0 = 0.5*SQRT(EPSVAC/MU)*1.
      AP = PI*RADCOT*RADCOT
      MYQABS = UABS_T/AP/I0*1.E-6
      
      PRINT *, ''
      PRINT *, 'Check full integral over sphere against Qabs'
      PRINT *, '************RESULTS*********************'
      PRINT *, "*  Fancy Qabs = ", QABS
      PRINT *, "*  My Qabs =    ", MYQABS
      PRINT *, '****************************************'
      !PRINT *, "Volume: ", 4./3.*pi*RADCOT**3.0D0
      
      
      NSHELLS = 10
      PRINT *, ''
      PRINT *, "Finally time to integrate some shells!"
      CALL IntegrateShells(RADCOR,RADCOT,WAVEL,
     &                   REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     &                   EPSVAC,OMEGA,NSHELLS)
      
      
      STOP
      END
      
C *******************************************************************
C  The goal of this whole mess:
C *********************************************************************
      subroutine IntegrateShells(RADCOR,RADCOT,WAVEL,
     &                   REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     &                   EPSVAC,OMEGA,NSHELLS)
     
      implicit none
     
      REAL*8 RADCOR,RADCOT,WAVEL,EPSVAC,OMEGA,RSTEP
      REAL*8 REFMED,REFRE1,REFIM1,REFRE2,REFIM2,RAD
      integer NSHELLS,I
      REAL*8 UABS(NSHELLS)
     
      PRINT *, 'This is where the magic happens'
      
      RSTEP = RADCOT/DBLE(NSHELLS-1)
      
      DO 100 I=1,NSHELLS
        RAD = DBLE(I-1)*RSTEP       
        CALL IntegrateUABSShell(RADCOR,RADCOT,WAVEL,
     &                         REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     &                         EPSVAC,OMEGA,UABS(I),RAD)
        
        
  100 CONTINUE
     
      PRINT *, "resulting field"
      DO 29 I=1,NSHELLS
        PRINT *, UABS(I)
   29 CONTINUE   
     
      return
      end
      
      
C *******************************************************************
C  This sets up and calls the cuba library to integrate over a shell
C *********************************************************************
  
      SUBROUTINE IntegrateUABSShell(RADCOR,RADCOT,WAVEL,
     &                         REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     &                         EPSVAC,OMEGA,UABS_T,RAD)

      implicit none

      integer ndim, ncomp, nvec, last, seed, mineval, maxeval
      cubareal epsrel, epsabs
      
      real*8 userdata(11),UABS_T
      
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

      REAL*8 RADCOR,RADCOT,WAVEL,EPSVAC,OMEGA
      REAL*8 REFMED,REFRE1,REFIM1,REFRE2,REFIM2,RAD
      
      external diffShell

      cubareal integral(ncomp), error(ncomp), prob(ncomp)
      integer verbose, neval, fail, nregions
      character*16 env

      integer c

      
      userdata(1) = WAVEL
      userdata(2) = REFMED
      userdata(3) = REFRE1
      userdata(4) = REFIM1
      userdata(5) = REFRE2
      userdata(6) = REFIM2
      userdata(7) = RADCOR
      userdata(8) = RADCOT
      userdata(9) = EPSVAC
      userdata(10) = OMEGA
      userdata(11) = RAD
   
      
 !     call getenv("CUBAVERBOSE", env)
 !     verbose = 2
 !     read(env, *, iostat=fail, end=999, err=999) verbose
 ! 999 continue

 !     print *, "-------------------- Vegas test --------------------"

      print *, "Integrating shell with radius ", RAD
 
      call Vegas(ndim, ncomp, diffShell, userdata, nvec, 
     &   epsrel, epsabs, verbose, seed, 
     &   mineval, maxeval, nstart, nincrease, nbatch, 
     &   gridno, statefile, spin, 
     &   neval, fail, integral, error, prob)

      !print *, "neval    =", neval
      !print *, "fail     =", fail
      print '(F20.12," +- ",F20.12,"   p = ",F8.3)', 
     &   (integral(c), error(c), prob(c), c = 1, ncomp)
      
      UABS_T = integral(1)

      RETURN
      END       

C *******************************************************************
C  This provides the integrand for an integral over a spherical shell
C *********************************************************************
        integer function diffShell(ndim, xx, ncomp, ff, userdata)
        implicit none
        integer ndim, ncomp, IWHERE
        cubareal xx(*), ff(*)
        real*8 userdata(11), pi
        REAL*8 RADCOR,RADCOT,WAVEL,XP(3),EPSVAC,OMEGA,RAD
        REAL*8 REFMED,REFRE1,REFIM1,REFRE2,REFIM2,EFSQ,UABS
        COMPLEX*16 EC(3),HC(3)
        parameter (pi = 3.14159265358979323846D0)
#define x xx(1)
#define y xx(2)
#define z xx(3)
#define f ff(1)
        
        WAVEL = userdata(1)
        REFMED = userdata(2)
        REFRE1 = userdata(3)
        REFIM1 = userdata(4)
        REFRE2 = userdata(5)
        REFIM2 = userdata(6)
        RADCOR = userdata(7)
        RADCOT = userdata(8)
        EPSVAC = userdata(9)
        OMEGA = userdata(10)
        RAD = userdata(11)
        
        !PRINT *, 'moopsy'
        !print *, 'WAVEL: ', WAVEL
        !print *, 'REFMED: ', REFMED
        !print *, 'REFRE1: ', REFRE1
        !print *, 'REFIM1: ', REFIM1
        !print *, 'REFRE2: ', REFRE2
        !print *, 'REFIM2: ', REFIM2
        !print *, 'RADCOR: ', RADCOR
        !print *, 'RADCOT: ', RADCOT
        !print *, 'RAD: ', RAD
                       
        XP(1) = RAD
        XP(2) = pi*y
        XP(3) = 2.0D0*pi*z
        
        !print *, XP(1),XP(2),XP(3)
        CALL FIELDVMW(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1                 RADCOR,RADCOT,XP,IWHERE,EC,HC)
        EFSQ=ABS(EC(1))**2.0D0+ABS(EC(2))**2.0D0+ABS(EC(3))**2.0D0
        IF(IWHERE.EQ.1) THEN
C
          UABS=EPSVAC*OMEGA*REFRE1*REFIM1*EFSQ
C
         ELSE IF(IWHERE.EQ.2) THEN
C
          UABS=EPSVAC*OMEGA*REFRE2*REFIM2*EFSQ
C
         ELSE
C
          UABS=0.0D0
C
         END IF

        f = UABS*2.0D0*pi*pi*x*x*sin(pi*y)


        diffShell = 0
        end

C *******************************************************************
C  This sets up and calls the cuba library to integrate over the sphere
C *********************************************************************
  
      SUBROUTINE IntegrateUABS(RADCOR,RADCOT,WAVEL,
     &                         REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     &                         EPSVAC,OMEGA,UABS_T)

      implicit none

      integer ndim, ncomp, nvec, last, seed, mineval, maxeval
      cubareal epsrel, epsabs
      
      real*8 userdata(10),UABS_T
      
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

      REAL*8 RADCOR,RADCOT,WAVEL,EPSVAC,OMEGA
      REAL*8 REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      
      external diffUABS

      cubareal integral(ncomp), error(ncomp), prob(ncomp)
      integer verbose, neval, fail, nregions
      character*16 env

      integer c

      
      userdata(1) = WAVEL
      userdata(2) = REFMED
      userdata(3) = REFRE1
      userdata(4) = REFIM1
      userdata(5) = REFRE2
      userdata(6) = REFIM2
      userdata(7) = RADCOR
      userdata(8) = RADCOT
      userdata(9) = EPSVAC
      userdata(10) = OMEGA
      
      PRINT *, 'moopsy'
      print *, 'WAVEL: ', WAVEL
      print *, 'REFMED: ', REFMED
      print *, 'REFRE1: ', REFRE1
      print *, 'REFIM1: ', REFIM1
      print *, 'REFRE2: ', REFRE2
      print *, 'REFIM2: ', REFIM2
      print *, 'RADCOR: ', RADCOR
      print *, 'RADCOT: ', RADCOT
      
      
      call getenv("CUBAVERBOSE", env)
      verbose = 2
      read(env, *, iostat=fail, end=999, err=999) verbose
  999 continue

      print *, "-------------------- Vegas test --------------------"

      call Vegas(ndim, ncomp, diffUABS, userdata, nvec, 
     &   epsrel, epsabs, verbose, seed, 
     &   mineval, maxeval, nstart, nincrease, nbatch, 
     &   gridno, statefile, spin, 
     &   neval, fail, integral, error, prob)

      print *, "neval    =", neval
      print *, "fail     =", fail
      print '(F20.12," +- ",F20.12,"   p = ",F8.3)', 
     &   (integral(c), error(c), prob(c), c = 1, ncomp)
      
      UABS_T = integral(1)

      RETURN
      END 
      
C************************************************************************

        integer function diffUABS(ndim, xx, ncomp, ff, userdata)
        implicit none
        integer ndim, ncomp, IWHERE
        cubareal xx(*), ff(*)
        real*8 userdata(10), pi
        REAL*8 RADCOR,RADCOT,WAVEL,XP(3),EPSVAC,OMEGA
        REAL*8 REFMED,REFRE1,REFIM1,REFRE2,REFIM2,EFSQ,UABS
        COMPLEX*16 EC(3),HC(3)
        parameter (pi = 3.14159265358979323846D0)
#define x xx(1)
#define y xx(2)
#define z xx(3)
#define f ff(1)
        
        WAVEL = userdata(1)
        REFMED = userdata(2)
        REFRE1 = userdata(3)
        REFIM1 = userdata(4)
        REFRE2 = userdata(5)
        REFIM2 = userdata(6)
        RADCOR = userdata(7)
        RADCOT = userdata(8)
        EPSVAC = userdata(9)
        OMEGA = userdata(10)
        
        !PRINT *, 'moopsy'
        !print *, 'WAVEL: ', WAVEL
        !print *, 'REFMED: ', REFMED
        !print *, 'REFRE1: ', REFRE1
        !print *, 'REFIM1: ', REFIM1
        !print *, 'REFRE2: ', REFRE2
        !print *, 'REFIM2: ', REFIM2
        !print *, 'RADCOR: ', RADCOR
        !print *, 'RADCOT: ', RADCOT
        
        XP(1) = RADCOT*x
        XP(2) = pi*y
        XP(3) = 2.0D0*pi*z
        
        !print *, XP(1),XP(2),XP(3)
        CALL FIELDVMW(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1                 RADCOR,RADCOT,XP,IWHERE,EC,HC)
        EFSQ=ABS(EC(1))**2.0D0+ABS(EC(2))**2.0D0+ABS(EC(3))**2.0D0
        IF(IWHERE.EQ.1) THEN
C
          UABS=EPSVAC*OMEGA*REFRE1*REFIM1*EFSQ
C
         ELSE IF(IWHERE.EQ.2) THEN
C
          UABS=EPSVAC*OMEGA*REFRE2*REFIM2*EFSQ
C
         ELSE
C
          UABS=0.0D0
C
         END IF

        f = UABS*RADCOT**3.0D0*2.0D0*pi*pi*x*x*sin(pi*y)


        diffUABS = 0
        end
      
      
  
      subroutine CubaTest()
      implicit none

      integer ndim, ncomp, nvec, last, seed, mineval, maxeval
      cubareal epsrel, epsabs!, userdata
      
      real*8 userdata
      
      parameter (ndim = 3)
      parameter (ncomp = 1)
      !parameter (userdata = 0)
      parameter (nvec = 1)
      parameter (epsrel = 1D-3)
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

      external integrand

      cubareal integral(ncomp), error(ncomp), prob(ncomp)
      integer verbose, neval, fail, nregions
      character*16 env

      integer c

      userdata = 29.0D0
      
      call getenv("CUBAVERBOSE", env)
      verbose = 2
      read(env, *, iostat=fail, end=999, err=999) verbose
  999 continue

      print *, "-------------------- Vegas test --------------------"

      call Vegas(ndim, ncomp, integrand, userdata, nvec, 
     &   epsrel, epsabs, verbose, seed, 
     &   mineval, maxeval, nstart, nincrease, nbatch, 
     &   gridno, statefile, spin, 
     &   neval, fail, integral, error, prob)

      print *, "neval    =", neval
      print *, "fail     =", fail
      print '(F20.12," +- ",F20.12,"   p = ",F8.3)', 
     &   (integral(c), error(c), prob(c), c = 1, ncomp)

      print *, " "
      print *, "-------------------- Suave test --------------------"

      call suave(ndim, ncomp, integrand, userdata, nvec,
     &   epsrel, epsabs, verbose + last, seed,
     &   mineval, maxeval, nnew, nmin, flatness,
     &   statefile, spin,
     &   nregions, neval, fail, integral, error, prob)

      print *, "nregions =", nregions
      print *, "neval    =", neval
      print *, "fail     =", fail
      print '(F20.12," +- ",F20.12,"   p = ",F8.3)',
     &   (integral(c), error(c), prob(c), c = 1, ncomp)
     
      end
      
      

      
      
C************************************************************************

        integer function integrand(ndim, xx, ncomp, ff, userdata)
        implicit none
        integer ndim, ncomp, IWHERE
        cubareal xx(*), ff(*)
        real*8 userdata
        !REAL*8 RADCOR,RADCOT,WAVEL,XP(3)
        !REAL*8 REFMED,REFRE1,REFIM1,REFRE2,REFIM2
        !COMPLEX*16 EC(3),HC(3)
#define x xx(1)
#define y xx(2)
#define z xx(3)
#define f ff(1)


        !CALL FIELDVMW(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
       !&            RADCOR,RADCOT,XP, IWHERE,EC,HC)
       ! EFSQ=ABS(EC(1))**2.0D0+ABS(EC(2))**2.0D0+ABS(EC(3))**2.0D0
       ! UABS=EPSVAC*OMEGA*REFRE1*REFIM1*EFSQ

        f = sin(2.*x)*cos(y)*exp(z)
        !PRINT *, x,y,z,f
        !print *, 'then...'
        !print *, f
        !print *, 'after'

        integrand = 0
        end
      
   
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     VMW rewrite of FIELD that receives and GIVES the E and H field in 
C     spherical coordinates
C     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE FIELDVMW(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1                 RADCOR,RADCOT,XP,IWHERE,EF,HF)
C
C     select field subroutines to calc E & H complex fields
C     (cartesian components)
C
      IMPLICIT NONE
      INTEGER IWHERE
      REAL*8 WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      REAL*8 RADCOR,RADCOT,R,THETA,PHI
      REAL*8 XP(3)
      COMPLEX*16 EF(3),HF(3)

#ifdef CHECK_TANGENTIAL_CONTINUITY
      INTEGER NCMAX
      PARAMETER (NCMAX=10)
      INTEGER IWHBEF,M,IBD
      INTEGER NC(3)
      REAL*8 RB,THETAB,PHIB,XPB(3)
      COMPLEX*16 EFB(3)
      SAVE IWHBEF,RB,THETAB,PHIB,NC,XPB,EFB
      DATA IWHBEF,NC,XPB/0,3*0,3*0.7772012D0/
#endif

      R = XP(1)
      THETA = XP(2)
      PHI = XP(3)
      
C
C avoid origin (R=0); division by rho (x) in subroutines
C use R=small but retain fake XP=YP=ZP=0 for vtk 3D plot
C
      IF(R.EQ.0.0D0) THEN
       R=MAX(RADCOR,RADCOT*1.0D-3)*1.0D-3
      END IF
C
      IF(R.LE.RADCOR) THEN
C core internal field
       IWHERE=1
       CALL FIELIN(WAVEL,REFRE1,REFIM1,R,THETA,PHI, EF,HF)
C
      ELSE IF(R.LE.RADCOT) THEN
C coating field
       IWHERE=2
       CALL FIELCT(WAVEL,REFRE2,REFIM2,
     1             R,THETA,PHI, EF,HF)
C
      ELSE
C external field = incident + scattered
       IWHERE=3
       CALL FIELEX(WAVEL,REFMED,R,THETA,PHI, EF,HF)
C
      END IF


C BH 59: check tangential continuity

#ifdef CHECK_TANGENTIAL_CONTINUITY

 1011 FORMAT('CHECK_TANGENTIAL_CONTINUITY: region[',
     1 I1,'] E(r),E(theta),E(phi) at (x y z=',3E15.8,
     2 ')(r theta phi=',E15.8,2F8.3,')->',
     3 3E15.8,' (Re,Im)',3('(',E15.8,',',E15.8,')'))

      IF((IWHERE.NE.IWHBEF).AND.(
     1 ((XP(1).EQ.XPB(1)).AND.(XP(2).EQ.XPB(2))).OR.
     2 ((XP(2).EQ.XPB(2)).AND.(XP(3).EQ.XPB(3))).OR.
     3 ((XP(3).EQ.XPB(3)).AND.(XP(1).EQ.XPB(1)))
     4 )) THEN

       IBD=IWHBEF+IWHERE-2
       NC(IBD)=NC(IBD)+1

       IF(NC(IBD).LT.NCMAX) THEN
        write(51,*) 'Crossing border:'
        write(51,1011) IWHBEF,(XPB(M),M=1,3),RB,THETAB,PHIB,
     1  (ABS(EFB(M)),M=1,3),(EFB(M),M=1,3)
        write(51,1011) IWHERE,(XP(M),M=1,3),R,THETA,PHI,
     1  (ABS(EF(M)),M=1,3),(EF(M),M=1,3)
       ELSE IF(NC(IBD).EQ.NCMAX) THEN
        write(51,*) 'etc ...'
       END IF

      END IF

      IWHBEF=IWHERE
      RB=R
      THETAB=THETA
      PHIB=PHI
      DO 200 M=1,3
       XPB(M)=XP(M)
       EFB(M)=EF(M)
  200 CONTINUE

#endif


      RETURN
      END   
      
      
      

      SUBROUTINE VMWTEST(IWHERE)
C     sanity check
      IMPLICIT NONE
      INTEGER IWHERE
C
      PRINT *, "hELLO, you asked for", IWHERE

      RETURN
      END
 
      
      SUBROUTINE RECMAX(IWHERE,XP,VALUE,VMAX)
C
C record max value and position for each domain
C
      IMPLICIT NONE
      INTEGER IWHERE,N
      REAL*8 VALUE
      REAL*8 XP(3),VMAX(3,4)
C
      IF(VALUE.GT.VMAX(IWHERE,4)) THEN
       VMAX(IWHERE,4)=VALUE
       DO 10 N=1,3
        VMAX(IWHERE,N)=XP(N)
   10  CONTINUE
      END IF
      RETURN
      END

      SUBROUTINE POL2CA(THETA,PHI,EP, EC)
C
C convert complex vector components based on spherical polar unit
C vectors (er, etheta, ephi) to those on cartesian (ex, ey, ez)
C
      IMPLICIT NONE
      REAL*8 THETA,PHI,ST,CT,SP,CP
      COMPLEX*16 EP(3),EC(3)
C
      ST=SIN(THETA)
      CT=COS(THETA)
      SP=SIN(PHI)
      CP=COS(PHI)
      EC(1)=ST*CP*EP(1)+CT*CP*EP(2)-SP*EP(3)
      EC(2)=ST*SP*EP(1)+CT*SP*EP(2)+CP*EP(3)
      EC(3)=CT*EP(1)-ST*EP(2)
      RETURN
      END

      SUBROUTINE ELLIPS(EC, EFMAJ,EFMIN,EPHI,ELLIPT,AZIM)
C
C compute electric field amplitudes [= major (a) & minor (b) axes 
C of vibration ellipse] from complex vector field 
C (Born & Wolf p.34-35 Eqs.62,66)
C max electric field amplitude = norm of EFMAJ
C PHASE [rad] = phase shift (epsilon) (-pi, pi]
C cf. -epsilon = eccentric angle between a and p; p = E(t=0) = Re(EC)
C EPHI [rad] = angle between p and a = atan(tan(epsilon)*b/a) (-pi, pi]
C ELLIPT = ellipticity norm(b)/norm(a)
C AZIM [rad] = azimuth angle between a and x-axis (tentative) [0, pi]
C (Bohren & Huffman p.46,50)
C
      IMPLICIT NONE
      INTEGER I
      REAL*8 P2,Q2,PQ,XT,YT,PHASE,EPHI,COSE,SINE,A,B,
     1       ELLIPT,AZIM
      REAL*8 P(3),Q(3),EFMAJ(3),EFMIN(3)
      COMPLEX*16 EC(3)
C
      DO 50 I=1,3
       P(I)=REAL(EC(I))
       Q(I)=AIMAG(EC(I))
   50 CONTINUE
      P2=P(1)*P(1)+P(2)*P(2)+P(3)*P(3)
      Q2=Q(1)*Q(1)+Q(2)*Q(2)+Q(3)*Q(3)
      PQ=P(1)*Q(1)+P(2)*Q(2)+P(3)*Q(3)
C     A2=(P2+Q2+SQRT((P2-Q2)*(P2-Q2)+4.0D0*PQ*PQ))/2.0D0
C     B2=(P2+Q2-SQRT((P2-Q2)*(P2-Q2)+4.0D0*PQ*PQ))/2.0D0

      XT=P2-Q2
      YT=2.0D0*PQ
      IF((XT.NE.0.0D0).OR.(YT.NE.0.0D0)) THEN
C atan2 result: (-pi, pi]
       PHASE=ATAN2(YT,XT)/2.0D0
      ELSE
       PHASE=0.0D0
      END IF

      COSE=COS(PHASE)
      SINE=SIN(PHASE)
      DO 60 I=1,3
       EFMAJ(I)=  P(I)*COSE +Q(I)*SINE
       EFMIN(I)=-(P(I)*SINE)+Q(I)*COSE
   60 CONTINUE
      A=SQRT(EFMAJ(1)*EFMAJ(1)+EFMAJ(2)*EFMAJ(2)+EFMAJ(3)*EFMAJ(3))
      B=SQRT(EFMIN(1)*EFMIN(1)+EFMIN(2)*EFMIN(2)+EFMIN(3)*EFMIN(3))

C     EPHI=ATAN2(B*TAN(PHASE),A)
      XT=A
      YT=B*TAN(PHASE)
      IF((XT.NE.0.0D0).OR.(YT.NE.0.0D0)) THEN
       EPHI=ATAN2(YT,XT)
      ELSE
       EPHI=0.0D0
      END IF

      IF(A.NE.0.0D0) THEN
       ELLIPT=B/A
C acos result: [0, pi]
       AZIM=ACOS(EFMAJ(1)/A)
      ELSE
       ELLIPT=0.0D0
       AZIM=0.0D0
      END IF

      RETURN
      END

      SUBROUTINE POYNTI(EC,HC, SF,RI,EHANG)
C
C compute optical intensity (irradiance) = magnitude of 
C time-averaged Poynting vector |<S>| = (1/2)|Re(Ec x Hc^*)|
C (Born & Wolf p.9-10)
C electric field E [V m-1] = EC * E0; EC[dimensionless]
C magnetic field H [A m-1] = HC * E0; HC[A V-1] = [(A m-1)(V m-1)-1]
C Poynting vector <S>[W m-2] = SF * (E0)^2; 
C SF[A V-1] = [(W m-2) (V m-1)-2]
C RI[A V-1]: optical intensity; norm of <S>
C EHANG: Re(E)-Re(H) angle [rad] [0, pi]
C
      IMPLICIT NONE
      INTEGER I
      REAL*8 RI,EHANG,ENORM,HNORM
      REAL*8 SF(3),E(3),H(3)
      COMPLEX*16 EC(3),HC(3),SC(3)
C
C complex Poynting vector
      SC(1)=(EC(2)*CONJG(HC(3))-EC(3)*CONJG(HC(2)))/2.0D0
      SC(2)=(EC(3)*CONJG(HC(1))-EC(1)*CONJG(HC(3)))/2.0D0
      SC(3)=(EC(1)*CONJG(HC(2))-EC(2)*CONJG(HC(1)))/2.0D0
C real Poynting vector
C avoid g77 buggy E-formatting: 0.12345-111 etc
      DO 50 I=1,3
       SF(I)=REAL(SC(I))
       IF(ABS(SF(I)).LT.1.0D-99) SF(I)=0.0D0
   50 CONTINUE
C optical intensity: norm of <S>
      RI=SQRT(SF(1)*SF(1)+SF(2)*SF(2)+SF(3)*SF(3))
C check angle between E & H
      DO 60 I=1,3
       E(I)=REAL(EC(I))
       H(I)=REAL(HC(I))
   60 CONTINUE
      ENORM=SQRT(E(1)*E(1)+E(2)*E(2)+E(3)*E(3))
      HNORM=SQRT(H(1)*H(1)+H(2)*H(2)+H(3)*H(3))
      EHANG=ACOS((E(1)*H(1)+E(2)*H(2)+E(3)*H(3))
     1      /(ENORM*HNORM))
      RETURN
      END















      SUBROUTINE HANDED(XC,SV, ANG,HAND)
C
C calc ellipsometric handedness (BH p.45) [tentative]
C XC complex vector field (EC or HC)
C SV Poynting vector
C ANG[deg] angle between Re(XC)xIm(XC) and SV
C (0 left-handed; 180 right-handed)
C HAND 1(right-handed) or -1(left-handed)
C if linearly polarized, return 0
C
      IMPLICIT NONE
      INTEGER N,HAND
      REAL*8 XYZ,XYN,ZPN,PI,ANG
      REAL*8 SV(3),XP(3),YP(3),ZP(3),XY(3)
      COMPLEX*16 XC(3)
C
C     PI=3.14159265D0
      PI=ACOS(-1.0D0)
C
      DO 10 N=1,3
       XP(N)=REAL(XC(N))
       YP(N)=AIMAG(XC(N))
       ZP(N)=SV(N)
   10 CONTINUE
C (X x Y).Z
      XY(1)=XP(2)*YP(3)-XP(3)*YP(2)
      XY(2)=XP(3)*YP(1)-XP(1)*YP(3)
      XY(3)=XP(1)*YP(2)-XP(2)*YP(1)
      XYN=SQRT(XY(1)*XY(1)+XY(2)*XY(2)+XY(3)*XY(3))
      ZPN=SQRT(ZP(1)*ZP(1)+ZP(2)*ZP(2)+ZP(3)*ZP(3))
      XYZ=XY(1)*ZP(1)+XY(2)*ZP(2)+XY(3)*ZP(3)
C
      IF(XYN.GT.0.0D0.AND.ZPN.GT.0.0D0) THEN
C acos -> [0, pi]
       ANG=ACOS(XYZ/(XYN*ZPN))*180.0D0/PI
C 'tranditional' convention
       IF(XYZ.GT.0.0D0) THEN
        HAND=-1
       ELSE IF(XYZ.EQ.0.0D0) THEN
        HAND=0
       ELSE
        HAND=1
       END IF
      ELSE
       ANG=0.0D0
       HAND=0
      END IF
C
      RETURN
      END


#ifdef CHECK_NUMERICAL_POYNTING

      SUBROUTINE CALDIV(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1                  RADCOR,RADCOT,XP,IWHERE,DH,SV, DIVS1,
     2                  DIVS2,DIVS3,SCMP)
C
C calc divS = dSx/dx + dSy/dy + dSz/dz
C Numerical differentiation BETA p.347
C y_k=y(x+kh)
C y' = (y_1 - y_0)/h (-> DIVS1)
C y' = (y_1 - y_-1)/2h (-> DIVS2)
C y' = (-y_2 + 8y_1 - 8y_-1 + y_-2)/12h
C y' = (y_3 - 9y_2 + 45y_1 - 45y_-1 + 9y_-2 - y_-3)/60h (-> DIVS3)

C numerical derivative is computed along orthonormal vectors 
C (DE) along S
C SCMP(N,M): Sx(x+kh) values on trial grids (for debug)
C
      IMPLICIT NONE
      INTEGER N,M,L,IDONE,IWHERE,IWHERE2
      REAL*8 WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1       RADCOR,RADCOT
      REAL*8 DH,DS1,DS2,DS3,DIVS1,DIVS2,DIVS3,DUM1,DUM2
      REAL*8 VNORM,SCZ
      REAL*8 XP(3),XPD(3),SV(3),SCD(3),SCP(3),SCM(3)
      REAL*8 DE(3,3),SCMP(3,8)
      COMPLEX*16 EC(3),HC(3)
C
C try to construct orthonormal vectors along S; this leads to 
C much better numerical results
C
      IDONE=0
C e1 = e(S)
      VNORM=SQRT(SV(1)*SV(1)+SV(2)*SV(2)+SV(3)*SV(3))
      IF(VNORM.EQ.0.0D0) GO TO 119
      DO 100 N=1,3
       DE(1,N)=SV(N)/VNORM
  100 CONTINUE
C e2 = e1 x ey
      DE(2,1)=-DE(1,3)
      DE(2,2)=0.0D0
      DE(2,3)=DE(1,1)
      VNORM=SQRT(DE(2,1)*DE(2,1)+DE(2,2)*DE(2,2)+DE(2,3)*DE(2,3))
      IF(VNORM.EQ.0.0D0) GO TO 119
      DO 105 N=1,3
       DE(2,N)=DE(2,N)/VNORM
  105 CONTINUE
C e3 = e1 x e2
      DE(3,1)=DE(1,2)*DE(2,3)-DE(1,3)*DE(2,2)
      DE(3,2)=DE(1,3)*DE(2,1)-DE(1,1)*DE(2,3)
      DE(3,3)=DE(1,1)*DE(2,2)-DE(1,2)*DE(2,1)
C
      IDONE=1
  119 CONTINUE
C if fails, just take ex, ey & ez
      IF(IDONE.EQ.0) THEN
       DO 120 N=1,3
        DO 125 M=1,3
         DE(N,M)=0.0D0
  125   CONTINUE
        DE(N,N)=1.0D0
  120  CONTINUE
      END IF
C
C numerical differentiation
C
      IDONE=0
C exclude origin (singular, fake in FIELD) in calc of divS
      IF((XP(1).EQ.0.0D0).AND.(XP(2).EQ.0.0D0).AND.
     1   (XP(3).EQ.0.0D0)) GO TO 148
      DIVS1=0.0D0
      DIVS2=0.0D0
      DIVS3=0.0D0

C loop: scan grids for each orthonormal vector e(N)

      DO 140 N=1,3
       SCZ=SV(1)*DE(N,1)+SV(2)*DE(N,2)+SV(3)*DE(N,3)
C
       DO 143 L=1,3
C XPD [um]
        DO 153 M=1,3
         XPD(M)=XP(M)+DBLE(L)*DH*DE(N,M)
  153   CONTINUE
        CALL FIELD(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1             RADCOR,RADCOT,XPD, IWHERE2,EC,HC)
        IF(IWHERE2.NE.IWHERE) GO TO 148
        CALL POYNTI(EC,HC, SCD,DUM1,DUM2)
C get component (Sx etc) in the coordination system (e1, e2, e3)
        SCP(L)=SCD(1)*DE(N,1)+SCD(2)*DE(N,2)+SCD(3)*DE(N,3)
  143  CONTINUE
C
       DO 146 L=1,3
        DO 156 M=1,3
         XPD(M)=XP(M)-DBLE(L)*DH*DE(N,M)
  156   CONTINUE
        CALL FIELD(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1             RADCOR,RADCOT,XPD, IWHERE2,EC,HC)
        IF(IWHERE2.NE.IWHERE) GO TO 148
        CALL POYNTI(EC,HC, SCD,DUM1,DUM2)
        SCM(L)=SCD(1)*DE(N,1)+SCD(2)*DE(N,2)+SCD(3)*DE(N,3)
  146  CONTINUE
C SC [A V-1] = [(W m-2) (V m-1)-2]
C DS, DIVS [(W m-3) (V m-1)-2]
       DS1=(SCP(1)-SCZ)/(DH*1.0D-6)
       DS2=(SCP(1)-SCM(1))/(2.0D0*DH*1.0D-6)
       DS3=(SCP(3)-9.0D0*SCP(2)+45.0D0*SCP(1)
     1     -SCM(3)+9.0D0*SCM(2)-45.0D0*SCM(1))
     2     /(60.0D0*DH*1.0D-6)

C record Sx(x+kh) etc values on trial grids (for debug)
       SCMP(N,1)=SCM(3)
       SCMP(N,2)=SCM(2)
       SCMP(N,3)=SCM(1)
       SCMP(N,4)=SCZ
       SCMP(N,5)=SCP(1)
       SCMP(N,6)=SCP(2)
       SCMP(N,7)=SCP(3)
       SCMP(N,8)=DS3
C
       DIVS1=DIVS1+DS1
       DIVS2=DIVS2+DS2
       DIVS3=DIVS3+DS3
  140 CONTINUE
C
      IDONE=1
C exclude boundary discontinuity (tentative)
  148 CONTINUE
      IF(IDONE.EQ.0) THEN
       DIVS1=0.0D0         
       DIVS2=0.0D0
       DIVS3=0.0D0
      END IF
      RETURN
      END

#endif /* CHECK_NUMERICAL_POYNTING */


      SUBROUTINE OPTCON(INIT,WAVEL,FNAME,WLFAC,
     1                  REFMED,REFRE1,REFIM1,REFRE2,REFIM2)
C
C get n,k values at specified wavelength
C
      IMPLICIT NONE
      INTEGER NDMAX
      PARAMETER (NDMAX=5000)
      INTEGER INIT,I,J
      REAL*8 WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      REAL*8 DATWL,DATRE,DATIM,RAT
      REAL*8 DOPT(3,3,NDMAX),FRE(3),FIM(3),WLFAC(3)
      INTEGER ND(3)
      CHARACTER*80 FNAME(3)
      CHARACTER*255 BUF
      SAVE DOPT,ND
C
      IF(INIT.EQ.1) THEN
       DO 550 I=1,3
        OPEN(1,FILE=FNAME(I),STATUS='OLD')
        J=0
  562   CONTINUE
        READ(1,'(A)',END=561) BUF
        IF(BUF(1:1).EQ.'#'.OR.BUF(1:1).EQ.';') GO TO 562
        IF(INDEX(BUF,'.').EQ.0) GO TO 562
        READ(BUF,*,ERR=562) DATWL,DATRE,DATIM
        IF(DATWL.EQ.0.0D0) GO TO 562
        J=J+1
        IF(J.GE.NDMAX) THEN
         WRITE(*,*) 'ndata > max for ',FNAME(I)
         STOP
        END IF
        DOPT(I,1,J)=DATWL*WLFAC(I)
        DOPT(I,2,J)=DATRE
        DOPT(I,3,J)=DATIM
        GO TO 562
  561   CONTINUE
        CLOSE(1)
        ND(I)=J
  570   FORMAT('Reading n,k(',I1,'): ndata = ',I4,' wavelength = ',
     *  E13.6,' - ',E13.6,' file = [',A,']')
        WRITE(*,570) I,ND(I),DOPT(I,1,1),DOPT(I,1,ND(I)),TRIM(FNAME(I))
        WRITE(51,570) I,ND(I),DOPT(I,1,1),DOPT(I,1,ND(I)),TRIM(FNAME(I))
  550  CONTINUE
      END IF
C
      DO 500 I=1,3
       IF(WAVEL.LT.DOPT(I,1,1).OR.WAVEL.GT.DOPT(I,1,ND(I))) THEN
        WRITE(*,*) 'wavelength out of range! ',WAVEL,' for data ',I
        STOP
       END IF
       DO 510 J=2,ND(I)
C wavelength is in the increasing order
        IF(WAVEL.LE.DOPT(I,1,J)) GO TO 511
  510  CONTINUE
  511  CONTINUE
C simple interpolation
       RAT=(WAVEL-DOPT(I,1,J-1))/(DOPT(I,1,J)-DOPT(I,1,J-1))
       FRE(I)=DOPT(I,2,J-1)+(DOPT(I,2,J)-DOPT(I,2,J-1))*RAT
       FIM(I)=DOPT(I,3,J-1)+(DOPT(I,3,J)-DOPT(I,3,J-1))*RAT
  500 CONTINUE
C
      REFMED=FRE(1)
C medium is assumed to be non-absorbing (k = 0)
      REFRE1=FRE(2)
      REFIM1=FIM(2)
      REFRE2=FRE(3)
      REFIM2=FIM(3)
C
      RETURN
      END

      SUBROUTINE BHCOAT (XG,YG,RFREL1G,RFREL2G,NSTOPF,
     1                   QEXT,QSCA,QBACK)
C       ************************************************
C       SUBROUTINE BHCOAT CALCULATES EFFICIENCIES FOR
C       EXTINCTION, TOTAL SCATTERING, AND BACKSCATTERING
C       FOR GIVEN SIZE PARAMETERS OF CORE AND COAT AND
C       RELATIVE REFRACTIVE INDICES
C       ALL BESSEL FUNCTIONS COMPUTED BY UPWARD RECURRENCE
C       ************************************************
C       caution: nonmagnetic meterial (mu = mu0) assumed p.183
C
C       input: X, Y (size parameters), RFREL1, RFREL2 (relative 
C       complex refractive indices), NSTOPF (max order for COEFF
C       necessary for field calc)
C       output: QEXT, QSCA, QBACK (efficiencies p.72), 
C       NORD,A,B,C,D,F,G,V,W (coeffs to be used in field calc 
C       p.93, 94, 182)
C
C       recurrence calc Bessel functions: j, D downward; y upward
C       p.128
C


      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=NMAX_COEFFS)
      REAL*8 UFTOL
      PARAMETER (UFTOL=1.0D-6)
      INTEGER NSTOP,N,IFLAG,NSTOPF,NORD
      REAL*8 XG,YG,QEXT,QSCA,QBACK,DEL,YSTOP,RNG,RCEQ1,RCEQ2
      COMPLEX*16 RFREL1G,RFREL2G,AN,BN,XBACK,CEQ1,CEQ2
      COMPLEX*16 AW,BW,CW,DW,FW,GW,VW,WW,TEMP
      COMPLEX*16 ANCW(NMAX),BNCW(NMAX)
      COMPLEX*16 PW(NMAX),QW(NMAX),RW(NMAX),SW(NMAX),TW(NMAX),UW(NMAX)
      COMMON /COEFF/A,B,C,D,F,G,V,W,NORD

      REAL*8 PSIBY,PSINY,ZERO,ONE,RN,X,Y
      COMPLEX*16 X1,X2,Y2,REFREL,YC,RFREL1,RFREL2,CI
      COMPLEX*16 CHIBX2,CHIBY2,CHINX2,CHINY2,DNX1,DNX2,DNY2
      COMPLEX*16 XXIBY,XXINY
      COMPLEX*16 PSINX1,PSINX2,PSINY2
      COMPLEX*16 CHIPX2,CHIPY2,ANCAP,BNCAP,DNBAR,GNBAR
      COMPLEX*16 CRACK,BRACK,AMESS1,AMESS2,AMESS3,AMESS4
      COMPLEX*16 P,Q,R,S,T,U
      COMPLEX*16 BJ(NMAX),BD(NMAX),BY(NMAX)
      COMPLEX*16 DX1(NMAX),DX2(NMAX),DY2(NMAX)
      COMPLEX*16 CHIX2(NMAX),CHIY2(NMAX)
      COMPLEX*16 PSIX1(NMAX),PSIX2(NMAX),PSIY2(NMAX),XXIY(NMAX)
      REAL*8 PSIY(NMAX),CHIY(NMAX)
      COMPLEX*16 A(NMAX),B(NMAX),C(NMAX),D(NMAX)
      COMPLEX*16 F(NMAX),G(NMAX),V(NMAX),W(NMAX)


#ifdef DEBUG_BESSEL
      INTEGER NP1,NP2
#endif

#ifdef DEBUG_BH183
      REAL*8 te2,te3,te4,te5
      COMPLEX*16 te1

      REAL*8 PSIPY
      COMPLEX*16 PSIPX1,PSIPX2,PSIPY2,PSIBX1,PSIBX2,PSIBY2,XXIPY

#endif

      COMMON /BJ/BJ
      COMMON /BD/BD
      COMMON /BY/BY


      ZERO=0.0D0
      ONE=1.0D0
      CI=DCMPLX(ZERO,ONE)


      X=XG
      Y=YG
      RFREL1=RFREL1G
      RFREL2=RFREL2G

C       *********************************************
C       DEL IS THE INNER SPHERE CONVERGENCE CRITERION
C       *********************************************
C disabled for now (see below)

      DEL=1.0D-8

      X1=RFREL1*X
      X2=RFREL2*X
      Y2=RFREL2*Y
      REFREL=RFREL2/RFREL1
C
C Wiscombe Criterion BH p.477,485; Barber-Hill p.196
C
      YSTOP=Y+4.0D0*Y**0.3333D0+2.0D0
      NSTOP=INT(YSTOP)
      NORD=MAX(NSTOP,NSTOPF)
      IF((NORD.GT.NMAX-1).OR.(NORD.LE.0)) THEN
       WRITE(*,*) 'NORD(',NORD,') > NMAX-1 (',NMAX-1,'): NSTOP=',
     1            NSTOP,' NSTOPF=',NSTOPF
       STOP
      END IF
C       ***********************************
C       SERIES TERMINATED AFTER NSTOP TERMS
C       ***********************************
      CALL BESJYD(X1,NORD)
      DO 140 N=1,NORD+1
       PSIX1(N)=X1*BJ(N)
  140 CONTINUE
      DO 100 N=1,NORD+1
       DX1(N)=BD(N)
  100 CONTINUE

      CALL BESJYD(X2,NORD)
      DO 150 N=1,NORD+1
       PSIX2(N)=X2*BJ(N)
       CHIX2(N)=-(X2*BY(N))
  150 CONTINUE
      DO 110 N=1,NORD+1
       DX2(N)=BD(N)
  110 CONTINUE

      CALL BESJYD(Y2,NORD)
      DO 160 N=1,NORD+1
       PSIY2(N)=Y2*BJ(N)
       CHIY2(N)=-(Y2*BY(N))
  160 CONTINUE
      DO 120 N=1,NORD+1
       DY2(N)=BD(N)
  120 CONTINUE
C
      YC=CMPLX(Y,ZERO)
      CALL BESJYD(YC,NORD)
      DO 130 N=1,NORD+1


       PSIY(N)=Y*REAL(BJ(N))
       CHIY(N)=-(Y*REAL(BY(N)))


       XXIY(N)=CMPLX(PSIY(N),-CHIY(N))
  130 CONTINUE
C
      IFLAG=0
C
      DO 200 N=1,NORD
C

       RN=DBLE(N)

       CHIBX2=CHIX2(N)
       CHINX2=CHIX2(N+1)
       CHIBY2=CHIY2(N)
       CHINY2=CHIY2(N+1)
       DNX1=DX1(N+1)
       DNX2=DX2(N+1)
       DNY2=DY2(N+1)

       PSIBY=PSIY(N)
       PSINY=PSIY(N+1)
       XXIBY=XXIY(N)
       XXINY=XXIY(N+1)
       PSINX1=PSIX1(N+1)
       PSINX2=PSIX2(N+1)
       PSINY2=PSIY2(N+1)

C BH 483
C      ANCAP = -(DNX2-REFREL*DNX1)*PSINX2/(REFREL*CHINX2*DNX1-CHIPX2)
C      BNCAP = -(REFREL*DNX2-DNX1)*PSINX2/(CHINX2*DNX1-REFREL*CHIPX2)
C BH 183 equivalently:
C      ANCAP = (RFREL2*PSINX2*PSIPX1-RFREL1*PSINX1*PSIPX2)
C     1       /(RFREL2*CHINX2*PSIPX1-RFREL1*CHIPX2*PSINX1)
C      BNCAP = (RFREL2*PSINX1*PSIPX2-RFREL1*PSINX2*PSIPX1)
C     1       /(RFREL2*CHIPX2*PSINX1-RFREL1*CHINX2*PSIPX1)
C BH 483 Wronskian relations: underflow prone!!!!!
C      PSINX2 = 1/(CHINX2*DNX2-CHIPX2)
C      PSINY2 = 1/(CHINY2*DNY2-CHIPY2)
C BH 484 recursive relations:
C      CHIPX2 = CHIBX2-RN*CHINX2/(RFREL2*X)
C      CHIPY2 = CHIBY2-RN*CHINY2/(RFREL2*Y)

       CHIPX2=CHIBX2-RN*CHINX2/X2
       CHIPY2=CHIBY2-RN*CHINY2/Y2

C       ANCAP=REFREL*DNX1-DNX2
C       ANCAP=ANCAP/(REFREL*DNX1*CHINX2-CHIPX2)
C       ANCAP=ANCAP/(CHINX2*DNX2-CHIPX2)
C       BRACK=ANCAP*(CHINY2*DNY2-CHIPY2)

       ANCAP=(REFREL*DNX1-DNX2)*PSINX2/(REFREL*CHINX2*DNX1-CHIPX2)
       BRACK=ANCAP/PSINY2

C       BNCAP=REFREL*DNX2-DNX1
C       BNCAP=BNCAP/(REFREL*CHIPX2-DNX1*CHINX2)
C       BNCAP=BNCAP/(CHINX2*DNX2-CHIPX2)
C       CRACK=BNCAP*(CHINY2*DNY2-CHIPY2)

       BNCAP=-(REFREL*DNX2-DNX1)*PSINX2/(CHINX2*DNX1-REFREL*CHIPX2)
       CRACK=BNCAP/PSINY2

       AMESS1=BRACK*CHIPY2
       AMESS2=BRACK*CHINY2
       AMESS3=CRACK*CHIPY2
       AMESS4=CRACK*CHINY2

C BH 483: underflow prone!!!!!
C      DNBAR = (DNY2-ANCAP*CHIPY2/PSINY2)/(1-ANCAP*CHINY2/PSINY2)
C      GNBAR = (DNY2-BNCAP*CHIPY2/PSINY2)/(1-BNCAP*CHINY2/PSINY2)

       DNBAR=DNY2-BRACK*CHIPY2
       DNBAR=DNBAR/(ONE-BRACK*CHINY2)
       GNBAR=DNY2-CRACK*CHIPY2
       GNBAR=GNBAR/(ONE-CRACK*CHINY2)


C Wronskian (BH p.483) underflow check

#ifdef CHECK_UNDERFLOW

       CALL UFCHECK2(CHIBX2,RN*CHINX2/X2,'CHIBX2-RN*CHINX2/X2')
       CALL UFCHECK2(CHIBY2,RN*CHINY2/Y2,'CHIBY2-RN*CHINY2/Y2')
       CALL UFCHECK2(CHINX2*DNX2,CHIPX2,'CHINX2*DNX2-CHIPX2')
       CALL UFCHECK2(CHINY2*DNY2,CHIPY2,'CHINY2*DNY2-CHIPY2')
       CALL EQCHECK(PSINX2,ONE/(CHINX2*DNX2-CHIPX2),
     1             'PSINX2=ONE/(CHINX2*DNX2-CHIPX2)')
       CALL EQCHECK(PSINY2,ONE/(CHINY2*DNY2-CHIPY2),
     1             'PSINY2=ONE/(CHINY2*DNY2-CHIPY2)')
C      CALL UFCHECK2(REFREL*DNX1,DNX2,'REFREL*DNX1-DNX2')
       CALL UFCHECK2(REFREL*CHINX2*DNX1,CHIPX2,
     1              'REFREL*CHINX2*DNX1-CHIPX2')
C      CALL UFCHECK2(REFREL*DNX2,DNX1,'REFREL*DNX2-DNX1')
       CALL UFCHECK2(CHINX2*DNX1,REFREL*CHIPX2,
     1              'CHINX2*DNX1-REFREL*CHIPX2')
       CALL UFCHECK2(DNY2,BRACK*CHIPY2,'DNY2-BRACK*CHIPY2')
       CALL UFCHECK2(CMPLX(ONE,ZERO),BRACK*CHINY2,'ONE-BRACK*CHINY2')
       CALL UFCHECK2(DNY2,CRACK*CHIPY2,'DNY2-CRACK*CHIPY2')
       CALL UFCHECK2(CMPLX(ONE,ZERO),CRACK*CHINY2,'ONE-CRACK*CHINY2')

#endif /* CHECK_UNDERFLOW */


C
C BH p.484 remedy for upward instability: unnecessary now?
C
C activated:
C      IF((IFLAG.EQ.1).OR.
C turned off:
       IF((IFLAG.EQ.1).AND.
     1    (ABS(AMESS1).LT.DEL*ABS(DNY2).AND.
     2     ABS(AMESS2).LT.DEL.AND.
     3     ABS(AMESS3).LT.DEL*ABS(DNY2).AND.
     4     ABS(AMESS4).LT.DEL)) THEN
        DNBAR=DNY2
        GNBAR=DNY2
        IFLAG=1
       END IF

C     AN=(DNBAR/RFREL2+RN/Y)*PSIY-PSI1Y
C     AN=AN/((DNBAR/RFREL2+RN/Y)*XXIY-XXI1Y)
C     BN=(RFREL2*GNBAR+RN/Y)*PSIY-PSI1Y
C     BN=BN/((RFREL2*GNBAR+RN/Y)*XXIY-XXI1Y)

C
C A(N) = a_n etc (n = 1, 2, ...)
C

C      A(N) = ((Y*DNBAR+RFREL2*RN)*PSINY-RFREL2*Y*PSIBY)/
C    1        ((Y*DNBAR+RFREL2*RN)*XXINY-RFREL2*Y*XXIBY)
C      B(N) = ((RFREL2*Y*GNBAR+RN)*PSINY-Y*PSIBY)/
C    1        ((RFREL2*Y*GNBAR+RN)*XXINY-Y*XXIBY)

C      A(N) = (PSINY*(DNBAR*Y+RFREL2*RN)-RFREL2*Y*PSIBY)
C    1       /(XXINY*(DNBAR*Y+RFREL2*RN)-RFREL2*Y*XXIBY)
C      B(N) = (PSINY*(GNBAR*RFREL2*Y+RN)-Y*PSIBY)
C    1       /(XXINY*(GNBAR*RFREL2*Y+RN)-Y*XXIBY)

       P = XXINY*(DNBAR*Y+RFREL2*RN)-RFREL2*Y*XXIBY
       Q = XXINY*(GNBAR*RFREL2*Y+RN)-Y*XXIBY

       A(N) = (PSINY*(DNBAR*Y+RFREL2*RN)-RFREL2*Y*PSIBY)/P
       B(N) = (PSINY*(GNBAR*RFREL2*Y+RN)-Y*PSIBY)/Q


#ifdef CHECK_UNDERFLOW

       CALL UFCHECK2(-DNBAR*Y,RFREL2*RN,'DNBAR*Y+RFREL2*RN')
       CALL UFCHECK2(PSINY*(DNBAR*Y+RFREL2*RN),RFREL2*Y*PSIBY,
     1              'PSINY*(DNBAR*Y+RFREL2*RN)-RFREL2*Y*PSIBY')
       CALL UFCHECK2(XXINY*(DNBAR*Y+RFREL2*RN),RFREL2*Y*XXIBY,
     1              'XXINY*(DNBAR*Y+RFREL2*RN)-RFREL2*Y*XXIBY')
       CALL UFCHECK2(-GNBAR*RFREL2*Y,CMPLX(RN,ZERO),'GNBAR*RFREL2*Y+RN')
       CALL UFCHECK2(PSINY*(GNBAR*RFREL2*Y+RN),CMPLX(Y*PSIBY,ZERO),
     1              'PSINY*(GNBAR*RFREL2*Y+RN)-Y*PSIBY')
       CALL UFCHECK2(XXINY*(GNBAR*RFREL2*Y+RN),Y*XXIBY,
     1              'XXINY*(GNBAR*RFREL2*Y+RN)-Y*XXIBY')

#endif /* CHECK_UNDERFLOW */


C      C(N) = Y*(RFREL2*X*DNX2*PSINX2-RFREL2*X*BNCAP*CHIBX2
C    1       +RN*BNCAP*CHINX2)*(PSINY*XXIBY-PSIBY*XXINY)
C    2       /(X*DNX1*PSINX1*(PSINY2-BNCAP*CHINY2)*(Y*XXIBY
C    3       -RFREL2*Y*GNBAR*XXINY-RN*XXINY))
C equivalently:
C      C(N) = RFREL1*Y*(PSINX2-BNCAP*CHINX2)*(PSINY*XXIBY-PSIBY*XXINY)
C    1       /(PSINX1*(PSINY2-BNCAP*CHINY2)*(Y*XXIBY-RFREL2*Y*GNBAR
C    2       *XXINY-RN*XXINY))
C      C(N) = RFREL2*Y*(PSINY*XXIBY-PSIBY*XXINY)/((CHINX2*DNX1
C    1       -REFREL*CHIPX2)*PSINX1*(PSINY2-BNCAP*CHINY2)*(Y*XXIBY
C    2       -RFREL2*Y*GNBAR*XXINY-RN*XXINY))
C
C underflow prone: PSINY2-BNCAP*CHINY2, PSINX2-BNCAP*CHINX2

C      C(N) = Y*(RFREL2*X*DNX2*PSINX2-RFREL2*X*BNCAP*CHIBX2
C    1       +RN*BNCAP*CHINX2)*(PSINY*XXIBY-PSIBY*XXINY)
C    2       /(X*DNX1*PSINX1*(PSINY2-BNCAP*CHINY2)
C    3       *(Y*XXIBY-RFREL2*Y*GNBAR*XXINY-RN*XXINY))

C BH 103,127 -> PSIBY*XXINY-PSINY*XXIBY = -i

       S = PSINY2-BNCAP*CHINY2
       U = RFREL2*X*DNX2*PSINX2-RFREL2*X*BNCAP*CHIBX2
     1    +RN*BNCAP*CHINX2

       C(N) = -CI*Y*U/(X*DNX1*PSINX1*S*Q)

       CEQ1 = RFREL1*Y*(PSINX2-BNCAP*CHINX2)*(PSINY*XXIBY
     1       -PSIBY*XXINY)/(PSINX1*(PSINY2-BNCAP*CHINY2)
     2       *(Y*XXIBY-RFREL2*Y*GNBAR*XXINY-RN*XXINY))

       CEQ2 = RFREL2*Y*(PSINY*XXIBY-PSIBY*XXINY)/((CHINX2*DNX1
     1       -REFREL*CHIPX2)*PSINX1*(PSINY2-BNCAP*CHINY2)
     2       *(Y*XXIBY-RFREL2*Y*GNBAR*XXINY-RN*XXINY))

C check equivalence of three C's as a round-off indicator
       RCEQ1=ABS((CEQ1-C(N))/C(N))
       RCEQ2=ABS((CEQ2-C(N))/C(N))
       IF((RCEQ1.GT.UFTOL).OR.(RCEQ2.GT.UFTOL)) THEN
        WRITE(51,*) 'c_n not equal! n,rdev(Ceq1),rdev(Ceq2)',
     1              N,RCEQ1,RCEQ2
       ENDIF


#ifdef CHECK_UNDERFLOW

       CALL UFCHECK2(PSINY*XXIBY,PSIBY*XXINY,
     1              'PSINY*XXIBY-PSIBY*XXINY')
       CALL UFCHECK2(PSINY2,BNCAP*CHINY2,'PSINY2-BNCAP*CHINY2')
       CALL UFCHECK3(RFREL2*X*DNX2*PSINX2,-RFREL2*X*BNCAP*CHIBX2,
     1               RN*BNCAP*CHINX2,
     2 'RFREL2*X*DNX2*PSINX2-RFREL2*X*BNCAP*CHIBX2+RN*BNCAP*CHINX2')
       CALL UFCHECK3(Y*XXIBY,-RFREL2*Y*GNBAR*XXINY,-RN*XXINY,
     1              'Y*XXIBY-RFREL2*Y*GNBAR*XXINY-RN*XXINY')
       CALL UFCHECK2(PSINX2,BNCAP*CHINX2,'PSINX2-BNCAP*CHINX2')
       CALL UFCHECK2(CHINX2*DNX1,REFREL*CHIPX2,
     1              'CHINX2*DNX1-REFREL*CHIPX2')
       CALL UFCHECK2(PSINY2,ANCAP*CHINY2,'PSINY2-ANCAP*CHINY2')
       CALL UFCHECK2(PSINX2,ANCAP*CHINX2,'PSINX2-ANCAP*CHINX2')

#endif /* CHECK_UNDERFLOW */


C      D(N) = RFREL2*Y*(PSINX2-ANCAP*CHINX2)*(PSINY*XXIBY-PSIBY*XXINY)
C    1       /(PSINX1*(PSINY2-ANCAP*CHINY2)*(RFREL2*Y*XXIBY
C    2       -Y*DNBAR*XXINY-RFREL2*RN*XXINY))
C      F(N) = RFREL2*Y*(PSINY*XXIBY-PSIBY*XXINY)/((PSINY2-BNCAP*CHINY2)
C    1       *(Y*XXIBY-RFREL2*Y*GNBAR*XXINY-RN*XXINY))
C      G(N) = RFREL2*Y*(PSINY*XXIBY-PSIBY*XXINY)/((PSINY2-ANCAP*CHINY2)
C    1       *(RFREL2*Y*XXIBY-Y*DNBAR*XXINY-RFREL2*RN*XXINY))
C      V(N) = RFREL2*Y*BNCAP*(PSINY*XXIBY-PSIBY*XXINY)/((PSINY2
C    1       -BNCAP*CHINY2)*(Y*XXIBY-RFREL2*Y*GNBAR*XXINY-RN*XXINY))
C      W(N) = RFREL2*Y*ANCAP*(PSINY*XXIBY-PSIBY*XXINY)/((PSINY2-ANCAP
C    1       *CHINY2)*(RFREL2*Y*XXIBY-Y*DNBAR*XXINY-RFREL2*RN*XXINY))

C      D(N) = RFREL2*Y*(PSINX2-ANCAP*CHINX2)*(PSINY*XXIBY
C    1       -PSIBY*XXINY)/(PSINX1*(PSINY2-ANCAP*CHINY2)
C    2       *(RFREL2*Y*XXIBY-Y*DNBAR*XXINY-RFREL2*RN*XXINY))

       R = PSINY2-ANCAP*CHINY2
       T = PSINX2-ANCAP*CHINX2

       D(N) = -CI*RFREL2*Y*T/(PSINX1*R*P)

C      F(N)  = RFREL2*Y*(PSINY*XXIBY-PSIBY*XXINY)/((PSINY2
C    1       -BNCAP*CHINY2)*(Y*XXIBY-RFREL2*Y*GNBAR*XXINY-RN*XXINY))

       F(N)  = -CI*RFREL2*Y/(S*Q)

C      G(N)  = RFREL2*Y*(PSINY*XXIBY-PSIBY*XXINY)/((PSINY2-ANCAP
C    1       *CHINY2)*(RFREL2*Y*XXIBY-Y*DNBAR*XXINY-RFREL2*RN*XXINY))

       G(N)  = -CI*RFREL2*Y/(R*P)

C      V(N) = RFREL2*Y*BNCAP*(PSINY*XXIBY-PSIBY*XXINY)/((PSINY2
C    1       -BNCAP*CHINY2)*(Y*XXIBY-RFREL2*Y*GNBAR*XXINY-RN*XXINY))
       V(N) = BNCAP*F(N)

C      W(N) = RFREL2*Y*ANCAP*(PSINY*XXIBY-PSIBY*XXINY)/((PSINY2-ANCAP
C    1       *CHINY2)*(RFREL2*Y*XXIBY-Y*DNBAR*XXINY-RFREL2*RN*XXINY))
       W(N) = ANCAP*G(N)

C no need to keep; just for write()
       ANCW(N)=ANCAP
       BNCW(N)=BNCAP
       PW(N)=P
       QW(N)=Q
       RW(N)=R
       SW(N)=S
       TW(N)=T
       UW(N)=U

C check BH 183 equations

#ifdef DEBUG_BH183

       if(n.le.10) then

C BH 127
        PSIBX1=PSIX1(N)
        PSIBX2=PSIX2(N)
        PSIBY2=PSIY2(N)

        PSIPX1=PSIBX1-RN*PSINX1/X1
        PSIPX2=PSIBX2-RN*PSINX2/X2
        PSIPY =PSIBY -RN*PSINY /Y
        PSIPY2=PSIBY2-RN*PSINY2/Y2
        XXIPY =XXIBY -RN*XXINY /Y
C BH 183
        write(51,*) ''
        te1=F(N)*RFREL1*PSINX2-V(N)*RFREL1*CHINX2-C(N)*RFREL2*PSINX1
        te2=abs(F(N)*RFREL1*PSINX2)
        te3=abs(V(N)*RFREL1*CHINX2)
        te4=abs(C(N)*RFREL2*PSINX1)
        write(51,*) 'BH 183-1',abs(te1),te2,te3,te4

        te1=W(N)*RFREL1*CHIPX2-G(N)*RFREL1*PSIPX2+D(N)*RFREL2*PSIPX1
        te2=abs(W(N)*RFREL1*CHIPX2)
        te3=abs(G(N)*RFREL1*PSIPX2)
        te4=abs(D(N)*RFREL2*PSIPX1)
        write(51,*) 'BH 183-2',abs(te1),te2,te3,te4

        te1=V(N)*CHIPX2-F(N)*PSIPX2+C(N)*PSIPX1
        te2=abs(V(N)*CHIPX2)
        te3=abs(F(N)*PSIPX2)
        te4=abs(C(N)*PSIPX1)
        write(51,*) 'BH 183-3',abs(te1),te2,te3,te4

        te1=G(N)*PSINX2-W(N)*CHINX2-D(N)*PSINX1
        write(51,*) 'BH 183-4',abs(te1)

        te1=RFREL2*PSIPY-A(N)*RFREL2*XXIPY-G(N)*PSIPY2+W(N)*CHIPY2
        te2=abs(RFREL2*PSIPY)
        te3=abs(A(N)*RFREL2*XXIPY)
        te4=abs(G(N)*PSIPY2)
        te5=abs(W(N)*CHIPY2)
        write(51,*) 'BH 183-5',abs(te1),te2,te3,te4,te5

        te1=RFREL2*B(N)*XXINY-RFREL2*PSINY+F(N)*PSINY2-V(N)*CHINY2
        write(51,*) 'BH 183-6',abs(te1)
        te1=PSINY-A(N)*XXINY-G(N)*PSINY2+W(N)*CHINY2
        write(51,*) 'BH 183-7',abs(te1)

        te1=B(N)*XXIPY-PSIPY+F(N)*PSIPY2-V(N)*CHIPY2
        te2=abs(B(N)*XXIPY)
        te3=abs(PSIPY)
        te4=abs(F(N)*PSIPY2)
        te5=abs(V(N)*CHIPY2)
        write(51,*) 'BH 183-8',abs(te1),te2,te3,te4,te5

        te1=PSIBY*XXINY-PSINY*XXIBY+CI
        te2=abs(PSIBY*XXINY)
        te3=abs(PSINY*XXIBY)
        write(51,*) 'PSIBY*XXINY-PSINY*XXIBY = -i',abs(te1),te2,te3

        te1=PSIBX1-RN*PSINX1/X1
        te2=abs(PSIBX1-RN*PSINX1/X1-DNX1*PSINX1)
        write(51,*) 'PSIPX1         ',abs(te1),te2
        te1=PSIBX2-RN*PSINX2/X2
        te2=abs(PSIBX2-RN*PSINX2/X2-DNX2*PSINX2)
        write(51,*) 'PSIPX2         ',abs(te1),te2
        te1=PSIBY2-RN*PSINY2/Y2
        te2=abs(PSIBY2-RN*PSINY2/Y2-DNY2*PSINY2)
        write(51,*) 'PSIPY2         ',abs(te1),te2

       end if

#endif /* DEBUG_BH183 */

C
  200 CONTINUE

  812 FORMAT(2E13.5,' ',7('(',E12.5,',',E12.5,')'))
      WRITE(51,*) 'BHCOAT pamareters: x,y, m1,m2,',
     1            'x1(=m1x),x2(=m2x),y2(=m2y),m(=m2/m1)'
      WRITE(51,812) XG,YG,RFREL1G,RFREL2G,RFREL1G*XG,RFREL2G*XG,
     1              RFREL2G*YG,RFREL2G/RFREL1G
  813 FORMAT(A3,I3,':',8('(',E12.5,',',E12.5,')'))
      WRITE(51,*) 'Scattering & field coeffs: n: a,b,c,d,f,g,v,w'
      DO 308 N=1,NORD
       AW=A(N)
       BW=B(N)
       CW=C(N)
       DW=D(N)
       FW=F(N)
       GW=G(N)
       VW=V(N)
       WW=W(N)
       WRITE(51,813) 'Sca',N,AW,BW,CW,DW,FW,GW,VW,WW
  308 CONTINUE

C 814 FORMAT(A3,I3,':',4('(',E12.5,',',E12.5,')'))
C     WRITE(51,*) 'Auxiliary: An,Bn'
C     DO 307 N=1,NORD
C      WRITE(51,814) 'Aux',N,ANCW(N),BNCW(N)
C 307 CONTINUE

  815 FORMAT(A3,I3,':',4E12.5,';',2E12.5,';',2E12.5)
      WRITE(51,*) 'Surface mode check: ',
     1'Abs(Denominator) P(a,d,g,w),Q(b,c,f,v),R(d,g,w),S(c,f,v);',
     2'Abs(Num/Denom) T/R(d),U/S(c); Abs(Auxiliary) An,Bn'
      DO 306 N=1,NORD
       WRITE(51,815) 'Aux',N,ABS(PW(N)),ABS(QW(N)),ABS(RW(N)),
     1 ABS(SW(N)),ABS(TW(N)/RW(N)),ABS(UW(N)/SW(N)),
     2 ABS(ANCW(N)),ABS(BNCW(N))
  306 CONTINUE

C Quick hack to check NaN
      TEMP=A(NORD)+F(NORD)+W(NORD)
      IF(ISNAN(REAL(TEMP))) THEN
       WRITE(*,*) 'NaN detected in scattering coeffs'
       STOP
      END IF

#ifdef DEBUG_BESSEL
      NP1=2
      NP2=NORD+1
      WRITE(51,*) 'BHCOAT Spherical Bessels for n = ',NP1-1,'&',NP2-1
  821 FORMAT(A19,2E13.5,' & ',2E13.5)
  822 FORMAT(A19,3('(',E12.5,',',E12.5,')'),' & ',
     1           3('(',E12.5,',',E12.5,')'))

      WRITE(51,821) 'PSIY,CHIY:',PSIY(NP1),CHIY(NP1),
     1                           PSIY(NP2),CHIY(NP2)
      WRITE(51,822) 'BJ,BD,BY:',BJ(NP1),BD(NP1),BY(NP1),
     1                          BJ(NP2),BD(NP2),BY(NP2)
      WRITE(51,822) 'PSIX1,PSIX2,PSIY2:',
     1               PSIX1(NP1),PSIX2(NP1),PSIY2(NP1),
     2               PSIX1(NP2),PSIX2(NP2),PSIY2(NP2)
      WRITE(51,822) 'CHIX2,CHIY2,XXIY:',
     1               CHIX2(NP1),CHIY2(NP1),XXIY(NP1),
     2               CHIX2(NP2),CHIY2(NP2),XXIY(NP2)
      WRITE(51,822) 'DX1,DX2,DY2:',DX1(NP1),DX2(NP1),DY2(NP1),
     1                             DX1(NP2),DX2(NP2),DY2(NP2)

#endif /* DEBUG_BESSEL */

C
C check dependence of QSCA etc on order NORD
C
      IF(NORD.GT.NSTOP) THEN
       QSCA=ZERO
       QEXT=ZERO
       XBACK=CMPLX(ZERO,ZERO)
       DO 300 N=1,NORD
        RNG=DBLE(N)
        AN=A(N)
        BN=B(N)
        QSCA=QSCA+(2.0D0*RNG+1.0D0)*(ABS(AN)*ABS(AN)+ABS(BN)*ABS(BN))
        QEXT=QEXT+(2.0D0*RNG+1.0D0)*(REAL(AN)+REAL(BN))
C        XBACK=XBACK+(2.0D0*RNG+1.0D0)*(-1.0D0)**RNG*(AN-BN)
C <arprec> bug: (-1.0D0)**RN -> MPLOG error!
        XBACK=XBACK+(2.0D0*RNG+1.0D0)*DBLE((-1)**N)*(AN-BN)
  300  CONTINUE
       QSCA=(2.0D0/(Y*Y))*QSCA
       QEXT=(2.0D0/(Y*Y))*QEXT
C      QBACK=XBACK*CONJG(XBACK)
       QBACK=REAL(XBACK*CONJG(XBACK))
       QBACK=(1.0D0/(Y*Y))*QBACK

 5152  FORMAT(/"BHCOAT: Check consistency for different sum orders")
       WRITE(51,5152)
       WRITE(51,*) '[NORD ] QSCA, QEXT, QBACK: ',NORD,QSCA,QEXT,QBACK
      END IF
C
C final result to be returned
C
      QSCA=ZERO
      QEXT=ZERO
      XBACK=CMPLX(ZERO,ZERO)
      DO 310 N=1,NSTOP
       RNG=DBLE(N)
       AN=A(N)
       BN=B(N)
       QSCA=QSCA+(2.0D0*RNG+1.0D0)*(ABS(AN)*ABS(AN)+ABS(BN)*ABS(BN))
       QEXT=QEXT+(2.0D0*RNG+1.0D0)*(REAL(AN)+REAL(BN))
C       XBACK=XBACK+(2.0D0*RNG+1.0D0)*(-1.0D0)**RNG*(AN-BN)
C <arprec> bug: (-1.0D0)**RN -> MPLOG error!
       XBACK=XBACK+(2.0D0*RNG+1.0D0)*DBLE((-1)**N)*(AN-BN)
  310 CONTINUE
      QSCA=(2.0D0/(Y*Y))*QSCA
      QEXT=(2.0D0/(Y*Y))*QEXT
C     QBACK=XBACK*CONJG(XBACK)
      QBACK=REAL(XBACK*CONJG(XBACK))
      QBACK=(1.0D0/(Y*Y))*QBACK

      WRITE(51,*) '[NSTOP] QSCA, QEXT, QBACK: ',NSTOP,QSCA,QEXT,QBACK

 5153 FORMAT("cf. x->Inf Qext->2 (BH107); Qback->R(0deg) (BH123)"/)
      WRITE(51,5153)
C
      RETURN
      END

      SUBROUTINE FIELD(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1                 RADCOR,RADCOT,XP, IWHERE,EC,HC)
C
C     select field subroutines to calc E & H complex fields
C     (cartesian components)
C     VMW: only takes a single point in cartesian space and returns local field
C           (also in cartesian space)
      IMPLICIT NONE
      INTEGER IWHERE
      REAL*8 WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      REAL*8 RADCOR,RADCOT,R,THETA,PHI
      REAL*8 XP(3)
      COMPLEX*16 EF(3),HF(3),EC(3),HC(3)

#ifdef CHECK_TANGENTIAL_CONTINUITY
      INTEGER NCMAX
      PARAMETER (NCMAX=10)
      INTEGER IWHBEF,M,IBD
      INTEGER NC(3)
      REAL*8 RB,THETAB,PHIB,XPB(3)
      COMPLEX*16 EFB(3)
      SAVE IWHBEF,RB,THETAB,PHIB,NC,XPB,EFB
      DATA IWHBEF,NC,XPB/0,3*0,3*0.7772012D0/
#endif

C     VMW: translate cartesian point to spherical
      R=SQRT(XP(1)*XP(1)+XP(2)*XP(2)+XP(3)*XP(3))
      IF(R.GT.0.0D0) THEN
       THETA=ACOS(XP(3)/R)
      ELSE
       THETA=0.0D0
      END IF
      IF((XP(1).NE.0.0D0).OR.(XP(2).NE.0.0D0)) THEN
       PHI=ATAN2(XP(2),XP(1))
      ELSE
       PHI=0.0D0
      END IF
C
C avoid origin (R=0); division by rho (x) in subroutines
C use R=small but retain fake XP=YP=ZP=0 for vtk 3D plot
C
      IF(R.EQ.0.0D0) THEN
       R=MAX(RADCOR,RADCOT*1.0D-3)*1.0D-3
      END IF
C     VMW: 
C
      IF(R.LE.RADCOR) THEN
C core internal field
       IWHERE=1
       CALL FIELIN(WAVEL,REFRE1,REFIM1,R,THETA,PHI, EF,HF)
C
      ELSE IF(R.LE.RADCOT) THEN
C coating field
       IWHERE=2
       CALL FIELCT(WAVEL,REFRE2,REFIM2,
     1             R,THETA,PHI, EF,HF)
C
      ELSE
C external field = incident + scattered
       IWHERE=3
       CALL FIELEX(WAVEL,REFMED,R,THETA,PHI, EF,HF)
C
      END IF


C BH 59: check tangential continuity

#ifdef CHECK_TANGENTIAL_CONTINUITY

 1011 FORMAT('CHECK_TANGENTIAL_CONTINUITY: region[',
     1 I1,'] E(r),E(theta),E(phi) at (x y z=',3E15.8,
     2 ')(r theta phi=',E15.8,2F8.3,')->',
     3 3E15.8,' (Re,Im)',3('(',E15.8,',',E15.8,')'))

      IF((IWHERE.NE.IWHBEF).AND.(
     1 ((XP(1).EQ.XPB(1)).AND.(XP(2).EQ.XPB(2))).OR.
     2 ((XP(2).EQ.XPB(2)).AND.(XP(3).EQ.XPB(3))).OR.
     3 ((XP(3).EQ.XPB(3)).AND.(XP(1).EQ.XPB(1)))
     4 )) THEN

       IBD=IWHBEF+IWHERE-2
       NC(IBD)=NC(IBD)+1

       IF(NC(IBD).LT.NCMAX) THEN
        write(51,*) 'Crossing border:'
        write(51,1011) IWHBEF,(XPB(M),M=1,3),RB,THETAB,PHIB,
     1  (ABS(EFB(M)),M=1,3),(EFB(M),M=1,3)
        write(51,1011) IWHERE,(XP(M),M=1,3),R,THETA,PHI,
     1  (ABS(EF(M)),M=1,3),(EF(M),M=1,3)
       ELSE IF(NC(IBD).EQ.NCMAX) THEN
        write(51,*) 'etc ...'
       END IF

      END IF

      IWHBEF=IWHERE
      RB=R
      THETAB=THETA
      PHIB=PHI
      DO 200 M=1,3
       XPB(M)=XP(M)
       EFB(M)=EF(M)
  200 CONTINUE

#endif


C convert to cartesian components
      CALL POL2CA(THETA,PHI,EF, EC)
      CALL POL2CA(THETA,PHI,HF, HC)
      RETURN
      END

      SUBROUTINE FIELIN(WAVEL,REFRE1,REFIM1,R,THETAG,PHIG, EFG,HFG)
C
C     core internal field
C     BH p.93 (4.40), 95 (4.50)
C     input units: WAVEL[um], R[um]
C     vector components use orthonormal basis vectors associated with 
C     spherical polar coordinate system (r, theta, phi), i.e.,
C     complex E-field Ec = EF(1) * er + EF(2) * etheta + EF(3) * ephi
C


      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=NMAX_COEFFS)
      INTEGER NSTOPF,NORD,I,N,N1
      REAL*8 WAVEL,REFRE1,REFIM1,R,THETAG,PHIG
      REAL*8 CC,OMEGA,MU1
      COMPLEX*16 EFG(3),HFG(3)
      COMMON /COEFF/A,B,C,D,F,G,V,W,NORD


      REAL*8 ZERO,ONE,TWO,PI,RN
      COMPLEX*16 CI,CREF1,WVK1,RHO1,BJN1,ENCAP,PSIP
      COMPLEX*16 HFFACT
      COMPLEX*16 BJ(NMAX)
      COMPLEX*16 A(NMAX),B(NMAX),C(NMAX),D(NMAX)
      COMPLEX*16 F(NMAX),G(NMAX),V(NMAX),W(NMAX)
      COMPLEX*16 VM1O1N(3),VM1E1N(3),VN1O1N(3),VN1E1N(3)
      COMPLEX*16 EF(3),HF(3)



      REAL*8 THETA,PHI
      REAL*8 PIN(NMAX),TAU(NMAX)

      COMMON /BJ/BJ
      COMMON /LEG/PIN,TAU


      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0


C     PI=3.14159265D0
      PI=ACOS(-ONE)

      THETA=THETAG
      PHI=PHIG

C
C Wiscombe Criterion BH p.477,485; Barber-Hill p.196
C
      CREF1=DCMPLX(REFRE1,REFIM1)
      WVK1=TWO*PI*CREF1/WAVEL
      RHO1=WVK1*R
      NSTOPF=AINT(ABS(RHO1)+4.05D0*ABS(RHO1)**.3333D0+2.0D0)
      IF((NSTOPF.GT.NMAX-1).OR.(NSTOPF.LE.0)) THEN
       WRITE(*,*) 'FIELIN: NSTOPF(',NSTOPF,') > NMAX-1 (',
     * NMAX-1,')'
       STOP
      ELSE IF(NSTOPF.GT.NORD) THEN
       WRITE(*,*) 'FIELIN: NSTOPF(',NSTOPF,') > NORD (',
     * NORD,') at R=',R
       STOP
      END IF
C
      CALL BESJYD(RHO1,NSTOPF)
      CALL LEGEND(THETA,NSTOPF)
C
      CI=CMPLX(ZERO,ONE)
      DO 10 I=1,3
       EF(I)=CMPLX(ZERO,ZERO)
       HF(I)=CMPLX(ZERO,ZERO)
   10 CONTINUE
C
      DO 20 N=1,NSTOPF
       N1=N+1
C      RN=DBLE(N)
       RN=N*ONE
C
       PSIP=RHO1*BJ(N1-1)-RN*BJ(N1)
C
       BJN1=BJ(N1)
       VM1O1N(1)= CMPLX(ZERO,ZERO)
       VM1O1N(2)= COS(PHI)*PIN(N1)*BJN1
       VM1O1N(3)=-(SIN(PHI)*TAU(N1)*BJN1)
       VM1E1N(1)= CMPLX(ZERO,ZERO)
       VM1E1N(2)=-(SIN(PHI)*PIN(N1)*BJN1)
       VM1E1N(3)=-(COS(PHI)*TAU(N1)*BJN1)
       VN1O1N(1)= SIN(PHI)*RN*(RN+ONE)*SIN(THETA)*PIN(N1)*BJN1
     1           /RHO1
       VN1O1N(2)= SIN(PHI)*TAU(N1)*PSIP/RHO1
       VN1O1N(3)= COS(PHI)*PIN(N1)*PSIP/RHO1
       VN1E1N(1)= COS(PHI)*RN*(RN+ONE)*SIN(THETA)*PIN(N1)*BJN1
     1           /RHO1
       VN1E1N(2)= COS(PHI)*TAU(N1)*PSIP/RHO1
       VN1E1N(3)=-(SIN(PHI)*PIN(N1)*PSIP/RHO1)
C
       ENCAP=CI**RN*(TWO*RN+ONE)/(RN*(RN+ONE))
       DO 30 I=1,3
        EF(I)=EF(I)+ENCAP*(C(N)*VM1O1N(I)-CI*D(N)*VN1E1N(I))
        HF(I)=HF(I)+ENCAP*(D(N)*VM1E1N(I)+CI*C(N)*VN1O1N(I))
   30  CONTINUE
   20 CONTINUE
C
C electric field E [V m-1] = EF * E0
C EF: dimensionless
C magnetic field H [A m-1] = HF * E0
C HF [A V-1] = [(A m-1) (V m-1)^-1]
C
C ASSUME NON-MAGNETIC (MU=MU0=const) [N A-2]
      MU1=4.0D0*PI*1.0D-7
C light speed [m s-1]
      CC=2.99792458D8
C angular frequency [s-1]
      OMEGA=2.0D0*PI*CC/(WAVEL*1.0D-6)
C k[um-1]; [N m s-1 A-1] = [V]
C H factor: [A V-1]
      HFFACT=-((WVK1*1.0D6)/(OMEGA*MU1))
      DO 40 I=1,3
       HF(I)=HFFACT*HF(I)
   40 CONTINUE
C
      DO 50 I=1,3
       EFG(I)=EF(I)
       HFG(I)=HF(I)
   50 CONTINUE

      RETURN
      END

      SUBROUTINE FIELEX(WAVEL,REFMED,R,THETAG,PHIG, EFG,HFG)
C
C     external scattering field = incident + scattered
C     BH p.92 (4.37), 94 (4.45), 95 (4.50)
C     assume: medium is non-absorbing; refim = 0; Uabs = 0
C


      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=NMAX_COEFFS)
      INTEGER NSTOPF,NORD,I,N,N1
      REAL*8 WAVEL,REFMED,R,THETAG,PHIG
      REAL*8 CC,OMEGA,MU
      COMPLEX*16 EFG(3),HFG(3)
      COMMON /COEFF/A,B,C,D,F,G,V,W,NORD


      REAL*8 ZERO,ONE,TWO,PI,RN
      COMPLEX*16 CI,CREF,WVK,RHO,ENCAP,XXIP,ZN
      COMPLEX*16 HFFACT
      COMPLEX*16 BJ(NMAX),BY(NMAX)
      COMPLEX*16 A(NMAX),B(NMAX),C(NMAX),D(NMAX)
      COMPLEX*16 F(NMAX),G(NMAX),V(NMAX),W(NMAX)
      COMPLEX*16 VM3O1N(3),VM3E1N(3),VN3O1N(3),VN3E1N(3)
      COMPLEX*16 EF(3),HF(3),EI(3),HI(3),ES(3),HS(3)


      REAL*8 THETA,PHI
      REAL*8 PIN(NMAX),TAU(NMAX)
      COMPLEX*16 HFFACTA,EIFAC,HIFAC


      COMMON /BJ/BJ
      COMMON /BY/BY
      COMMON /LEG/PIN,TAU


      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0

C     PI=3.14159265D0
      PI=ACOS(-ONE)

      THETA=THETAG
      PHI=PHIG

C
C Wiscombe Criterion BH p.477,485; Barber-Hill p.196
C
      CREF=DCMPLX(REFMED,0.0D0)
      WVK=TWO*PI*CREF/WAVEL
      RHO=WVK*R
      NSTOPF=AINT(ABS(RHO)+4.05D0*ABS(RHO)**.3333D0+2.0D0)
      IF((NSTOPF.GT.NMAX-1).OR.(NSTOPF.LE.0)) THEN
       WRITE(*,*) 'FIELEX: NSTOPF(',NSTOPF,') > NMAX-1 (',
     * NMAX-1,')'
       STOP
      ELSE IF(NSTOPF.GT.NORD) THEN
       WRITE(*,*) 'FIELEX: NSTOPF(',NSTOPF,') > NORD (',
     * NORD,') at R=',R
       STOP
      END IF
C
      CALL BESJYD(RHO,NSTOPF)
      CALL LEGEND(THETA,NSTOPF)
C
      CI=CMPLX(ZERO,ONE)
      DO 10 I=1,3
       ES(I)=CMPLX(ZERO,ZERO)
       HS(I)=CMPLX(ZERO,ZERO)
   10 CONTINUE
C
      DO 20 N=1,NSTOPF
       N1=N+1
C      RN=DBLE(N)
       RN=N*ONE
C
       ZN=BJ(N1)+CI*BY(N1)
       XXIP=RHO*(BJ(N1-1)+CI*BY(N1-1))-RN*ZN
C
       VM3O1N(1)= CMPLX(ZERO,ZERO)
       VM3O1N(2)= COS(PHI)*PIN(N1)*ZN
       VM3O1N(3)=-(SIN(PHI)*TAU(N1)*ZN)
       VM3E1N(1)= CMPLX(ZERO,ZERO)
       VM3E1N(2)=-(SIN(PHI)*PIN(N1)*ZN)
       VM3E1N(3)=-(COS(PHI)*TAU(N1)*ZN)
       VN3O1N(1)= SIN(PHI)*RN*(RN+ONE)*SIN(THETA)*PIN(N1)*ZN
     1           /RHO
       VN3O1N(2)= SIN(PHI)*TAU(N1)*XXIP/RHO
       VN3O1N(3)= COS(PHI)*PIN(N1)*XXIP/RHO
       VN3E1N(1)= COS(PHI)*RN*(RN+ONE)*SIN(THETA)*PIN(N1)*ZN
     1           /RHO
       VN3E1N(2)= COS(PHI)*TAU(N1)*XXIP/RHO
       VN3E1N(3)=-(SIN(PHI)*PIN(N1)*XXIP/RHO)

C
C scattered field: BH p.94 (4.45)
       ENCAP=CI**RN*(TWO*RN+ONE)/(RN*(RN+ONE))
       DO 30 I=1,3
        ES(I)=ES(I)+ENCAP*(CI*A(N)*VN3E1N(I)-B(N)*VM3O1N(I))
        HS(I)=HS(I)+ENCAP*(CI*B(N)*VN3O1N(I)+A(N)*VM3E1N(I))
   30  CONTINUE
   20 CONTINUE
C
C incident E field: BH p.89 (4.21); cf. p.92 (4.37), p.93 (4.38)
C basis unit vectors = er, etheta, ephi
C
C     EIFAC=CEXP(CI*WVK*R*COS(THETA))
      EIFAC=EXP(CI*WVK*R*COS(THETA))
      EI(1)= EIFAC*SIN(THETA)*COS(PHI)
      EI(2)= EIFAC*COS(THETA)*COS(PHI)
      EI(3)=-(EIFAC*SIN(PHI))
C
C electric field E [V m-1] = EF * E0
C
      DO 35 I=1,3
       EF(I)=EI(I)+ES(I)
   35 CONTINUE
C
C magnetic field
C
C ASSUME NON-MAGNETIC (MU=MU0=const) [N A-2]
      MU=4.0D0*PI*1.0D-7
C light speed [m s-1]
      CC=2.99792458D8
      OMEGA=2.0D0*PI*CC/(WAVEL*1.0D-6)
      HFFACT=(WVK*1.0D6)/(OMEGA*MU)
      DO 40 I=1,3
       HS(I)=HFFACT*HS(I)
   40 CONTINUE
C
C incident H field: BH p.26 (2.43), p.89 (4.21)
C
      HFFACTA=HFFACT
      HIFAC=EIFAC*HFFACTA
      HI(1)=HIFAC*SIN(THETA)*SIN(PHI)
      HI(2)=HIFAC*COS(THETA)*SIN(PHI)
      HI(3)=HIFAC*COS(PHI)
C
      DO 45 I=1,3
       HF(I)=HI(I)+HS(I)
   45 CONTINUE
C
      DO 50 I=1,3
       EFG(I)=EF(I)
       HFG(I)=HF(I)
   50 CONTINUE

      RETURN
      END

      SUBROUTINE FIELCT(WAVEL,REFRE2,REFIM2,
     1                  R,THETAG,PHIG, EFG,HFG)
C
C     coating field
C     BH p.182
C


      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=NMAX_COEFFS)
      INTEGER NSTOPF,NORD,I,N,N1
      REAL*8 WAVEL,REFRE2,REFIM2,R,THETAG,PHIG
      REAL*8 CC,OMEGA,MU2
      COMPLEX*16 EFG(3),HFG(3)
      COMMON /COEFF/A,B,C,D,F,G,V,W,NORD


      REAL*8 ZERO,ONE,TWO,PI,RN
      COMPLEX*16 CI,CREF2,WVK2,RHO2,BJN1,BYN1,ENCAP,PSIP,CHIP
      COMPLEX*16 HFFACT
      COMPLEX*16 BJ(NMAX),BY(NMAX)
      COMPLEX*16 A(NMAX),B(NMAX),C(NMAX),D(NMAX)
      COMPLEX*16 F(NMAX),G(NMAX),V(NMAX),W(NMAX)
      COMPLEX*16 VM1O1N(3),VM1E1N(3),VN1O1N(3),VN1E1N(3)
      COMPLEX*16 VM2O1N(3),VM2E1N(3),VN2O1N(3),VN2E1N(3)
      COMPLEX*16 EF(3),HF(3)



      REAL*8 THETA,PHI
      REAL*8 PIN(NMAX),TAU(NMAX)


      COMMON /BJ/BJ
      COMMON /BY/BY
      COMMON /LEG/PIN,TAU

#ifdef CHECK_SURFACE_MODE
      REAL*8 ESMBIG
      PARAMETER (ESMBIG=10.0)
      INTEGER MODEI,MODEP,MONMAX,MOIMAX,MOPMAX
      REAL*8 TEMP,EABS,EFCMP,EFCMAX
      COMPLEX*16 VM1OMX,VN1EMX,VM2OMX,VN2EMX
      COMPLEX*16 FMAX,GMAX,VMAX,WMAX
      REAL*8 ABSEF(3),ABSPR(4)

#ifdef CHECK_UNDERFLOW
      REAL*8 UFRAT
#endif

      EFCMAX=0.0
      MODEI=0
      MODEP=0
#endif /* CHECK_SURFACE_MODE */


      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0


C     PI=3.14159265D0
      PI=ACOS(-ONE)

      THETA=THETAG
      PHI=PHIG

C
C Wiscombe Criterion BH p.477,485; Barber-Hill p.196
C
      CREF2=DCMPLX(REFRE2,REFIM2)
      WVK2=TWO*PI*CREF2/WAVEL
      RHO2=WVK2*R
      NSTOPF=AINT(ABS(RHO2)+4.05D0*ABS(RHO2)**.3333D0+2.0D0)
      IF((NSTOPF.GT.NMAX-1).OR.(NSTOPF.LE.0)) THEN
       WRITE(*,*) 'FIELCT: NSTOPF(',NSTOPF,') > NMAX-1 (',
     * NMAX-1,')'
       STOP
      ELSE IF(NSTOPF.GT.NORD) THEN
       WRITE(*,*) 'FIELCT: NSTOPF(',NSTOPF,') > NORD (',
     * NORD,') at R=',R
       STOP
      END IF
C
      CALL BESJYD(RHO2,NSTOPF)
      CALL LEGEND(THETA,NSTOPF)
C
      CI=CMPLX(ZERO,ONE)
      DO 10 I=1,3
       EF(I)=CMPLX(ZERO,ZERO)
       HF(I)=CMPLX(ZERO,ZERO)
   10 CONTINUE
C
      DO 20 N=1,NSTOPF
       N1=N+1
C      RN=DBLE(N)
       RN=N*ONE
C
       PSIP=RHO2*BJ(N1-1)-RN*BJ(N1)
C
       BJN1=BJ(N1)
       VM1O1N(1)= CMPLX(ZERO,ZERO)
       VM1O1N(2)= COS(PHI)*PIN(N1)*BJN1
       VM1O1N(3)=-(SIN(PHI)*TAU(N1)*BJN1)
       VM1E1N(1)= CMPLX(ZERO,ZERO)
       VM1E1N(2)=-(SIN(PHI)*PIN(N1)*BJN1)
       VM1E1N(3)=-(COS(PHI)*TAU(N1)*BJN1)
       VN1O1N(1)= SIN(PHI)*RN*(RN+ONE)*SIN(THETA)*PIN(N1)*BJN1
     1           /RHO2
       VN1O1N(2)= SIN(PHI)*TAU(N1)*PSIP/RHO2
       VN1O1N(3)= COS(PHI)*PIN(N1)*PSIP/RHO2
       VN1E1N(1)= COS(PHI)*RN*(RN+ONE)*SIN(THETA)*PIN(N1)*BJN1
     1           /RHO2
       VN1E1N(2)= COS(PHI)*TAU(N1)*PSIP/RHO2
       VN1E1N(3)=-(SIN(PHI)*PIN(N1)*PSIP/RHO2)
C
       CHIP=RHO2*BY(N1-1)-RN*BY(N1)
C
       BYN1=BY(N1)
       VM2O1N(1)= CMPLX(ZERO,ZERO)
       VM2O1N(2)= COS(PHI)*PIN(N1)*BYN1
       VM2O1N(3)=-(SIN(PHI)*TAU(N1)*BYN1)
       VM2E1N(1)= CMPLX(ZERO,ZERO)
       VM2E1N(2)=-(SIN(PHI)*PIN(N1)*BYN1)
       VM2E1N(3)=-(COS(PHI)*TAU(N1)*BYN1)
       VN2O1N(1)= SIN(PHI)*RN*(RN+ONE)*SIN(THETA)*PIN(N1)*BYN1
     1           /RHO2
       VN2O1N(2)= SIN(PHI)*TAU(N1)*CHIP/RHO2
       VN2O1N(3)= COS(PHI)*PIN(N1)*CHIP/RHO2
       VN2E1N(1)= COS(PHI)*RN*(RN+ONE)*SIN(THETA)*PIN(N1)*BYN1
     1           /RHO2
       VN2E1N(2)= COS(PHI)*TAU(N1)*CHIP/RHO2
       VN2E1N(3)=-(SIN(PHI)*PIN(N1)*CHIP/RHO2)
C
       ENCAP=CI**RN*(TWO*RN+ONE)/(RN*(RN+ONE))
       DO 30 I=1,3
        EF(I)=EF(I)+ENCAP*(F(N)*VM1O1N(I)-CI*G(N)*VN1E1N(I)
     1                    +V(N)*VM2O1N(I)-CI*W(N)*VN2E1N(I))
        HF(I)=HF(I)+ENCAP*(G(N)*VM1E1N(I)+CI*F(N)*VN1O1N(I)
     1                    +W(N)*VM2E1N(I)+CI*V(N)*VN2O1N(I))
   30  CONTINUE

#ifdef CHECK_SURFACE_MODE

C seek origins (dominant/surface modes) of huge EFSQ in coating
C BH p.100, 326, 329

       DO 770 I=1,3
        ABSEF(I)=ABS(EF(I))
  770  CONTINUE
       EABS=SQRT(ABSEF(1)**2.0D0+ABSEF(2)**2.0D0+ABSEF(3)**2.0D0)
       IF(EABS.GT.ESMBIG) THEN

C       write(88,*) 'Huge EF! Seeking origins:',(EF(I),I=1,3)
C       write(88,*) 'ENCAP :',ENCAP
C       write(88,*) 'VM1O1N:',(VM1O1N(I),I=1,3)
C       write(88,*) 'VN1E1N:',(VN1E1N(I),I=1,3)
C       write(88,*) 'VM2O1N:',(VM2O1N(I),I=1,3)
C       write(88,*) 'VN2E1N:',(VN2E1N(I),I=1,3)
C       write(88,*) 'VM2O1N(2)= COS(PHI)*PIN(N1)*BYN1',
C    1              COS(PHI),PIN(N1),BYN1
C       write(88,*) 'VN2E1N(2)= COS(PHI)*TAU(N1)*CHIP/RHO2',
C    1              COS(PHI),TAU(N1),CHIP,RHO2
C       write(88,*) 'RHO1,RHO2',RHO1,RHO2
C       write(88,*) 'N FGVW:',N,F(N),G(N),V(N),W(N)

        TEMP=MAX(ABSEF(1),ABSEF(2),ABSEF(3))
        DO 771 I=1,3
         IF(TEMP.EQ.ABSEF(I)) MODEI=I
  771   CONTINUE
        ABSPR(1)=ABS(F(N)*VM1O1N(MODEI))
        ABSPR(2)=ABS(G(N)*VN1E1N(MODEI))
        ABSPR(3)=ABS(V(N)*VM2O1N(MODEI))
        ABSPR(4)=ABS(W(N)*VN2E1N(MODEI))
        TEMP=MAX(ABSPR(1),ABSPR(2),ABSPR(3),ABSPR(4))
        DO 772 I=1,4
         IF(TEMP.EQ.ABSPR(I)) MODEP=I
  772   CONTINUE
        EFCMP=ABS(ENCAP)*ABSPR(MODEP)
        IF(EFCMP.GT.EFCMAX) THEN
         EFCMAX=EFCMP
         MONMAX=N
         MOIMAX=MODEI
         MOPMAX=MODEP
         VM1OMX=VM1O1N(MODEI)
         VN1EMX=VN1E1N(MODEI)
         VM2OMX=VM2O1N(MODEI)
         VN2EMX=VN2E1N(MODEI)
         FMAX=F(N)
         GMAX=G(N)
         VMAX=V(N)
         WMAX=W(N)
        END IF

       END IF

#endif /* CHECK_SURFACE_MODE */

   20 CONTINUE

C
C electric field E [V m-1] = EF * E0
C magnetic field
C
C ASSUME NON-MAGNETIC (MU=MU0=const) [N A-2]
      MU2=4.0D0*PI*1.0D-7
C light speed [m s-1]
      CC=2.99792458D8
      OMEGA=2.0D0*PI*CC/(WAVEL*1.0D-6)
      HFFACT=-((WVK2*1.0D6)/(OMEGA*MU2))
      DO 40 I=1,3
       HF(I)=HFFACT*HF(I)
   40 CONTINUE
C
      DO 50 I=1,3
       EFG(I)=EF(I)
       HFG(I)=HF(I)
   50 CONTINUE

#ifdef CHECK_SURFACE_MODE

 9118 FORMAT(A28,1X,4('(',E12.5,',',E12.5,')'))
 9119 FORMAT(A28,1X,4('(',E12.5,13X,      ')'))

      IF(EFCMAX.GT.0.0) THEN

#ifdef CHECK_UNDERFLOW
       UFRAT=ABS(EABS/EFCMAX)
       CALL UFCHECK1(UFRAT,'FIELCT:CHECK_SURFACE_MODE')
#endif

       write(51,*) ''
       write(51,*) 'FIELCT - Large EF (Surface mode?) EABS = ',
     1 EABS,'(some En >',ESMBIG,') at (R,THETA,PHI):',
     2 R,THETAG,PHIG
       write(51,*) 'ABS(EF(R)),ABS(EF(THETA)),ABS(EF(PHI))',
     1             ABS(EFG(1)),ABS(EFG(2)),ABS(EFG(3))
       write(51,*) 'Underflow ratio (sum/max_term) =',EABS/EFCMAX
       write(51,*) 'at n, n_axis(e_r,e_theta,e_phi), ',
     1 'n_term(f*VM1,g*VN1,v*VM2,w*VN2)',MONMAX,MOIMAX,MOPMAX
       write(51,9118) 'VM1O1N,VN1E1N,VM2O1N,VN2E1N',
     1                 VM1OMX,VN1EMX,VM2OMX,VN2EMX
       write(51,9119) 'ABS(VM1,VN1,VM2,VN2)       ',
     1 ABS(VM1OMX),ABS(VN1EMX),ABS(VM2OMX),ABS(VN2EMX)
       write(51,9118) 'F,     G,     V,     W     ',
     1                 FMAX,GMAX,VMAX,WMAX
       write(51,9119) 'ABS(f,g,v,w)               ',
     1                 ABS(FMAX),ABS(GMAX),ABS(VMAX),ABS(WMAX)

      END IF

#endif /* CHECK_SURFACE_MODE */

      RETURN
      END

      SUBROUTINE BESJYD(X,NORD)
C       ************************************************
C       SPHERICAL BESSELS: jn(x), yn(x), Dn(x)
C       use SBESJH() for consistency

C        jn(x) (n = 0 up to nord)
C       BJ(N) = j_(N-1) (N=1 to NORD+1)
C       X: COMPLEX
C       CAUTION: jn(x) is NOT upward recurrence stable, 
C       but downward recurrence stable.  
C       upward r. causes serious error for x ~ 0.
C       BH p. 87 (4.13), 88, 128
C       Barber-Hill p.251
C       Abramowitz-Stegun p.452 (Miller Algorithm)

C        yn(x) (n = 0 up to nord)
C       BY(N) = y_(N-1) (N=1 to NORD+1)
C       X: COMPLEX
C       yn(x) is usually upward recurrence stable.
C       Pitfall: it is NEITHER upward NOR downward stable
C       if imag(X) is large, where as a function of n, yn 
C       decreases first, down to a small value, then increases
C       rapidly again, so round-offs are critical.  
C       note: yn(0) is negatively infinite but casts no 
C       problem, as it is not called for core field.
C       BH p. 87 (4.13), p.128

C       LOGARITHMIC DERIVATIVE D(J) CALCULATED BY DOWNWARD
C       RECURRENCE BEGINNING WITH INITIAL VALUE 0.0 + I*0.0
C       AT J = NMX
C        D_n(x) = PSI_n'(x)/PSI_n(x) (n = 0 up to nord)
C       BD(N) = D_(N-1) (N=1 to NORD+1)
C       X: COMPLEX
C       CAUTION: D_n(x) is NOT upward recurrence stable, 
C       but downward recurrence stable.  
C       BH p.127 (4.89), 128, 478
C       Barber-Hill p.250
C       ************************************************



      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=NMAX_COEFFS)
      REAL*8 UFTOL
      PARAMETER (UFTOL=1.0D-6)
      INTEGER NORD,N
      INTEGER LMAX,IFAIL
      PARAMETER (LMAX=NMAX_COEFFS)


      REAL*8 ZERO,ONE
      COMPLEX*16 X,CI,BJ(NMAX),BY(NMAX),BD(NMAX)
      COMPLEX*16 XJ(0:LMAX),XJP(0:LMAX),
     1           XH1(0:LMAX),XH1P(0:LMAX)


      COMMON /BJ/BJ
      COMMON /BY/BY
      COMMON /BD/BD


      ZERO=0.0D0
      ONE=1.0D0
      CI=DCMPLX(ZERO,ONE)


C
      IF((NORD.GT.NMAX-1).OR.(NORD.LE.0)) THEN
       WRITE(*,*) 'NORD(',NORD,') > NMAX-1 (',NMAX-1,')'
       STOP
      END IF

C Thompson's (case 2): seems OK even with standard Fortran!
      CALL SBESJH(X,NORD, XJ,XJP,XH1,XH1P,IFAIL)

      DO 10 N=1,NORD+1
       BJ(N)=XJ(N-1)
   10 CONTINUE

      DO 20 N=1,NORD+1
       BY(N)=(XH1(N-1)-XJ(N-1))/CI
   20 CONTINUE

      DO 30 N=1,NORD+1
       BD(N)=XJP(N-1)/XJ(N-1)+ONE/X
   30 CONTINUE

#ifdef DEBUG_BESSEL
      CALL BESDEB(X,NORD)
#endif

      RETURN
      END


#ifdef DEBUG_BESSEL

      SUBROUTINE BESDEB(X,NORD)
C       ************************************************
C       Debug Bessels: compare several routines
C       ************************************************


      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=NMAX_COEFFS)
      REAL*8 UFTOL
      PARAMETER (UFTOL=1.0D-6)
      INTEGER NORD,N,NST,NMX

      REAL*8 RN,ZERO,ONE,TWO,SMALL
      COMPLEX*16 X,CNORM,TEMP,TEM2
      COMPLEX*16 BJ(NMAX),BY(NMAX),BD(NMAX),T(3)
      COMPLEX*16 AJ(NMAX),AY(NMAX),AD(NMAX),AP(NMAX)

      COMMON /BJ/BJ
      COMMON /BY/BY
      COMMON /BD/BD

      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0
      SMALL=1.0D-35

C jn(x)

C check limits

      IF(ABS(X).LT.1.0D-3) THEN
       WRITE(51,*) 'JBESSE: ABS(X) < 1.0D-3; check downward ',
     1 'recurrence; expected -> j_0(0) = 1, j_n(0) = 0 (n > 0)'


       WRITE(51,*) 'X = ',X
       DO 60 N=1,NORD+1
        WRITE(51,*) 'j_',N-1,'(X) = ',BJ(N)
   60  CONTINUE


      END IF

C check other algorithms

C Case 0: upward recurrence NG
C
C     BJ(1)=SIN(X)/X
C     BJ(2)=BJ(1)/X-COS(X)/X
C     DO 10 N=3,NORD+1
C      RN=N-1
C      BJ(N)=(2.0*RN-1.0)*BJ(N-1)/X-BJ(N-2)
C  10 CONTINUE
C     DO 20 N=1,NORD+1
C      WRITE(55,*) 'j',N-1,'(',X,')=',BJ(N)
C  20 CONTINUE

C Case 1: downward recurrence


      NST=NORD+1+INT((101.0D0+REAL(X))**0.5D0)


      T(3)=CMPLX(ZERO,ZERO)
      T(2)=CMPLX(SMALL,ZERO)
      DO 30 N=NST-1,NORD,-1


       RN=DBLE(N)


       T(1)=(TWO*RN+ONE)*T(2)/X-T(3)
       T(3)=T(2)
       T(2)=T(1)
   30 CONTINUE
C T(2) = j_(nord-1) = BJ(NORD)
      AJ(NORD+1)=T(3)
      AJ(NORD)=T(2)
      DO 40 N=NORD-1,1,-1


       RN=DBLE(N)


       AJ(N)=(TWO*RN+ONE)*AJ(N+1)/X-AJ(N+2)
   40 CONTINUE
C normalize by j_0(x) = sin(x)/x
      CNORM=SIN(X)/X/AJ(1)
      DO 50 N=1,NORD+1
       AJ(N)=AJ(N)*CNORM
   50 CONTINUE
C
      DO 99 N=1,NORD+1
       TEMP=(AJ(N)-BJ(N))/BJ(N)
       IF(ABS(TEMP).GT.UFTOL) WRITE(88,*) 'NE! X,N,rdev(down-BJ)',


     1 X,N,ABS(TEMP)


   99 CONTINUE

C Case 2: Cantrell's

      AJ(1)=SIN(X)/X
      DO 11 N=2,NORD+1


       RN=DBLE(N-1)


       AJ(N)=AJ(N-1)/(BD(N)+RN/X)
       TEMP=(AJ(N)-BJ(N))/BJ(N)
       IF(ABS(TEMP).GT.UFTOL) WRITE(88,*) 'NE! X,N,rdev(Cantrell-BJ)',


     1 X,N,ABS(TEMP)


   11 CONTINUE

C check derivative recurrence consistency: PSIP=PSIB-RN*PSIN/X
#ifdef DEBUG_BESSEL_DERIV
      do 19 N=1,3
C      RL=L*ONE -> no need
       RN=N

       C=X*BJ(N)-RN*X*BJ(N+1)/X
       D=BJ(N+1)+X*XJP(N)
       te1=abs(C)
       te2=abs(C-D)
       write(51,*) 'JBESSE: L,PSIP1,dif',N,te1,te2
   19 continue
#endif /* DEBUG_BESSEL_DERIV */

C yn(x)

C limits

      IF(ABS(X).LT.1.0D-3) THEN
       WRITE(51,*) 'YBESSE: ABS(X) < 1.0D-3; check upward ',
     1 'recurrence; expected -> y_n(0) = -INF'


       WRITE(51,*) 'X = ',X
       DO 20 N=1,NORD+1
        WRITE(51,*) 'y_',N-1,'(X) = ',BY(N)
   20  CONTINUE


      END IF

C Case 1: upward reccurence (fails)


      TWO=2.0D0

      AY(1)=-(COS(X)/X)
      AY(2)=AY(1)/X-SIN(X)/X
      DO 100 N=3,NORD+1


       RN=DBLE(N-1)

       AY(N)=(TWO*RN-ONE)*AY(N-1)/X-AY(N-2)
       TEMP=(AY(N)-BY(N))/BY(N)
       IF(ABS(TEMP).GT.UFTOL) WRITE(88,*) 'NE! X,N,Upward-BY,rdev',


     1 X,N,TEM2,ABS(TEMP)

  100 CONTINUE

C Case 2: Thompson's (case 1): still subject to round-offs
      AY(1)=-(COS(X)/X)
      AP(1)=(SIN(X)+COS(X)/X)/X
      DO 110 N=2,NORD+1


       RN=DBLE(N-1)


       AY(N)=AY(N-1)*(RN-ONE)/X-AP(N-1)
       AP(N)=AY(N-1)-AY(N)*(RN+ONE)/X
       TEMP=(AY(N)-BY(N))/BY(N)
       IF(ABS(TEMP).GT.UFTOL) WRITE(88,*) 'NE! X,N,Thomp1-BY,rdev',


     1 X,N,AY(N),ABS(TEMP)


  110 CONTINUE

C output y_n(x) for round-off checking (im(x) != 0)

#ifdef OUTPUT_YBESSE

      IF(ABS(AIMAG(X)).GT.1.0D-3) THEN
       WRITE(88,*) 'OUTPUT_YBESSE: y_n(x) output beg -----------'


       WRITE(88,*) 'X = ',X
       DO 25 N=1,NORD+1
        WRITE(88,*) 'y_',N-1,'(X) = ',BY(N)
   25  CONTINUE


       WRITE(88,*) 'OUTPUT_YBESSE: y_n(x) output end -----------'
      END IF

#endif /* OUTPUT_YBESSE */

C Dn(x)

C limits

      IF(ABS(X).LT.1.0D-3) THEN
       WRITE(51,*) 'DBESSE: ABS(X) < 1.0D-3; check downward ',
     1 'recurrence; expected -> D_n(0) = INF'


       WRITE(51,*) 'X = ',X
       DO 220 N=1,NORD+1
        WRITE(51,*) 'D_',N-1,'(X) = ',BD(N)
  220  CONTINUE


      END IF

      NMX=AINT(MAX(DBLE(NORD),ABS(X)))+15

      IF((NMX.GT.NMAX-1).OR.(NMX.LE.0)) THEN
       WRITE(*,*) 'NMX(',NMX,') > NMAX-1 (',NMAX-1,')'
       STOP
      END IF
C
      AD(NMX+1)=CMPLX(ZERO,ZERO)
      DO 10 N=NMX,1,-1


       RN=DBLE(N)


       AD(N)=(RN/X)-(ONE/(AD(N+1)+RN/X))

   10 CONTINUE
C
      DO 15 N=1,NORD+1
       TEMP=(AD(N)-BD(N))/BD(N)
       IF(ABS(TEMP).GT.UFTOL) WRITE(88,*) 'NE! X,N,rdev(downw-BD)',


     1 X,N,ABS(TEMP)


   15 CONTINUE

C
      RETURN
      END

#endif /* DEBUG_BESSEL */


      SUBROUTINE LEGEND(TH,NORD)
C       ************************************************
C       ASSOCIATED LEGENDRE (m=1) COMPUTED BY UPWARD RECURRENCE
C        pi_n = P1n(cos(theta))/sin(theta)
C        tau_n = dP1n(cos(theta))/d(theta)
C       PIN(N) = pi_(N-1)
C       TAU(N) = tau_(N-1)
C       "No particular computational problems."
C       BH p. 94 (4.46) <- I believe you, Craig...
C       ************************************************



      IMPLICIT NONE
      INTEGER NMAX
      PARAMETER (NMAX=NMAX_COEFFS)
      INTEGER NORD,N


      REAL*8 TH,U,RN,ZERO,ONE,TWO
      REAL*8 PIN(NMAX),TAU(NMAX)


      COMMON /LEG/PIN,TAU


      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0


C
      IF((NORD.GT.NMAX-1).OR.(NORD.LE.0)) THEN
       WRITE(*,*) 'NORD(',NORD,') > NMAX-1 (',NMAX-1,')'
       STOP
      END IF
C
      U=COS(TH)
      PIN(1)=ZERO
      PIN(2)=ONE
      TAU(1)=ZERO
      TAU(2)=U
      DO 10 N=3,NORD+1
C      RN=DBLE(N-1)
       RN=(N-1)*ONE
       PIN(N)=((TWO*RN-ONE)*U*PIN(N-1)-RN*PIN(N-2))/(RN-ONE)
       TAU(N)=RN*U*PIN(N)-(RN+ONE)*PIN(N-1)
   10 CONTINUE

C     DO 20 N=1,NORD+1
C      WRITE(51,*) 'pi ',N-1,'(',TH,')=',PIN(N)
C      WRITE(51,*) 'tau',N-1,'(',TH,')=',TAU(N)
C  20 CONTINUE

      RETURN
      END

C Spherical Besssel routine that works well
C Ref:
C I.J. THOMPSON & A.R. BARNETT, 
C Computer Physics Communications 47 (1987) 245 - 257

      SUBROUTINE SBESJH(X,LMAX,XJ,XJP,XH1,XH1P,IFAIL)
C ***                                                       I.J.Thompson
C ***                                                       31 May 1985.
C ***  COMPLEX SPHERICAL BESSEL FUNCTIONS from l=0 to l=LMAX
C ***    for X in the UPPER HALF PLANE ( Im(X) > -3)
C ***
C ***    XJ(l)   = j/l(x)          regular solution: XJ(0)=sin(x)/x
C ***    XJP(l)  = d/dx j/l(x)
C ***    XH1(l)  = h(1)/l(x)       irregular Hankel function:
C ***    XH1P(l) = d/dx h(1)/l(x)            XH1(0) = j0(x) + i. y0(x)
C ***                                               =(sin(x)-i.cos(x))/x
C ***                                               = -i.exp(i.x)/x
C ***  Using complex CF1, and trigonometric forms for l=0 solutions.
C ***

      IMPLICIT NONE
      INTEGER LIMIT,IFAIL,LMAX,L
      COMPLEX*16 XDB,XH1LDB,XH1BDB,XW
      PARAMETER (LIMIT=20000)

#ifdef DEBUG_BESSEL_DERIV
      real*8 te1,te2
#endif

      COMPLEX*16 X,CI,XI,W,PL,F,B,D,C,DEL,XJ0
      COMPLEX*16 XJ(0:LMAX),XJP(0:LMAX),XH1(0:LMAX),
     1           XH1P(0:LMAX)
      REAL*8 ZERO,ONE,ACCUR,TM30,ABSC




      ZERO=0.0D0
      ONE=1.0D0
      ACCUR=1.0D-12
      TM30=1D-30

      CI=CMPLX(ZERO,ONE)

C      ABSC(W) = ABS(REAL(W)) + ABS(IMAG(W))
C      ABSC(W) = ABS(mpreal(W)) + ABS(AIMAG(W))
      IFAIL= -1
C      IF(ABSC(X).LT.ACCUR .OR. IMAG(X).LT.-3.0) GO TO 5


      ABSC=ABS(REAL(X))+ABS(AIMAG(X))


      IF(ABSC.LT.ACCUR .OR. AIMAG(X).LT.-3.0D0) GO TO 5
      XI = ONE/X
      W  = XI + XI


      PL = LMAX*XI


      F = PL + XI
      B  = F + F + XI
      D  = ZERO
      C  = F
      DO 1 L=1,LIMIT
       D  = B - D
       C  = B - ONE/C
C         IF(ABSC(D).LT.TM30) D = TM30
C         IF(ABSC(C).LT.TM30) C = TM30


       ABSC=ABS(REAL(D))+ABS(AIMAG(D))


       IF(ABSC.LT.TM30) D = TM30


       ABSC=ABS(REAL(C))+ABS(AIMAG(C))


       IF(ABSC.LT.TM30) C = TM30
       D = ONE / D
       DEL= D * C
       F = F * DEL
       B = B + W
C    1 IF(ABSC(DEL-ONE).LT.ACCUR) GO TO 2


       ABSC=ABS(REAL(DEL-ONE))+ABS(AIMAG(DEL-ONE))


       IF(ABSC.LT.ACCUR) GO TO 2
    1 CONTINUE
        IFAIL = -2
        GO TO 5
C
    2 XJ(LMAX)   = TM30
      XJP(LMAX)  = F * XJ(LMAX)
C
C *** Downward recursion to l=0 (N.B.  Coulomb Functions)
C
      DO 3 L = LMAX-1,0,-1
       XJ(L) = PL*XJ(L+1) + XJP(L+1)
       XJP(L)= PL*XJ(L)   - XJ(L+1)
    3 PL = PL - XI

C *** Calculate the l=0 Bessel Functions
      XJ0  = XI * SIN(X)
      XH1(0) = EXP(CI*X) * XI * (-CI)
      XH1P(0)= XH1(0) * (CI - XI)
C
C *** Rescale XJ, XJP,  converting to spherical Bessels.
C *** Recur   XH1,XH1P             AS spherical Bessels.
C
      W = ONE/XJ(0)
      PL = XI
      DO 4 L = 0,LMAX
       XJ(L)  =  XJ0*(W*XJ(L))
       XJP(L) =  XJ0*(W*XJP(L)) - XI*XJ(L)
       IF(L.EQ.0) GO TO 4
       XH1(L) = (PL-XI) * XH1(L-1) - XH1P(L-1)

C check if hankel is increasing (upward stable)
       IF(ABS(XH1(L)).LT.ABS(XH1(L-1))) THEN
        XDB=X
        XH1LDB=XH1(L)
        XH1BDB=XH1(L-1)
        WRITE(88,*) 'SBESJH: hankel not increasing; x ',
     1  XDB,',',L-1,XH1BDB,'->',L,XH1LDB
       END IF

       PL = PL + XI
C       XH1P(L)=- PL     * XH1(L)   + XH1(L-1)
       XH1P(L)=- (PL     * XH1(L))   + XH1(L-1)
    4 CONTINUE

C check derivative recurrence consistency: PSIP=PSIB-RN*PSIN/X
#ifdef DEBUG_BESSEL_DERIV
      do 99 L=1,3
C      RL=L*ONE -> no need
       ABSC=L

       C=X*XJ(L-1)-ABSC*X*XJ(L)/X
       D=XJ(L)+X*XJP(L)
       te1=abs(C)
       te2=abs(C-D)
       write(51,*) 'SBESJH: L,PSIP1,dif',L,te1,te2
   99 continue
#endif /* DEBUG_BESSEL_DERIV */

C success
      IFAIL = 0
      RETURN

C failure
    5 CONTINUE
      XW=X
      WRITE(6,*) 'SBESJH(X = ',XW,'): IFAIL = ',IFAIL
      stop
      RETURN
      END


#ifdef CHECK_UNDERFLOW

      SUBROUTINE EQCHECK(CX,CY,STR)
C check if CX.EQ.CY; STR is a comment for debugging.



      IMPLICIT NONE
      REAL*8 UFTOL
      PARAMETER (UFTOL=1.0D-6)
      REAL*8 R
      CHARACTER(LEN=*) STR


      COMPLEX*16 CX,CY,CMAX


      IF(ABS(CX).GT.ABS(CY)) THEN
       CMAX=CX
      ELSE
       CMAX=CY
      END IF
      R=ABS((CX-CY)/CMAX)
      IF(R.GT.UFTOL) WRITE(88,*) 'EQCHECK : rdev  = ',R,STR

      RETURN
      END

      SUBROUTINE UFCHECK1(R,STR)
C check underflow & round-off (if R is too small)
C STR is a comment for debugging.

      IMPLICIT NONE
      REAL*8 UFTOL
      PARAMETER (UFTOL=1.0D-6)
      REAL*8 R,RDMIN
      CHARACTER*80 RDSTR
      CHARACTER(LEN=*) STR
      common /RDMIN/RDMIN,RDSTR

      IF(R.LT.UFTOL) THEN
C      WRITE(88,*) 'UFCHECK1: ratio = ',R,STR
       IF(R.LT.RDMIN) THEN
        RDMIN=R
        RDSTR=STR
       END IF
      END IF

      RETURN
      END

      SUBROUTINE UFCHECK2(CX,CY,STR)
C check underflow & round-off (if CX-CY is too small)
C STR is a comment for debugging.


      IMPLICIT NONE
      REAL*8 R
      CHARACTER(LEN=*) STR


      COMPLEX*16 CX,CY,CMAX


      IF(ABS(CX).GT.ABS(CY)) THEN
       CMAX=CX
      ELSE
       CMAX=CY
      END IF
      R=ABS((CX-CY)/CMAX)
      CALL UFCHECK1(R,STR)
      RETURN
      END

      SUBROUTINE UFCHECK3(CX,CY,CZ,STR)
C check underflow & round-off (if CX+CY+CZ is too small)
C STR is a comment for debugging.

      IMPLICIT NONE
      REAL*8 R
      CHARACTER(LEN=*) STR

      COMPLEX*16 CX,CY,CZ,CMAX


      IF(ABS(CX).GT.ABS(CY)) THEN
       CMAX=CX
      ELSE
       CMAX=CY
      END IF
      IF(ABS(CZ).GT.ABS(CMAX)) THEN
       CMAX=CZ
      END IF
      R=ABS((CX+CY+CZ)/CMAX)
      CALL UFCHECK1(R,STR)
      RETURN
      END

#endif /* CHECK_UNDERFLOW */

