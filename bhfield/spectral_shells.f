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


      PROGRAM test_main

      IMPLICIT NONE
      REAL*8 UFTOL
      PARAMETER (UFTOL=1.0D-6)
      
      INTEGER I,NSTOPF,NSHELLS
      
      REAL*8 WLFAC(3),UABS_T
      
      REAL*8 RADCOR,RADCOT,WAVEL,PI,CC,EPSVAC,MU,OMEGA
      REAL*8 REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      REAL*8 X,Y,Y1,Y2,Y3,Y4,YMAX,EXTMAX,RMAX
      REAL*8 QEXT,QSCA,QBACK,QABS,MYQABS
      REAL*8 AP
      
      COMPLEX*16 RFREL1,RFREL2
      
      REAL*8 XPMIN(3),XPMAX(3)
      REAL*8 XCMIN(3),XCMAX(3)
      
      CHARACTER*80 DIRNAM,FILNAM(3),FNAME(3)
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
      WAVEL = 0.507D0
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
      PRINT *, "using this nk data for reference medium: ", FNAME(1)
      PRINT *, "got: ",REFMED," + i*",0.00," for ",WAVEL," micron"
      PRINT *, "using this nk data for core medium: ", FNAME(2)
      PRINT *, "got: ",REFRE1," + i*",REFIM1," for ",WAVEL," micron"
      PRINT *, "using this nk data for shell medium: ", FNAME(3)
      PRINT *, "got: ",REFRE2," + i*",REFIM2," for ",WAVEL," micron"
      
C     get relative refractive indices
      RFREL1=DCMPLX(REFRE1,REFIM1)/REFMED
      RFREL2=DCMPLX(REFRE2,REFIM2)/REFMED
  
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
      
      CALL IntegrateUABS(RADCOR,RADCOT,WAVEL,
     &                   REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     &                   EPSVAC,OMEGA,MU,UABS_T)
      
      PRINT *, "UABS_T = ", UABS_T
      
      AP = PI*RADCOT*RADCOT
      MYQABS = UABS_T/AP*1.E-6
      
      PRINT *, ''
      PRINT *, 'Check full integral over sphere against Qabs'
      PRINT *, '************RESULTS*********************'
      PRINT *, "*  Fancy Qabs = ", QABS
      PRINT *, "*  My Qabs =    ", MYQABS
      PRINT *, '****************************************'
      
      
      NSHELLS = 50
      PRINT *, ''
      PRINT *, "Finally time to integrate some shells!"
      CALL IntegrateShells(RADCOR,RADCOT,WAVEL,
     &                   REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     &                   EPSVAC,OMEGA,MU,NSHELLS)
      
      
      STOP
      END
