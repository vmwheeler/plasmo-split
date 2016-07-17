      PROGRAM BHFIELD


C       *******************************************************
C       COAT IS THE CALLING PROGRAM FOR BHCOAT, THE SUBROUTINE
C       THAT CALCULATES EFFICIENCIES FOR A COATED SPHERE.
C       FOR GIVEN RADII AND REFRACTIVE INDICES OF INNER AND
C       OUTER SPHERES, REFRACTIVE INDEX OF SURROUNDING
C       MEDIUM, AND FREE SPACE WAVELENGTH, COAT CALCULATES SIZE
C       PARAMETERS AND RELATIVE REFRACTIVE INDICES
C       *******************************************************
C
C                        **********CAUTION**********
C
C       BHCOAT SHOULD NOT BE USED FOR LARGE, HIGHLY ABSORBING
C       COATED SPHERES
C       X*REFIM1, X*REFIM2, AND Y*REFIM2 SHOULD BE LESS THAN ABOUT 30
C
C                        **********CAUTION**********
C
C
      IMPLICIT NONE
      REAL*8 UFTOL
      PARAMETER (UFTOL=1.0D-6)
      INTEGER I,J,K,NSTOPF,N,IWHERE,IW,NCOUNT,NTOTAL,
     1        EHAND,HHAND
      REAL*8 RADCOR,RADCOT,FCOR,WAVEL,PI,AVOGAD,QTOEP2
      REAL*8 REFMED,REFRE1,REFIM1,REFRE2,REFIM2,RNORM
      REAL*8 EPSDYE,CDYE,ALPHA,EXTMAX,RMAX,YMAX,Y1,Y2,Y3,Y4
      REAL*8 CC,OMEGAP,VF,OMEGA,CORR,EPSRE,EPSIM,PATH,EPSIMCOR
      REAL*8 X,Y,QEXT,QSCA,QBACK,QABS,EPSEXT,EPSABS
      REAL*8 EFSQ,HFSQ,UABS,EPSVAC,MU,RI,RIPLAN,EANGH,HANGH
      REAL*8 EFPHI,HFPHI,DPHASE,ELLIPT,EAZIM,HLLIPT,HAZIM,
     1       EHANG,DEVEHA,DEVINT,EANGHD,HANGHD
      REAL*8 DHFAC,DH,DIVS1,DIVS2,DIVS3,DIVSR,VFRAC
      INTEGER NGRID(3),NFPT(3)
      REAL*8 WLFAC(3),XP(3),SV(3),XPMIN(3),XPMAX(3),XPSTEP(3)
      REAL*8 MXESQ(3,4),MXUAB(3,4),MXRIN(3,4),MXEHA(3,4),MXINT(3,4),
     1       MXDVU(3,4),MXDVS(3,4),MXELL(3,4),MXEAZ(3,4),MXEAH(3,4),
     2       MXHLL(3,4),MXHAZ(3,4),MXHAH(3,4),MXDPH(3,4)
      REAL*8 VMXDIV(7)
      REAL*8 EFMAJ(3),EFMIN(3),HFMAJ(3),HFMIN(3)
      COMPLEX*16 RFREL1,RFREL2,EPSC,EPSCM,EPSC1,EPSC2
      COMPLEX*16 EC(3),HC(3)
      CHARACTER*10 CASE,SFIELD(4)
      CHARACTER*24 STRTIM
      CHARACTER*50 SARG
      CHARACTER*50 FNLOGF,FNZAXS,FNESQ(4),FNUAB(4),FNVER(4),
     1             FNVEL(4),FNVHR(4),FNVHL(4),FNVSV(4)
      CHARACTER*80 DIRNAM,FILNAM(3),FNAME(3)

#ifdef CHECK_NUMERICAL_POYNTING
      INTEGER M
      REAL*8 DEVUAB,DEVDVS
      REAL*8 SCMP(3,8),VMXDSX(3,8)
#endif

#ifdef CHECK_UNDERFLOW
      INTEGER MPDIGITR
      REAL*8 RDMIN
      CHARACTER*80 RDSTR
      common /RDMIN/RDMIN,RDSTR

      RDMIN=UFTOL
      RDSTR='Default'
#endif

C
C factor for numerical differentiation in calc of divS
C
      DHFAC=1.0D-3
C
C optical constants data: water (wl/um), silica (wl/A), silver (wl/A)
C wavelength[um] = wavelength_data * WLFAC
C

C     DIRNAM='C:/Windoc/Programming/Fortran/Scatter/src/'
C     DIRNAM='./'

#ifdef DATA_DIR
      DIRNAM=DATA_DIR
      PRINT *, DATA_DIR
#else
      DIRNAM='./'
#endif

      FILNAM(1)='Segelstein.txt'
      WLFAC(1)=1.0D0
      FILNAM(2)='SiO2_palik.nk'
      WLFAC(2)=1.0D-4
      FILNAM(3)='Ag_palik.nk'
      WLFAC(3)=1.0D-4
      DO 901 I=1,3
       FNAME(I)=DIRNAM(1:MAX(INDEX(DIRNAM,' ')-1,1))//
     1          FILNAM(I)(1:MAX(INDEX(FILNAM(I),' ')-1,1))
  901 CONTINUE
C
C output file names
C
      FNLOGF='bhfield.log'
      FNZAXS  ='EU_zax.txt'
      FNESQ(4)='E_0allf.dat'
      FNESQ(1)='E_1core.dat'
      FNESQ(2)='E_2coat.dat'
      FNESQ(3)='E_3exte.dat'
      FNUAB(4)='U_0allf.dat'
      FNUAB(1)='U_1core.dat'
      FNUAB(2)='U_2coat.dat'
C non-absorbing medium: external UABS = 0
      FNUAB(3)='U_3exte.dat'
      FNVER(4)='V_0Ereim.dat'
      FNVER(1)='V_1Ereim.dat'
      FNVER(2)='V_2Ereim.dat'
      FNVER(3)='V_3Ereim.dat'
      FNVEL(4)='V_0Eelli.dat'
      FNVEL(1)='V_1Eelli.dat'
      FNVEL(2)='V_2Eelli.dat'
      FNVEL(3)='V_3Eelli.dat'
      FNVHR(4)='V_0Hreim.dat'
      FNVHR(1)='V_1Hreim.dat'
      FNVHR(2)='V_2Hreim.dat'
      FNVHR(3)='V_3Hreim.dat'
      FNVHL(4)='V_0Helli.dat'
      FNVHL(1)='V_1Helli.dat'
      FNVHL(2)='V_2Helli.dat'
      FNVHL(3)='V_3Helli.dat'
      FNVSV(4)='V_0Poynt.dat'
      FNVSV(1)='V_1Poynt.dat'
      FNVSV(2)='V_2Poynt.dat'
      FNVSV(3)='V_3Poynt.dat'
C
      I=0


      IF((IARGC().NE.14).AND.(IARGC().NE.19)) THEN
       WRITE(*,*)
     1 'Usage: bhfield wl[um] rad_core[um] rad_coat[um]',
     2 'n_grid_x xspan_min[um] xspan_max[um]',
     3 'n_grid_y yspan_min[um] yspan_max[um]',
     4 'n_grid_z zspan_min[um] zspan_max[um]',
     5 'case Kreibig',
     6 '[ref_med ref_re1 ref_im1 ref_re2 ref_im2 (case=other)]'
       WRITE(*,*) 'case = nanoshell/liposome/HPC/barber/other'
       WRITE(*,*) 'Kreibig = 0.0 - 1.0 (used for Ag)'
       STOP
      END IF


      I=I+1
      CALL GETARG(I,SARG)
      READ(SARG,*) WAVEL
      I=I+1
      CALL GETARG(I,SARG)
      READ(SARG,*) RADCOR
      I=I+1
      CALL GETARG(I,SARG)
      READ(SARG,*) RADCOT

      DO 9876 N=1,3
       I=I+1
       CALL GETARG(I,SARG)
       READ(SARG,*) NGRID(N)
       I=I+1
       CALL GETARG(I,SARG)
       READ(SARG,*) XPMIN(N)
       I=I+1
       CALL GETARG(I,SARG)
       READ(SARG,*) XPMAX(N)
 9876 CONTINUE

      I=I+1
      CALL GETARG(I,SARG)
      READ(SARG,*) CASE
      I=I+1
      CALL GETARG(I,SARG)
      READ(SARG,*) FCOR
C
      IF(CASE(1:5).EQ.'other') THEN
       I=I+1
       CALL GETARG(I,SARG)
       READ(SARG,*) REFMED
       I=I+1
       CALL GETARG(I,SARG)
       READ(SARG,*) REFRE1
       I=I+1
       CALL GETARG(I,SARG)
       READ(SARG,*) REFIM1
       I=I+1
       CALL GETARG(I,SARG)
       READ(SARG,*) REFRE2
       I=I+1
       CALL GETARG(I,SARG)
       READ(SARG,*) REFIM2
      END IF

C avoid radcor=0; bessel diverges
      IF(RADCOR.LE.0.0D0) THEN
       WRITE(*,*) 
     1'Core radius must be positive (to avoid Bessel divergence)!'
       STOP
      END IF

      OPEN(51,FILE=FNLOGF,STATUS='UNKNOWN')

#if defined(DEBUG_BESSEL) || defined(CHECK_UNDERFLOW)
      OPEN(88,FILE='bhdebug.log',STATUS='UNKNOWN')


      WRITE(88,*) 'standard version'


#endif /* DEBUG_BESSEL or CHECK_UNDERFLOW */

      CALL FDATE(STRTIM)
      WRITE(51,*) 'Started: ',STRTIM
      WRITE(51,*) ''

      WRITE(51,*) 'Command line:'


      IF(IARGC().EQ.14) THEN
       WRITE(51,*) 'bhfield ',WAVEL,RADCOR,RADCOT,
     1             (NGRID(I),XPMIN(I),XPMAX(I),I=1,3),CASE,FCOR
      ELSE
       WRITE(51,*) 'bhfield ',WAVEL,RADCOR,RADCOT,
     1             (NGRID(I),XPMIN(I),XPMAX(I),I=1,3),CASE,FCOR,
     2             REFMED,REFRE1,REFIM1,REFRE2,REFIM2
      END IF
   11 FORMAT (/"COATED SPHERE SCATTERING: bhfield (version: ",A,
     1        1X,A,")"/)
      WRITE(51,11) BHFIELD_VERSION,'standard'


C       ******************************************************
C       REFMED = (REAL) REFRACTIVE INDEX OF SURROUNDING MEDIUM
C       ******************************************************
C water at 1.064um; k assumed to be zero
C     REFMED=1.320506
C       *********************************************
C       REFRACTIVE INDEX OF CORE = REFRE1 + I*REFIM1
C       REFRACTIVE INDEX OF COAT = REFRE2 + I*REFIM2
C       *********************************************
C silica at 1.064um
C     REFRE1=1.53413
C     REFIM1=0.0
C Ag at 1.064um
C     REFRE2=0.234
C     REFIM2=7.21
C
C  12 FORMAT ("# REFMED = ",F8.4/"# REFRE1 =",E14.6,
C    13X,"REFIM1 =",E14.6/"# REFRE2 =",E14.6,3X,"REFIM2 =",E14.6)
C     WRITE(*,12) REFMED,REFRE1,REFIM1,REFRE2,REFIM2
C       ********************************
C       RADCOR = RADIUS OF CORE
C       RADCOT = RADIUS OF COAT
C       RADCOR, RADCOT, WAVEL SAME UNITS
C       ********************************
C nanoshell = silica (D = 100nm) + Ag (20nm thick)
C     RADCOR=0.050
C     RADCOT=0.070
C
C single wavelength = 1064nm
C     WAVEL=1.064
C Kreibig correction factor: 0.0(none) - 1.0(full)
C     FCOR=1.0
C  13 FORMAT (5X,"CORE RADIUS =",F7.3,3X,"COAT RADIUS =",F7.3/
C    15X,"WAVELENGTH =",F7.4)
C     WRITE(*,13) RADCOR,RADCOT,WAVEL
   13 FORMAT ("Input parameters:"/
     1"WAVELENGTH [um] =",F7.4,3X,
     2"CORE RADIUS[um] =",F7.3,3X,"COAT RADIUS[um] =",F7.3/
     3"X-Span: ","Ngrid =",I4," min[um] =",F7.3," max[um] =",F7.3/
     4"Y-Span: ","Ngrid =",I4," min[um] =",F7.3," max[um] =",F7.3/
     5"Z-Span: ","Ngrid =",I4," min[um] =",F7.3," max[um] =",F7.3/
     6"CASE =",A,"Kreibig = ",F7.4)
      WRITE(51,13) WAVEL,RADCOR,RADCOT,
     1             (NGRID(I),XPMIN(I),XPMAX(I),I=1,3),CASE,FCOR
C
C constants
C
C     PI=3.14159265D0
      PI=ACOS(-1.0D0)
C light speed [m s-1]
      CC=2.99792458D8
C eps0[F m-1]
      EPSVAC=1.0D7/(4.0D0*PI*CC*CC)
C assume non-magnetic (MU=MU0=const) [N A-2]
      MU=4.0D0*PI*1.0D-7
C angular frequency [s-1]
      OMEGA=2.0D0*PI*CC/(WAVEL*1.0D-6)
C
C epsilon(ext) = Qext * PI * a * a * NA * 10^(-11) / ln(10)
C a [um], epsilon [M-1 cm-1]
      AVOGAD=6.02214199D23
      QTOEP2=PI*AVOGAD*1.0D-11/LOG(10.0D0)
C
C case: nanoshell; silica/Ag/water
C
      IF(CASE(1:9).EQ.'nanoshell') THEN
C get optical constants for water, silica & Ag
       CALL OPTCON(1,WAVEL,FNAME,WLFAC, 
     1             REFMED,REFRE1,REFIM1,REFRE2,REFIM2)
       IF(FCOR.EQ.0.0D0) GO TO 2001
C
C Kreibig mean free path correction (Ag): BH p.337; Kreibig74
C
C Plasma frequency of silver [rad s-1] (Kreibig)
       OMEGAP=1.38D16
C Fermi velocity of silver [m s-1] (Kreibig)
       VF=1.4D6
C correction coeff [um]: Kreibig assumes path(L) = radius(R) (eq 15)
       CORR=FCOR*OMEGAP*OMEGAP*VF/(OMEGA*OMEGA*OMEGA)*1.0D6
C complex dielectric func (relative permittivity) of silver (coating)
       EPSRE=REFRE2*REFRE2-REFIM2*REFIM2
       EPSIM=2.0D0*REFRE2*REFIM2
C correction: path = radcot - radcor
       PATH=RADCOT-RADCOR
       EPSIMCOR=EPSIM+CORR/PATH
   28  FORMAT("KREIBIG: EPSRE =",F7.2," EPSIM =",F5.2,
     1 " CORR[nm] =",F5.2," PATH[nm] =",F7.2," EPSIM(COR) =",F5.2)
       WRITE(51,28) EPSRE,EPSIM,CORR*1.0D3,PATH*1.0D3,EPSIMCOR
C complex refractive index, corrected
   29  FORMAT(A, ": REFRE2 = ",F5.2," REFIM2 = ",F5.2)
       WRITE(51,29) "n,k (coat) before Kreibig corr: ",REFRE2,REFIM2
       REFRE2=SQRT((SQRT(EPSRE*EPSRE+EPSIMCOR*EPSIMCOR)+EPSRE)/2.0D0)
       REFIM2=SQRT((SQRT(EPSRE*EPSRE+EPSIMCOR*EPSIMCOR)-EPSRE)/2.0D0)
       WRITE(51,29) "n,k (coat) after  Kreibig corr: ",REFRE2,REFIM2
 2001  CONTINUE
C
C case: liposome; water/lipid bilayer(=water+dye)/water
C to do: calc refre2 by Drude or Kramers-Kronig
C
      ELSE IF(CASE(1:8).EQ.'liposome') THEN
C get optical constants for water; only REFMED needed
       CALL OPTCON(1,WAVEL,FNAME,WLFAC,
     1             REFMED,REFRE1,REFIM1,REFRE2,REFIM2)
       REFRE1=REFMED
       REFIM1=0.0D0
       REFRE2=REFMED
C IR-165: abs_epsilon = 1.05e5 # {molar extinction coeff/M-1cm-1}
       EPSDYE=1.05D5
C conc = 0.4e-3   {average? dye conc/M}
C conc = 6.5e-3 # {dye conc inside bilayer/M}
C {calculated from vpc=1700A^3, dye/lipid = 1/150, Note p89} check?
       CDYE=6.5D-3
C BH p.29
C {linear absorption coefficient [m-1]: cf. CRC 12-126}
       ALPHA=EPSDYE*CDYE*LOG(10.0D0)*1.0D2
       REFIM2=ALPHA*WAVEL*1.0D-6/(4.0D0*PI)
C
C case: HPC; HPC(n=1.43)/HPC/water
C calc for simple nanoparticle: coating shell causes no artifact, ok
C cf. HPC/water/water shows artifact NG
C
      ELSE IF(CASE(1:3).EQ.'HPC') THEN
C get optical constants for water; only REFMED needed
       CALL OPTCON(1,WAVEL,FNAME,WLFAC,
     1             REFMED,REFRE1,REFIM1,REFRE2,REFIM2)
       REFRE1=1.43D0
       REFIM1=0.0D0
       REFRE2=1.43D0
       REFIM2=0.0D0
C
C case: barber; compare with Barber-Hill results for testing
C
      ELSE IF(CASE(1:6).EQ.'barber') THEN
       REFMED=1.0D0
       REFRE1=1.5D0
       REFIM1=0.0D0
       REFRE2=1.5D0
       REFIM2=0.0D0
C
C case: other; specify optical constants via command line args
C
      ELSE IF(CASE(1:6).EQ.'other') THEN
C
      ELSE
       WRITE(*,*) 'Error: unknown case',CASE
       STOP
      END IF
C
C      RFREL1=COMPLEX(REFRE1,REFIM1)/REFMED
C      RFREL2=COMPLEX(REFRE2,REFIM2)/REFMED
      RFREL1=DCMPLX(REFRE1,REFIM1)/REFMED
      RFREL2=DCMPLX(REFRE2,REFIM2)/REFMED
C
   12 FORMAT ("Ref Index (medium) = ",F8.4/
     1"Ref Index (core) n, k = ",E14.6,", ",E14.6/
     2"Ref Index (coat) n, k = ",E14.6,", ",E14.6/
     3"Wavelength in the medium [um] = ",E14.6/
     4"Relative Ref Index (core) n, k = ",E14.6,", ",E14.6/
     5"Relative Ref Index (coat) n, k = ",E14.6,", ",E14.6)
      WRITE(51,12) REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1 WAVEL/REFMED,
     2 REAL(RFREL1),AIMAG(RFREL1),REAL(RFREL2),AIMAG(RFREL2)

C Reflectance at normal incidence BH 32; to compare with Qback BH 123

 2113 FORMAT("Reflectance (at normal incidence) ",A,"->",A,
     1": R = ",F8.3)
      RNORM=ABS((1.0D0-RFREL1)/(1.0D0+RFREL1))**2.0D0
      WRITE(51,2113) 'med ','core',RNORM
      RNORM=ABS((1.0D0-RFREL2)/(1.0D0+RFREL2))**2.0D0
      WRITE(51,2113) 'med ','coat',RNORM
      RNORM=ABS((1.0D0-(RFREL2/RFREL1))/(1.0D0+(RFREL2/RFREL1)))
     1      **2.0D0
      WRITE(51,2113) 'core','coat',RNORM

C
      X=2.0D0*PI*RADCOR*REFMED/WAVEL
      Y=2.0D0*PI*RADCOT*REFMED/WAVEL

C dielectric functions
      EPSCM=DCMPLX(REFRE1*REFRE1-REFIM1*REFIM1,2.0D0*REFRE1*REFIM1)
      EPSC1=DCMPLX(REFRE2*REFRE2-REFIM2*REFIM2,2.0D0*REFRE2*REFIM2)
      EPSC2=DCMPLX(REFMED*REFMED,0.0D0)

      WRITE(51,*) 'Surface mode condition check (BH p.326-330):'
      WRITE(51,*) 'Small sphere: Close to -(n+1)/n (n=1,2,...)?  ',
     1'  m1**2      = ',RFREL1**2,', m2**2      = ',RFREL2**2,
     2', (m2/m1)**2 = ',(RFREL2/RFREL1)**2

      WRITE(51,*) 'Finite-size: Close to -1?  ',
     1'eps1/((2+12*x*x/5)*epsm) = ',
     2 EPSC1/((2.0D0+12.0D0*X*X/5.0D0)*EPSCM),
     3'eps2/((2+12*y*y/5)*epsm) = ',
     4 EPSC2/((2.0D0+12.0D0*Y*Y/5.0D0)*EPSCM)

      VFRAC=(RADCOR/RADCOT)**3
      WRITE(51,*) 'Small coated sphere: Close to -1?  ',
     1'2*f*(eps2-epsm)*(eps1-eps2)/((eps2+2*epsm)*(eps1+2*eps2)) = ',
     2 2*VFRAC*(EPSC2-EPSCM)*(EPSC1-EPSC2)
     3 /((EPSC2+2*EPSCM)*(EPSC1+2*EPSC2))

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

      EXTMAX=MAX(ABS(XPMIN(1)),ABS(XPMIN(2)),ABS(XPMIN(3)),
     1           ABS(XPMAX(1)),ABS(XPMAX(2)),ABS(XPMAX(3)))

C     RMAX=EXTMAX*SQRT(2.0D0)
      RMAX=EXTMAX*SQRT(3.0D0)
      Y4=2.0D0*PI*RMAX*REFMED/WAVEL
      YMAX=MAX(Y1,Y2,Y3,Y4)
      NSTOPF=INT(YMAX+4.05D0*YMAX**0.3333D0+2.0D0)
C
   14 FORMAT ("CORE SIZE PARAM = ",F8.3,", COAT SIZE",
     1" PARAM = ",F8.3,", NSTOP(estim) = ",I3)
      WRITE(51,14) X,Y,NSTOPF
C
C calculate coefficients
C
      CALL BHCOAT(X,Y,RFREL1,RFREL2,NSTOPF, QEXT,QSCA,QBACK)
C
C  67 FORMAT (/,1X,"QSCA =",E13.6,3X,"QEXT =",E13.6,3X,
C    1"QBACK =",E13.6//)
C     WRITE(*,67) QSCA,QEXT,QBACK
C
      QABS=QEXT-QSCA
      EPSEXT=QEXT*QTOEP2*RADCOT*RADCOT
      EPSABS=QABS*QTOEP2*RADCOT*RADCOT
C
   18 FORMAT('wavelength[nm],epsilon(ext),epsilon(abs),Qext,Qabs,',
     1'Qsca,Qback')
      WRITE(51,18)
   16 FORMAT('',F7.2,6E11.4)
      WRITE(51,16) WAVEL*1.0D3,EPSEXT,EPSABS,QEXT,QABS,QSCA,QBACK
C

  950 FORMAT('# ',A8,' field: ',A/'# x[um] y[um] z[um] ',A)
  955 FORMAT('# Fields on z-axis: z[um] ',A)

      SFIELD(1)='Core    '
      SFIELD(2)='Coat    '
      SFIELD(3)='External'
      SFIELD(4)='All     '

      DO 80 IW=1,4
       OPEN(10+IW,FILE=FNESQ(IW),STATUS='UNKNOWN')
       OPEN(20+IW,FILE=FNUAB(IW),STATUS='UNKNOWN')
       OPEN(30+IW,FILE=FNVER(IW),STATUS='UNKNOWN')
       OPEN(40+IW,FILE=FNVEL(IW),STATUS='UNKNOWN')
       OPEN(60+IW,FILE=FNVHR(IW),STATUS='UNKNOWN')
       OPEN(70+IW,FILE=FNVHL(IW),STATUS='UNKNOWN')
       OPEN(80+IW,FILE=FNVSV(IW),STATUS='UNKNOWN')
       WRITE(10+IW,950) SFIELD(IW),'EFSQ','EFSQ[dimensionless]'
       WRITE(20+IW,950) SFIELD(IW),'Uabs','UABS[F m-1 s-1]'
       WRITE(30+IW,950) SFIELD(IW),
     1 'vector electric field: snapshots [Re (t=0), Im (t=period/4)]',
     2 'ReEx ReEy ReEz ImEx ImEy ImEz'
       WRITE(40+IW,950) SFIELD(IW),
     1 'vector electric field: vibration ellipse (major & minor axes), '
     2 //'ellipticity, azimuth[deg], p-a angle (phi)[deg], '
     3 //'handedness angle[deg] & handedness',
     4 'E-majx majy majz minx miny minz ellipt azim phi ah hand'
       WRITE(60+IW,950) SFIELD(IW),
     1 'vector magnetic field: snapshots [Re (t=0), Im (t=period/4)]',
     2 'ReHx ReHy ReHz ImHx ImHy ImHz'
       WRITE(70+IW,950) SFIELD(IW),
     1 'vector magnetic field: vibration ellipse (major & minor axes), '
     2 //'ellipticity, azimuth[deg], p-a angle (phi)[deg], '
     3 //'E-&H-phase dif[deg], handedness angle[deg] & handedness',
     4 'H-majx majy majz minx miny minz ellipt azim phi pdif ah hand'
       WRITE(80+IW,950) SFIELD(IW),
     1 'Poynting vector <S>, EH angle, optical irradiance (intensity) '
     2 //'(norm<S>), I(plane), -div<S> (1st-3rd), UABS & DIVSR',
     3 'Sx Sy Sz EHang[deg] RI RIPLAN -DIVS1 -DIVS2 -DIVS3 UABS DIVSR'
   80 CONTINUE
      OPEN(15,FILE=FNZAXS,STATUS='UNKNOWN')
      WRITE(15,955) 'EFSQ[dimensionless] UABS[F m-1 s-1]'

C
C scan X-Z (polarized paralell), Y-Z (perpendicular), or X-Y-Z space
C and compute electric & magnetic fields and absorbed energy density
C
C arrays for recording max values & positions
      DO 90 IW=1,3
       NFPT(IW)=0
      DO 90 N=1,4
       MXESQ(IW,N)=0.0D0
       MXUAB(IW,N)=0.0D0
       MXRIN(IW,N)=0.0D0
       MXEHA(IW,N)=0.0D0
       MXINT(IW,N)=0.0D0
       MXDVU(IW,N)=0.0D0
       MXDVS(IW,N)=0.0D0
       MXELL(IW,N)=0.0D0
       MXEAZ(IW,N)=0.0D0
       MXEAH(IW,N)=0.0D0
       MXHLL(IW,N)=0.0D0
       MXHAZ(IW,N)=0.0D0
       MXHAH(IW,N)=0.0D0
       MXDPH(IW,N)=0.0D0
   90 CONTINUE
      DO 92 N=1,7
       VMXDIV(N)=0.0D0
   92 CONTINUE

C scan steps

      DO 9874 I=1,3
       IF(NGRID(I).GT.1) THEN
        XPSTEP(I)=(XPMAX(I)-XPMIN(I))/DBLE(NGRID(I)-1)
       ELSE
        XPSTEP(I)=0.0D0
       END IF
 9874 CONTINUE

C DH [um] = difference to calc divS derivatives
C     DH=(EXTMAX/DBLE(IMAX))*DHFAC

      DH=1.0D99
      DO 9872 I=1,3
       IF(XPSTEP(I).NE.0.0D0) THEN
        DH=MIN(DH,ABS(XPSTEP(I)))
       END IF
 9872 CONTINUE
      DH=DH*DHFAC

      NTOTAL=NGRID(1)*NGRID(2)*NGRID(3)
      NCOUNT=0

C
C big loop: scan x, y, z (x runs first)
C

      DO 100 I=1,NGRID(3)
       XP(3)=XPMIN(3)+DBLE(I-1)*XPSTEP(3)
       DO 110 J=1,NGRID(2)
        XP(2)=XPMIN(2)+DBLE(J-1)*XPSTEP(2)
        DO 120 K=1,NGRID(1)
         XP(1)=XPMIN(1)+DBLE(K-1)*XPSTEP(1)
C
         CALL FIELD(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1              RADCOR,RADCOT,XP, IWHERE,EC,HC)

         NFPT(IWHERE)=NFPT(IWHERE)+1
C
         NCOUNT=NCOUNT+1
         WRITE(*,"(I9,'/',I9,': Grid ',3I4,'[Region ',I1,']')")
     1   NCOUNT,NTOTAL,K,J,I,IWHERE
C
C vector electric field: snapshots [Re (t=0), Im (t=period/4)]
C electric field strength E*(E^*) [(V m-1)^2] = EFSQ * E0**2
C EFSQ: dimensionless
C
         EFSQ=ABS(EC(1))**2.0D0+ABS(EC(2))**2.0D0+ABS(EC(3))**2.0D0
C
C vector magnetic field: snapshots [Re (t=0), Im (t=period/4)]
C magnetic field strength H*(H^*) [(A m-1)^2] = HFSQ * E0**2
C HFSQ[(A V-1)^2] = [(A m-1)^2 (V m-1)^-2]  (unused)
C
         HFSQ=ABS(HC(1))**2.0D0+ABS(HC(2))**2.0D0+ABS(HC(3))**2.0D0
C
C absorbed energy per unit volume and time Ua [W m-3] = UABS * E0**2
C eps0[F m-1], omega[s-1], E0[V m-1]; [F]=[C V-1], [W]=[C V s-1]
C UABS: [F m-1 s-1] = [(W m-3) (V m-1)^-2]
C Our formula: Uabs = EPSVAC*OMEGA*REFRE*REFIM*EFSQ (nonmagnetic)
C
         IF(IWHERE.EQ.1) THEN
C
          UABS=EPSVAC*OMEGA*REFRE1*REFIM1*EFSQ
          EPSC=EPSC1
C
         ELSE IF(IWHERE.EQ.2) THEN
C
          UABS=EPSVAC*OMEGA*REFRE2*REFIM2*EFSQ
          EPSC=EPSC2
C
         ELSE
C
          UABS=0.0D0
          EPSC=EPSCM
C
         END IF
C
C optional: intensity for homogeneous plane wave: BH p.29
C mu[N A-2], eps0[F m-1], EPSC complex dielectric func (relative permittivity)
C I(plane) [W m-2], RIPLAN[W V-2] = [(W m-2)(V m-1)^-2]
C
         RIPLAN=0.5D0*REAL(SQRT(EPSC*EPSVAC/MU))*EFSQ
C
C Poynting vector and intensity (irradiance)
C
         CALL POYNTI(EC,HC, SV,RI,EHANG)
         EHANG=EHANG*180.0D0/PI
C check deviation from right angle
         DEVEHA=ABS(EHANG-90.0D0)
C check I(true) vs I(plane)
         IF(RI.NE.0.0D0) THEN
          DEVINT=ABS((RI-RIPLAN)/RI)
         ELSE
          DEVINT=0.0D0
         END IF

C
C optional: numerical derivative for Uabs = -div<S> (slow)
C

#ifdef CHECK_NUMERICAL_POYNTING

         CALL CALDIV(WAVEL,REFMED,REFRE1,REFIM1,REFRE2,REFIM2,
     1               RADCOR,RADCOT,XP,IWHERE,DH,SV, DIVS1,
     2               DIVS2,DIVS3,SCMP)
C relative evaluation to check round-off errors in numerical deriv
C divSr [dimensionless] = delta(S) / norm(S) = div(S) * DH / norm(S)
C should be large enough (e.g. 1d-8) compared with real*8 digits (15)
         IF(RI.NE.0.0D0) THEN
          DIVSR=DIVS1*DH*1.0D-6/RI
         ELSE
          DIVSR=0.0D0
         END IF
         IF(DIVSR.NE.0.0D0.AND.ABS(DIVSR).LT.1.0D-11) THEN
          WRITE(51,*) 'WARNING - possible round-off: ',
     1    'abs(DIVSR) < 1D-11',DIVSR,'at',(XP(N),N=1,3)
         END IF
C        IF(ABS(DIVS1).GT.ABS(VMXDIV(4))) THEN
         IF(ABS(DIVS3).GT.ABS(VMXDIV(6))) THEN
          VMXDIV(4)=-DIVS1
          VMXDIV(5)=-DIVS2
          VMXDIV(6)=-DIVS3
          VMXDIV(7)=DIVSR
          VMXDIV(1)=XP(1)
          VMXDIV(2)=XP(2)
          VMXDIV(3)=XP(3)
          DO 150 N=1,3
           DO 160 M=1,8
            VMXDSX(N,M)=SCMP(N,M)
  160      CONTINUE
  150     CONTINUE
         END IF
C check UABS vs -div<S>: should be equal
C absolute deviation
         DEVUAB=ABS(UABS+DIVS3)
C relative deviation
         IF(UABS.NE.0.0D0.AND.DIVS3.NE.0.0D0) THEN
          DEVDVS=DEVUAB/UABS
         ELSE
          DEVDVS=0.0D0
         END IF

         CALL RECMAX(IWHERE,XP,DEVUAB,MXDVU)
         CALL RECMAX(IWHERE,XP,DEVDVS,MXDVS)

#else /* not CHECK_NUMERICAL_POYNTING */

         DIVS1=0.0D0
         DIVS2=0.0D0
         DIVS3=0.0D0
         DIVSR=0.0D0

#endif /* CHECK_NUMERICAL_POYNTING */


C
C vector electric field: ellipse magnitudes and phase
C
         CALL ELLIPS(EC, EFMAJ,EFMIN,EFPHI,ELLIPT,EAZIM)
         EFPHI=EFPHI*180.0D0/PI
         EAZIM=EAZIM*180.0D0/PI
C fold angle to fit within 0-90 deg
         EAZIM=90.0D0-ABS(EAZIM-90.0D0)
         CALL HANDED(EC,SV, EANGH,EHAND)
C angle of handedness: near 0 or 180 deg; record deviation
         EANGHD=90.0D0-ABS(EANGH-90.0D0)
C
C vector magnetic field: ellipse magnitudes and phase
C
         CALL ELLIPS(HC, HFMAJ,HFMIN,HFPHI,HLLIPT,HAZIM)
         HFPHI=HFPHI*180.0D0/PI
         HAZIM=HAZIM*180.0D0/PI
         HAZIM=90.0D0-ABS(HAZIM-90.0D0)
         CALL HANDED(HC,SV, HANGH,HHAND)
         HANGHD=90.0D0-ABS(HANGH-90.0D0)
C optional: DPHASE = (E-phase - H-phase) diff [deg]
C phase factor (phi) = angle between major axis and Re(E or H)
C (Born & Wolf p.35)
C acos -> [0, pi]; fit angle to 0-180 deg
         DPHASE=ACOS(COS((EFPHI-HFPHI)*PI/180.0D0))*180.0D0/PI
C phase shift (phi) may be near 0, 90 or 180 deg; record deviation
         IF(DPHASE.GT.135.0D0) THEN
          DPHASE=180.0D0-DPHASE
         ELSE IF(DPHASE.GT.45.0D0) THEN
          DPHASE=ABS(DPHASE-90.0D0)
         END IF

C keep records
         CALL RECMAX(IWHERE,XP,EFSQ,  MXESQ)
         CALL RECMAX(IWHERE,XP,UABS,  MXUAB)
         CALL RECMAX(IWHERE,XP,RI,    MXRIN)
         CALL RECMAX(IWHERE,XP,DEVEHA,MXEHA)
         CALL RECMAX(IWHERE,XP,DEVINT,MXINT)
         CALL RECMAX(IWHERE,XP,ELLIPT,MXELL)
         CALL RECMAX(IWHERE,XP,EAZIM, MXEAZ)
         CALL RECMAX(IWHERE,XP,EANGHD,MXEAH)
         CALL RECMAX(IWHERE,XP,HLLIPT,MXHLL)
         CALL RECMAX(IWHERE,XP,HAZIM, MXHAZ)
         CALL RECMAX(IWHERE,XP,HANGHD,MXHAH)
         CALL RECMAX(IWHERE,XP,DPHASE,MXDPH)

C
C output data
C

C each domain (core, coat, or external) & all fields

  701    FORMAT(3E13.5,E13.5)
  706    FORMAT(3E13.5,6E13.5)
  711    FORMAT(3E13.5,6E13.5,F6.3,3F6.1,I3)
  712    FORMAT(3E13.5,6E13.5,F6.3,4F6.1,I3)
  713    FORMAT(3E13.5,3E13.5,F7.2,6E13.5,E11.3)

         DO 130 IW=1,4
          IF(IW.EQ.4.OR.IW.EQ.IWHERE) THEN
           WRITE(10+IW,701) (XP(N),N=1,3),EFSQ
           WRITE(20+IW,701) (XP(N),N=1,3),UABS
           WRITE(30+IW,706) (XP(N),N=1,3),
     1     (REAL(EC(N)),N=1,3),(AIMAG(EC(N)),N=1,3)
           WRITE(40+IW,711) (XP(N),N=1,3),
     1     (EFMAJ(N),N=1,3),(EFMIN(N),N=1,3),ELLIPT,
     2     EAZIM,EFPHI,EANGH,EHAND
           WRITE(60+IW,706) (XP(N),N=1,3),
     1     (REAL(HC(N)),N=1,3),(AIMAG(HC(N)),N=1,3)
           WRITE(70+IW,712) (XP(N),N=1,3),
     1     (HFMAJ(N),N=1,3),(HFMIN(N),N=1,3),HLLIPT,
     2     HAZIM,HFPHI,DPHASE,HANGH,HHAND
           WRITE(80+IW,713) (XP(N),N=1,3),(SV(N),N=1,3),
     1     EHANG,RI,RIPLAN,-DIVS1,-DIVS2,-DIVS3,UABS,DIVSR
          ELSE
           WRITE(10+IW,701) (XP(N),N=1,3),0.0D0
           WRITE(20+IW,701) (XP(N),N=1,3),0.0D0
           WRITE(30+IW,706) (XP(N),N=1,3),(0.0D0,N=1,6)
           WRITE(40+IW,711) (XP(N),N=1,3),(0.0D0,N=1,10),0
           WRITE(60+IW,706) (XP(N),N=1,3),(0.0D0,N=1,6)
           WRITE(70+IW,712) (XP(N),N=1,3),(0.0D0,N=1,11),0
           WRITE(80+IW,713) (XP(N),N=1,3),(0.0D0,N=1,11)
          END IF
  130    CONTINUE

C
C field on z-axis
C
         IF((ABS(XP(1)).LT.1.0D-12).AND.(ABS(XP(2)).LT.1.0D-12)) THEN
          WRITE(15,*) XP(3),EFSQ,UABS
         END IF

C end of scanning loop

  120   CONTINUE
  110  CONTINUE
  100 CONTINUE

C
   32 FORMAT(/'Definitions & units: '/
     1 'EFSQ [dimensionless] ',
     2 '= squared amplitude of complex E-field Ec*(Ec^*) / E0**2'/
     3 'HFSQ [(A V-1)^2] = [(A m-1)^2 (V m-1)^-2] ',
     4 '= squared amplitude of complex H-field Hc*(Ec^*) / E0**2'/
     5 'UABS [F m-1 s-1] = [A V-1 m-1] = [(W m-3) (V m-1)^-2] ',
     6 '= absorbed energy per unit volume and time Ua [W m-3] / E0**2'/
     7 'SV [A V-1] = [(W m-2) (V m-1)^-2] ',
     8 '= Poynting vector <S> [W m-2] / E0**2'/
     9 'RI [A V-1] = [(W m-2) (V m-1)^-2] ',
     A '= optical intensity (irradiance) (= norm of <S>) / E0**2'/
     B 'DIVS [A V-1 m-1] = [(W m-3) (V m-1)^-2] ',
     C '= div <S> / E0**2 (1st, 2nd & 3rd order approx.)'/
     D 'DIVSR [dimensionless] = relative divergence (round-off check) ',
     E '= delta(S) / norm(S) = div <S> * dx / norm <S>'/)
      WRITE(51,32)

   33 FORMAT(A48,3(' max at (',3E11.3,') =',E11.3:' ; '))
   36 FORMAT(A48,3(' max at (',3E11.3,') =',F11.3:' ; '))

      DO 172 IW=1,3
       WRITE(51,"(/'Max val & pos in ',A8,' field (',I9,' points):'/)") 
     1 SFIELD(IW),NFPT(IW)
       WRITE(51,33) 'EFSQ             ',(MXESQ(IW,N),N=1,4)
       WRITE(51,33) 'UABS             ',(MXUAB(IW,N),N=1,4)
       WRITE(51,33) 'Irradiance:      ',(MXRIN(IW,N),N=1,4)
       WRITE(51,36) 'Dev of ReE-ReH angle from 90 [deg]:',
     1              (MXEHA(IW,N),N=1,4)
       WRITE(51,36) 'Rel dev I(true) vs I(plane):',(MXINT(IW,N),N=1,4)

       WRITE(51,36) 'E-Ellipticity:   ',(MXELL(IW,N),N=1,4)
       WRITE(51,36) 'H-Ellipticity:   ',(MXHLL(IW,N),N=1,4)
       WRITE(51,36) 'E-Azimuth [deg]: ',(MXEAZ(IW,N),N=1,4)
       WRITE(51,36) 'H-Azimuth [deg]: ',(MXHAZ(IW,N),N=1,4)
       WRITE(51,36) 'Dev E-Angle handedness [deg]:',(MXEAH(IW,N),N=1,4)
       WRITE(51,36) 'Dev H-Angle handedness [deg]:',(MXHAH(IW,N),N=1,4)
       WRITE(51,36) 'Dev of E- & H-phi dif from 0, 90 or 180 [deg]:',
     1              (MXDPH(IW,N),N=1,4)

#ifdef CHECK_NUMERICAL_POYNTING
       WRITE(51,33) 'Abs dev Uabs vs -div<S> [A V-1 m-1]:',
     1              (MXDVU(IW,N),N=1,4)
       WRITE(51,36) 'Rel dev Uabs vs -div<S>:',(MXDVS(IW,N),N=1,4)
#endif

  172 CONTINUE

#ifdef CHECK_NUMERICAL_POYNTING

   34 FORMAT(/'Max ABS(divS1-3): at (',3E11.3,') =',3E11.3,' [DIVSR =',
     1       E10.2,', DHFAC =',E10.2,', DH [um] =',E10.2,']')
   35 FORMAT('Axis e',I1,': Sx [dx = -3DH to 3DH] =',7E11.3,
     1       ' -> dSx(3rd) =',E11.3)

      WRITE(51,34) (VMXDIV(N),N=1,7),DHFAC,DH
      WRITE(51,*) 'divS numerical derivatives - slope check:'
      DO 170 N=1,3
       WRITE(51,35) N,(VMXDSX(N,M),M=1,8)
  170 CONTINUE

#endif

      CALL FDATE(STRTIM)
      WRITE(51,*) ''
      WRITE(51,*) 'Ended: ',STRTIM

C
      DO 180 IW=1,4
       CLOSE(10+IW)
       CLOSE(20+IW)
       CLOSE(30+IW)
       CLOSE(40+IW)
       CLOSE(60+IW)
       CLOSE(70+IW)
       CLOSE(80+IW)
  180 CONTINUE
      CLOSE(15)
      CLOSE(51)

#ifdef CHECK_UNDERFLOW

      MPDIGITR=-INT(LOG10(RDMIN))+4
      WRITE(88,*) 'Underflow assessment: min ratio = ',RDMIN,
     1'(',RDSTR(1:MAX(INDEX(RDSTR,' ')-1,1)),')'
      WRITE(88,*) 'ARP mpdigit should be (at least) > ', 
     1            MPDIGITR

#endif /* CHECK_UNDERFLOW */

#if defined(DEBUG_BESSEL) || defined(CHECK_UNDERFLOW)
      CLOSE(88)
#endif

      STOP
      END
