      SUBROUTINE BHCOAT(XX,YY,RRFRL1,RRFRL2,QQEXT,QQSCA,QBACK,GSCA)
      IMPLICIT NONE

! Arguments:

      REAL GSCA,QBACK,QQEXT,QQSCA,XX,YY
      COMPLEX RRFRL1,RRFRL2

! Local variables:

      DOUBLE COMPLEX II
      PARAMETER(II=(0.0,1.D0))
      DOUBLE PRECISION DEL
      PARAMETER(DEL=1.D-8)

      INTEGER IFLAG,N,NSTOP
      DOUBLE PRECISION
     &   CHI0Y,CHI1Y,CHIY,EN,PSI0Y,PSI1Y,PSIY,QEXT,QSCA,RN,X,Y,YSTOP
      DOUBLE COMPLEX
     &   AMESS1,AMESS2,AMESS3,AMESS4,AN,AN1,ANCAP,
     &   BN,BN1,BNCAP,BRACK,
     &   CHI0X2,CHI0Y2,CHI1X2,CHI1Y2,CHIX2,CHIPX2,CHIPY2,CHIY2,CRACK,
     &   D0X1,D0X2,D0Y2,D1X1,D1X2,D1Y2,DNBAR,GNBAR,
     &   REFREL,RFREL1,RFREL2,
     &   XBACK,XI0Y,XI1Y,XIY,
     &   X1,X2,Y2

!***********************************************************************
!
! Subroutine BHCOAT calculates Q_ext, Q_sca, Q_back, g=<cos> 
! for coated sphere.
! All bessel functions computed by upward recurrence.
! Input:
!        X = 2*PI*RCORE*REFMED/WAVEL
!        Y = 2*PI*RMANT*REFMED/WAVEL
!        RFREL1 = REFCOR/REFMED
!        RFREL2 = REFMAN/REFMED 
! where  REFCOR = complex refr.index of core)
!        REFMAN = complex refr.index of mantle)
!        REFMED = real refr.index of medium)
!        RCORE = radius of core
!        RMANT = radius of mantle
!        WAVEL = wavelength of light in ambient medium

! returns:
!        QQEXT = C_ext/pi*rmant^2
!        QQSCA = C_sca/pi*rmant^2
!        QBACK = 4*pi*(dQ_sca/dOmega)
!              = "radar backscattering efficiency factor"
!        GSCA  = <cos(theta)> for scattered power
!
! Routine BHCOAT is taken from Bohren & Huffman (1983)
! extended by Prof. Francis S. Binkowski of The University of North
! Carolina at Chapel Hill to evaluate GSCA=<cos(theta)>
! History:
! 92.11.24 (BTD) Explicit declaration of all variables
! 00.05.16 (BTD) Added IMPLICIT NONE
! 12.04.10 (FSB) Modified by Prof. Francis S. Binkowski of
!                The University of North Carolina at Chapel Hill
!                to evaluate GSCA=<cos(theta)>
! 12.06.15 (BTD) Cosmetic changes
!***********************************************************************

      X=XX
      Y=YY
      RFREL1=RRFRL1
      RFREL2=RRFRL2
!         -----------------------------------------------------------
!              del is the inner sphere convergence criterion
!         -----------------------------------------------------------
      X1=RFREL1*X
      X2=RFREL2*X
      Y2=RFREL2*Y
      YSTOP=Y+4.*Y**0.3333+2.0
      REFREL=RFREL2/RFREL1
      NSTOP=YSTOP
!         -----------------------------------------------------------
!              series terminated after nstop terms
!         -----------------------------------------------------------
      D0X1=COS(X1)/SIN(X1)
      D0X2=COS(X2)/SIN(X2)
      D0Y2=COS(Y2)/SIN(Y2)
      PSI0Y=COS(Y)
      PSI1Y=SIN(Y)
      CHI0Y=-SIN(Y)
      CHI1Y=COS(Y)
      XI0Y=PSI0Y-II*CHI0Y
      XI1Y=PSI1Y-II*CHI1Y
      CHI0Y2=-SIN(Y2)
      CHI1Y2=COS(Y2)
      CHI0X2=-SIN(X2)
      CHI1X2=COS(X2)
      QSCA=0.0
      QEXT=0.0
      XBACK=(0.0,0.0)
      IFLAG=0
      DO N=1,NSTOP
         RN=N
         EN=RN
         PSIY=(2.0*RN-1.)*PSI1Y/Y-PSI0Y
         CHIY=(2.0*RN-1.)*CHI1Y/Y-CHI0Y
         XIY=PSIY-II*CHIY
         D1Y2=1.0/(RN/Y2-D0Y2)-RN/Y2
         IF(IFLAG.EQ.0)THEN

! calculate inner sphere ancap, bncap
!           and brack and crack

            D1X1=1.0/(RN/X1-D0X1)-RN/X1
            D1X2=1.0/(RN/X2-D0X2)-RN/X2
            CHIX2=(2.0*RN-1.0)*CHI1X2/X2-CHI0X2
            CHIY2=(2.0*RN-1.0)*CHI1Y2/Y2-CHI0Y2
            CHIPX2=CHI1X2-RN*CHIX2/X2
            CHIPY2=CHI1Y2-RN*CHIY2/Y2
            ANCAP=REFREL*D1X1-D1X2
            ANCAP=ANCAP/(REFREL*D1X1*CHIX2-CHIPX2)
            ANCAP=ANCAP/(CHIX2*D1X2-CHIPX2)
            BRACK=ANCAP*(CHIY2*D1Y2-CHIPY2)
            BNCAP=REFREL*D1X2-D1X1
            BNCAP=BNCAP/(REFREL*CHIPX2-D1X1*CHIX2)
            BNCAP=BNCAP/(CHIX2*D1X2-CHIPX2)
            CRACK=BNCAP*(CHIY2*D1Y2-CHIPY2)

! calculate convergence test expressions for inner sphere
! see pp 483-485 of Bohren & Huffman for definitions

            AMESS1=BRACK*CHIPY2
            AMESS2=BRACK*CHIY2
            AMESS3=CRACK*CHIPY2
            AMESS4=CRACK*CHIY2

         ENDIF ! test on iflag.eq.0

! now test for convergence for inner sphere
! all four criteria must be satisfied
! see p 484 of Bohren & Huffman

         IF(ABS(AMESS1).LT.DEL*ABS(D1Y2).AND.
     &      ABS(AMESS2).LT.DEL.AND.
     &      ABS(AMESS3).LT.DEL*ABS(D1Y2).AND.
     &      ABS(AMESS4).LT.DEL)THEN

! convergence for inner sphere

            BRACK=(0.,0.)
            CRACK=(0.,0.)
            IFLAG=1
         ELSE

! no convergence yet

            IFLAG=0

         ENDIF
         DNBAR=D1Y2-BRACK*CHIPY2
         DNBAR=DNBAR/(1.0-BRACK*CHIY2)
         GNBAR=D1Y2-CRACK*CHIPY2
         GNBAR=GNBAR/(1.0-CRACK*CHIY2)

! store previous values of an and bn for use in computation of 
! g=<cos(theta)>

         IF(N.GT.1)THEN
            AN1=AN
            BN1=BN
         ENDIF

! update an and bn

         AN=(DNBAR/RFREL2+RN/Y)*PSIY-PSI1Y
         AN=AN/((DNBAR/RFREL2+RN/Y)*XIY-XI1Y)
         BN=(RFREL2*GNBAR+RN/Y)*PSIY-PSI1Y
         BN=BN/((RFREL2*GNBAR+RN/Y)*XIY-XI1Y)

! calculate sums for qsca,qext,xback

         QSCA=QSCA+(2.0*RN+1.0)*(ABS(AN)*ABS(AN)+ABS(BN)*ABS(BN))
         XBACK=XBACK+(2.0*RN+1.0)*(-1.)**N*(AN-BN)
         QEXT=QEXT+(2.0*RN+1.0)*(DBLE(AN)+DBLE(BN))

! (FSB) calculate the sum for the asymmetry factor

         GSCA=GSCA+((2.0*EN+1.)/(EN*(EN+1.0)))*
     &        (REAL(AN)*REAL(BN)+IMAG(AN)*IMAG(BN))
         IF(N.GT.1)THEN
            GSCA=GSCA+((EN-1.)*(EN+1.)/EN)*
     &           (REAL(AN1)*REAL(AN)+IMAG(AN1)*IMAG(AN)+
     &            REAL(BN1)*REAL(BN)+IMAG(BN1)*IMAG(BN))
         ENDIF

! continue update for next iteration

         PSI0Y=PSI1Y
         PSI1Y=PSIY
         CHI0Y=CHI1Y
         CHI1Y=CHIY
         XI1Y=PSI1Y-II*CHI1Y
         CHI0X2=CHI1X2
         CHI1X2=CHIX2
         CHI0Y2=CHI1Y2
         CHI1Y2=CHIY2
         D0X1=D1X1
         D0X2=D1X2
         D0Y2=D1Y2
      ENDDO

! have summed sufficient terms
! now compute QQSCA,QQEXT,QBACK, and GSCA
      
      QQSCA=(2.0/(Y*Y))*QSCA
      QQEXT=(2.0/(Y*Y))*QEXT
      QBACK=(ABS(XBACK))**2
      QBACK=(1.0/(Y*Y))*QBACK
      GSCA=2.0*GSCA/QSCA
      RETURN
      END
