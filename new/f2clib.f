C--   G.W. 02/12/2008 Standalone code for structure functions used in
C--   MSTW 2008 fits.  Code copied directly from fitting program.
C--   Only neutral-current structure functions in this version,
C--   not charged-current structure function (to be done).
C--   Comments to Graeme Watt <watt(at)hep.ucl.ac.uk>
C--   or Robert Thorne <thorne(at)hep.ucl.ac.uk>.

C--   Input variables for MSTWNC:
C--     x = Bjorken-x value.
C--     q = sqrt(Q2) in GeV.
C--     IPN = 1, 2 for p or n.
C--     IORD = 0, 1, 2 for LO, NLO, NNLO (pass in COMMON block).
C--   Output variables: f2,f2c,f2b,fl,flc,flb.
      SUBROUTINE MSTWNC(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION FTEMP(5)
      INTEGER alphaSorder,alphaSnfmax
      DOUBLE PRECISION distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ
      COMMON/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
      COMMON/GRPTHY/FLAVOR
      COMMON/DYLAMB/xlam,S0
      COMMON/iordCommon/iord
      COMMON/GAUS96/XI(96),WI(96),XX(97),NTERMS ! G.W. 15/02/2007
!      DATA PI,PI2/3.14159,9.8696/

!      IF (IORD.EQ.2) THEN
!         CALL MSTWNCnnlo(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
!         RETURN
!      ELSE IF (IORD.NE.0.AND.IORD.NE.1) THEN
!         WRITE(6,*) "Error in MSTWNC, IORD = ",IORD
!         STOP
!      END IF
! 
!      Q2=q**2.d0                 ! photon virtuality
!      WW2=(1.-X)*Q2/X+0.88      ! 0.88 is square of proton mass 
! 
!      xlam = 0.3D0
!      Q02 = 1.D0
!      S0=dLOG(Q02/xlam**2.d0)
!      S=dLOG(dLOG(Q2/xlam**2.d0)/S0)
!      qsdt = 4.D0*mCharm**2
!      qsct = 4.D0*mBottom**2
!      epsc4=qsdt/q2
!      fpsc4=qsdt/ww2
!      epsc=epsc4/4.D0
!      epsb4=qsct/q2
!      fpsb4=qsct/ww2
!      epsb=epsb4/4.D0
!      FLAVOR=3.D0
!      FAC=FLAVOR
!      CF=4.D0/3.D0
!      enf=flavor
!      AL1=LOG(1.-X)
!      T=S0*EXP(S)
!      AL=ALPHA(T)/(4.D0*PI)
!      Schm=LOG(LOG(0.25*Qsdt/xlam**2)/S0)
!      Sbot=LOG(LOG(0.25*Qsct/xlam**2)/S0)
!      tchm=s0*exp(schm)
!      tbot=s0*exp(sbot)
!      alchm=alpha(tchm)/(4.D0*PI)
!      albot=alpha(tbot)/(4.D0*PI)
!      ca=3.D0
!      al39=al
! 
!
! 
!      CALL FETCH(X,S,IPN,FTEMP)
!      theory=ftemp(1)
!      fx=ftemp(1)
!      fg=ftemp(2)
!      fc=ftemp(3) !4/9*(c+barc)
!      fb=ftemp(4)
!      fsing=ftemp(5)
!      sfac=2./9.
!      ffx=0.
!      ffc=0.
!      ffb=0.
!      fflx=0.
!      fflc=0.
!      fflb=0.
!      if(epsc.gt.1.) fc=0.
!      if(epsb.gt.1.) fb=0.
! 
!      IF3=0
!
!      FAC=6./9.
!      facc=4./9.
!      facb=1./9.
! 
!      IF (IORD.EQ.0) THEN
!      ffx=FX
!      ELSE
!      ffx=FX+FX*AL39*CF*(-9.-2.*PI2/3.+AL1*(-3.+2.*AL1))
!c$$$      fflx=fflx+fx*al**2*(-0.012)
!      fflx=fflx+fx*al**2*CLNN2C(x) ! G.W. 02/11/2007
!      END IF
!      
! 
!      DO 23 I=1,NTERMS
!      Y=0.5*(1.-X)*XI(I)+0.5*(1.+X)
!      XY=X/Y
!      AL1=LOG(1.-Y)
!      CALL FETCH(XY,S,IPN,FTEMP)
!      if(epsc.gt.1.) ftemp(3)=0.
!      if(epsb.gt.1.) ftemp(4)=0.
!      IF (IORD.NE.0) THEN
!      C22=CF*(6.+4.*Y-2.*(1.+Y*Y)/(1.-Y)*LOG(Y)-2.*(1.+Y)*AL1
!     2-IF3*2.*(1.+Y))
!      C23=CF*(-3.+4.*AL1)/(1.-Y)
!      CG2=2.*FAC*(-1.+8.*Y*(1.-Y)+(1.-2.*Y+2.*Y*Y)*LOG(1./Y-1.))
!      END IF
!      f1lq=4.*cf*y
!      f1lg=8.*fac*y*(1.-y)
!      IF (IORD.NE.0) THEN
!      ffx=ffx+.5*(1.-X)*WI(I)*AL39*
!     2(C22*FTEMP(1)+C23*(FTEMP(1)-FX))
!      ffx=ffx+.5*(1.-X)*WI(I)*AL39*CG2*FTEMP(2)
!      END IF
!      fflx=fflx+.5*(1.-x)*wi(i)*al*f1lq*ftemp(1)
!      fflx=fflx+.5*(1.-x)*wi(i)*al*f1lg*ftemp(2)
!
!      IF (IORD.NE.0) THEN
!      Y1=1.-Y
!      DL=LOG(Y)
!      DL2=DL*DL
!      DLM1=LOG(Y1)
!      DLM2=DLM1*DLM1
!      DLM3=DLM2*DLM1
!      DLM4=DLM3*DLM1
!c$$$      FNS2LQ=128./9.d0*y*DLM1**2-46.50*y*DLM1-84.094*DL*DLM1
!c$$$     x-37.338 +89.53*y
!c$$$     X+33.82*y**2+y*DL*(32.90+18.41*DL)-128./9.d0*DL
!c$$$     X+16./27.d0*flavor*(6.*y*DLM1-12*y*DL-25.*y+6.)
!      FNS2LQ = CLNN2A(Y,INT(flavor)) ! G.W. 02/11/2007
!      FS2LQ=fac*((15.94-5.212*y)*Y1*Y1*DLM1+(0.421+1.520*y)*DL*DL
!     X+28.09*Y1*DL-(2.370/Y-19.27)*Y1**3)
!      F2LG=fac*((94.74-49.20*y)*y1*DLM2+864.8*Y1*DLM1
!     x+1161.*y*DLM1*DL
!     X+60.06*y*DL*DL+39.66*Y1*DL-5.333*(1./Y-1.))
!      fflx=fflx+0.5*(1.-x)*WI(I)*AL*AL*
!     X(FNS2LQ*ftemp(1)+FS2LQ*ftemp(5)+F2LG*ftemp(2))
!      END IF
!
!   23 CONTINUE
!   21 CONTINUE
!
!      xcmax=1./(1.+epsc4)
!      if(xcmax.le.x) go to 321
!      xcmup=x/xcmax
!      CALL FETCH(Xcmup,S,IPN,FTEMP)
!      if(epsc.gt.1.) ftemp(3)=0.
!      fc=ftemp(3) 
!      IF (IORD.EQ.0) THEN
!      ffc=fc
!      ELSE
!      AL1c=LOG(1.-xcmup) 
!c      fc = 4./9.*(1.-xcmup)**3. !originally not here. toy pdf
!      ffc=Fc+Fc*AL*CF*(-9.-2.*PI2/3.+AL1c*(-3.+2.*AL1c)) !Contains the delta of C2q^1 and the terms of the + prescription (see my thesis)
!      END IF
!      
! 
!      DO 323 I=1,NTERMS
!      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
!      XY=X/Y
!      CALL FETCH(XY,S,IPN,FTEMP)
!      if(epsc.gt.1.) ftemp(3)=0.
!      gluxy=ftemp(2)
!      fcxy=ftemp(3) !its 4/9*(c+barc)!
!      ycmup=y/xcmax
!c      if(ycmup.gt.0.9999) ycmup=0.9999
!      IF (IORD.NE.0) THEN
!      c0c=1.
!      p0qg=ycmup**2+(1.-ycmup)**2
!      C22c=CF*(6.+4.*Ycmup-2.*(1.+Ycmup*Ycmup)/(1.-Ycmup)*
!     .LOG(Ycmup)-2.*(1.+Ycmup)*log(1.-ycmup)
!     2-IF3*2.*(1.+Ycmup))  !it's part of C2q^1 see my thesis
!      C23c=CF*(-3.+4.*log(1.-ycmup))/(1.-Ycmup) !it's part of C2q^1 (see my thesis)
!      if(epsc.gt.1.d0) c0c=0.d0
!      cg21c=2.*facc*cheavy(1,y,epsc)
!      cg22c=2.*facc*c0c*p0qg*log(1./epsc)
!      END IF
!      clg2c=2.*facc*cheavy(7,y,epsc)
!      f1lq=cheavy(8,ycmup,epsc)
!      IF (IORD.NE.0) THEN
!c      gluxy=(1.-x/y)**3. !originally not here. toy pdf
!c      fcxy =4./9.*(1.-x/y)**3. !originally not here. its 4/9*(c+barc).toy pdf
!      ffc=ffc+0.5/xcmax*(xcmax-X)*WI(I)*AL*(C22c*fcxy+C23c*(fcxy-fc)) 
!      ffc=ffc+0.5*(xcmax-x)*wi(i)*al*(cg21c-cg22c/xcmax)*gluxy 
!      END IF
!      CALL FETCH(X,S,IPN,FTEMP)
!c      gluxy=(1.-x/y)**3. !originally not here. toy pdf
!c      fcxy =4./9.*(1.-x/y)**3. !originally not here. its 4/9*(c+barc).toy pdf
!      fflc=fflc+0.5/xcmax*(xcmax-x)*wi(i)*al*f1lq*fcxy
!      fflc=fflc+0.5*(xcmax-x)*wi(i)*al*clg2c*gluxy 
!  323 CONTINUE
!      print*, x, Q2,fflc
!
!      if(epsc.gt.1.) then
!      xcmax=1./(1.+epsc4)
!      eps=epsc
!      if(xcmax.le.x) go to 321
!      DO 324 I=1,NTERMS
!      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
!      XY=X/Y
!      CALL FETCH(XY,S,IPN,FTEMP)
!      gluxy=ftemp(2)
!      IF (IORD.EQ.0) THEN
!      cg21c=2.*facc*cheavy(1,y,eps)
!      ffc=ffc+0.5*(xcmax-x)*wi(i)*al*(cg21c)*gluxy
!      ELSE
!      singxy=ftemp(5)
!      cgff2=facc*(c2gffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .c2gffnsh(y,eps)*0.5*exp(1-eps**2))
!      cqff2=facc*(c2qffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .c2qffnsh(y,eps)*0.5*exp(1-eps**2))
!C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
!      ffc=ffc+0.5*(xcmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*
!     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) )) !originally not commented
!      END IF
!
!  324 CONTINUE
!
!      else
!      xcmax=1.D0/(1.D0+ 4.D0)
!      eps=1.
!      if(xcmax.le.x) go to 321
!      DO 325 I=1,NTERMS
!      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
!      XY=X/Y
!      CALL FETCH(XY,Schm,IPN,FTEMP)
!      gluxy=ftemp(2)
!      IF (IORD.EQ.0) THEN
!      cg21c=2.*facc*cheavy(1,y,eps)
!      ffc=ffc+0.5*(xcmax-x)*wi(i)*alchm*(cg21c)*gluxy
!      ELSE
!      singxy=ftemp(5)
!      cgff2=facc*(c2gffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .c2gffnsh(y,eps)*0.5*exp(1-eps**2))
!      cqff2=facc*(c2qffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .c2qffnsh(y,eps)*0.5*exp(1-eps**2))
!C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
!      ffc=ffc+0.5*(xcmax-x)*wi(i)*alchm**2.*(cgff2*gluxy+cqff2*
!     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) )) !originally not commented
!      END IF
!
!  325 CONTINUE
!      endif
!
!
!      IF (IORD.NE.0) THEN
!      xcmax=1./(1.+epsc4)
!      if(xcmax.le.x) go to 321
!      xcmup=x/xcmax
!      CALL FETCH(XCMUP,S,IPN,FTEMP)
!      if(epsc.gt.1.) ftemp(3)=0. ! G.W. 05/11/2007
!c      fflc=fflc+ftemp(3)*(AL**2*CLNN2C(xcmup)) 
!c     &     *1.25*(1/(1+4.*epsc)-0.2) ! originally not commented
!      fflc=fflc !originally not here
!
! 529  continue
!      if(epsc.gt.1.) then
!      xcmax=1./(1.+epsc4)
!      eps=epsc
!      if(xcmax.le.x) go to 321
!      DO 524 I=1,NTERMS
!      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
!      XY=X/Y
!      CALL FETCH(XY,S,IPN,FTEMP)
!      gluxy=ftemp(2)
!      singxy=ftemp(5)
!      cgffl=facc*(clgffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .clgffnsh(y,eps)*0.5*exp(1-eps**2))
!      cqffl=facc*(clqffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .clqffnsh(y,eps)*0.5*exp(1-eps**2))
!C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
!c      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl* 
!c     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) )) !originally not commented
!       fflc=fflc !originally not here
!
!  524 CONTINUE
!      else
!      xcmax=(1./(1.+epsc4))
!      eps=epsc
!      if(xcmax.le.x) go to 321
!      DO 525 I=1,NTERMS
!      Y=0.5*(xcmax-X)*XI(I)+0.5*(xcmax+X)
!      XY=X/Y
!      CALL FETCH(XY,S,IPN,FTEMP)
!      if(epsc.gt.1.) ftemp(3)=0.
!      gluxy=ftemp(2)
!      fcxy=ftemp(3)
!      singxy=ftemp(5)
!      ymul=y*(1+epsc4)
!      Y1mul=1.-Ymul
!      DL=LOG(Ymul)
!      DL2=DL*DL
!      DLM1=LOG(Y1mul)
!      DLM2=DLM1*DLM1
!      DLM3=DLM2*DLM1
!      DLM4=DLM3*DLM1
!c$$$      FNS2LQmul=128./9.d0*ymul*DLM1**2-46.50*ymul*DLM1-84.094*DL*DLM1
!c$$$     x-37.338 +89.53*ymul
!c$$$     X+33.82*ymul**2+ymul*DL*(32.90+18.41*DL)-128./9.d0*DL
!c$$$     X+16./27.d0*enf*(6.*ymul*DLM1-12*ymul*DL-25.*ymul+6.)
!      FNS2LQmul = CLNN2A(YMUL,INT(enf))   ! G.W. 02/11/2007
!      FS2LQmul=((15.94-5.212*ymul)*Y1mul*Y1mul*DLM1+(0.421+1.520*ymul)
!     x*DL*DL+28.09*Y1mul*DL-(2.370/Ymul-19.27)*Y1mul**3)
!      cgvfl=facc*((clgffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
!     .clgffnsl(y,eps)*0.5*exp(1-1/eps**2))-clgvfsub(ymul,eps)/xcmax)
!      cqvfl=facc*((clqffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
!     .clqffnsl(y,eps)*0.5*exp(1-1/eps**2)))
!c      fflc=fflc+0.5/xcmax*(xcmax-x)*wi(i)*al**2.*(fcxy*FNS2LQmul+
!c     xfcxy*FS2LQmul)*1.25*(1/(1+4.*eps)-0.2) !originally not commented
!c      fflc=fflc+0.5*(xcmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*
!c     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) )) !originally not commented
!      fflc=fflc !originally not here
!  525 CONTINUE
!      endif
!      END IF
!
!
!  321 CONTINUE
!
!c      IF (IORD.NE.0) THEN
!c      if(ffc.lt.0.) ffc=0.
!c      END IF
!
!      xbmax=1./(1.+epsb4)
!      if(xbmax.le.x+0.00001) go to 421
!      xbmup=x/xbmax
!      CALL FETCH(Xbmup,S,IPN,FTEMP)
!      if(epsb.gt.1.) ftemp(4)=0.
!      fb=ftemp(4)
!      IF (IORD.EQ.0) THEN       ! G.W. 05/11/2008
!      ffb=fb                    ! G.W. 05/11/2008
!      ELSE                      ! G.W. 05/11/2008
!      AL1b=LOG(1.-xbmup)
!      ffb=Fb+Fb*AL*CF*(-9.-2.*PI2/3.+AL1b*(-3.+2.*AL1b))
!      END IF
!
!      DO 423 I=1,NTERMS
!      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
!      XY=X/Y
!      CALL FETCH(XY,S,IPN,FTEMP)
!      if(epsb.gt.1.) ftemp(4)=0.
!      gluxy=ftemp(2)
!      fbxy=ftemp(4)
!      ybmup=y/xbmax
!      if(ybmup.gt.0.99999d0) ybmup=0.99999d0
!      IF (IORD.NE.0) THEN
!      c0b=1.
!      p0qg=ybmup**2+(1.-ybmup)**2
!      C22b=CF*(6.+4.*Ybmup-2.*(1.+Ybmup*Ybmup)/(1.-Ybmup)*
!     .LOG(Ybmup)-2.*(1.+Ybmup)*log(1.-ybmup)
!     2-IF3*2.*(1.+Ybmup))
!      C23b=CF*(-3.+4.*log(1.-ybmup))/(1.-Ybmup)
!      if(epsb.gt.1.d0) c0b=0.d0
!      cg21b=2.*facb*cheavy(1,y,epsb)
!      cg22b=2.*facb*c0b*p0qg*log(1./epsb)
!      ffb=ffb+0.5/xbmax*(xbmax-X)*WI(I)*AL*(C22b*fbxy+C23b*(fbxy-fb))
!      ffb=ffb+0.5*(xbmax-x)*wi(i)*al*(cg21b-cg22b/xbmax)*gluxy
!      END IF
!      clg2b=2.*facb*cheavy(7,y,epsb)
!      f1lq=cheavy(8,ybmup,epsb)
!      fflb=fflb+0.5/xbmax*(xbmax-x)*wi(i)*al*f1lq*fbxy
!      fflb=fflb+0.5*(xbmax-x)*wi(i)*al*clg2b*gluxy
!
!
!  423 CONTINUE
!
!      if(epsb.gt.1.) then
!      xbmax=1./(1.+epsb4)
!      eps=epsb
!      if(xbmax.le.x) go to 421
!      DO 424 I=1,NTERMS
!      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
!      XY=X/Y
!      CALL FETCH(XY,S,IPN,FTEMP)
!      gluxy=ftemp(2)
!      IF (IORD.EQ.0) THEN
!      cg21b=2.*facb*cheavy(1,y,eps)
!      ffb=ffb+0.5*(xbmax-x)*wi(i)*al*(cg21b)*gluxy
!      ELSE
!c$$$      singxy=ftemp(5)+9./8.*ftemp(3) ! Why 9/8 and not 9/4?
!      singxy=ftemp(5)           ! G.W. 27/07/2007
!      cgff2=facb*(c2gffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .c2gffnsh(y,eps)*0.5*exp(1-eps**2))
!      cqff2=facb*(c2qffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .c2qffnsh(y,eps)*0.5*exp(1-eps**2))
!c$$$      ffb=ffb+0.5*(xbmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*singxy)
!C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
!      ffb=ffb+0.5*(xbmax-x)*wi(i)*al**2.*(cgff2*gluxy+cqff2*
!     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))
!      END IF
!
!  424 CONTINUE
!
!      else
!c$$$      xbmax=1./(1.+ 4.)
!      xbmax=1.D0/(1.D0+ 4.D0)
!      eps=1.
!      if(xbmax.le.x) go to 421
!      DO 425 I=1,NTERMS
!      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
!      XY=X/Y
!      CALL FETCH(XY,Sbot,IPN,FTEMP)
!      gluxy=ftemp(2)
!      IF (IORD.EQ.0) THEN
!      cg21b=2.*facb*cheavy(1,y,eps)
!      ffb=ffb+0.5*(xbmax-x)*wi(i)*albot*(cg21b)*gluxy
!      ELSE
!c$$$      singxy=ftemp(5)+9./8.*ftemp(3) ! Why 9/8 and not 9/4?
!      singxy=ftemp(5)           ! G.W. 27/07/2007
!      cgff2=facb*(c2gffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .c2gffnsh(y,eps)*0.5*exp(1-eps**2))
!      cqff2=facb*(c2qffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .c2qffnsh(y,eps)*0.5*exp(1-eps**2))
!c$$$      ffb=ffb+0.5*(xbmax-x)*wi(i)*albot**2.*(cgff2*gluxy+cqff2*singxy)
!C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
!      ffb=ffb+0.5*(xbmax-x)*wi(i)*albot**2.*(cgff2*gluxy+cqff2*
!     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))
!      END IF
!
!  425 CONTINUE
!      endif
!
!      IF (IORD.NE.0) THEN
!
!C--   G.W. 05/11/2007 This contribution was missing for b, but included for c.
!      xbmax=1./(1.+epsb4)
!      if(xbmax.le.x) go to 421
!      xbmup=x/xbmax
!      CALL FETCH(XBMUP,S,IPN,FTEMP)
!      if(epsb.gt.1.) ftemp(4)=0. ! G.W. 05/11/2007
!      fflb=fflb+ftemp(4)*(AL**2*CLNN2C(xbmup))
!     &     *1.25*(1/(1+4.*epsb)-0.2) ! G.W. 05/11/2007
!
! 629  continue
!      if(epsb.gt.1.) then
!      xbmax=1./(1.+epsb4)
!      eps=epsb
!      if(xbmax.le.x) go to 421
!      DO 624 I=1,NTERMS
!      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X)
!      XY=X/Y
!      CALL FETCH(XY,S,IPN,FTEMP)
!      gluxy=ftemp(2)
!      singxy=ftemp(5)
!      cgffl=facb*(clgffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .clgffnsh(y,eps)*0.5*exp(1-eps**2))
!      cqffl=facb*(clqffnsl(y,eps)*(1-0.5*exp(1-eps**2))+
!     .clqffnsh(y,eps)*0.5*exp(1-eps**2))
!c$$$      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*singxy)
!C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
!      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgffl*gluxy+cqffl*
!     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))
!
!  624 CONTINUE
!      else
!      xbmax=(1./(1.+epsb4))
!      eps=epsb
!      if(xbmax.le.x) go to 421
!      DO 625 I=1,NTERMS
!      Y=0.5*(xbmax-X)*XI(I)+0.5*(xbmax+X) 
!      XY=X/Y
!      CALL FETCH(XY,S,IPN,FTEMP)
!      if(epsb.gt.1.) ftemp(3)=0.
!      gluxy=ftemp(2)
!      fbxy=ftemp(4)
!      singxy=ftemp(5)
!      ymul=y*(1+epsb4)
!      Y1mul=1.-Ymul
!      DL=LOG(Ymul)
!      DL2=DL*DL
!      DLM1=LOG(Y1mul)
!      DLM2=DLM1*DLM1
!      DLM3=DLM2*DLM1
!      DLM4=DLM3*DLM1
!c$$$      FNS2LQmul=128./9.d0*ymul*DLM1**2-46.50*ymul*DLM1-84.094*DL*DLM1
!c$$$     x-37.338 +89.53*ymul
!c$$$     X+33.82*ymul**2+ymul*DL*(32.90+18.41*DL)-128./9.d0*DL
!c$$$     X+16./27.d0*enf*(6.*ymul*DLM1-12*ymul*DL-25.*ymul+6.)
!      FNS2LQmul = CLNN2A(YMUL,INT(ENF)) ! G.W. 02/11/2007
!      FS2LQmul=((15.94-5.212*ymul)*Y1mul*Y1mul*DLM1+(0.421+1.520*ymul)
!     x*DL*DL+28.09*Y1mul*DL-(2.370/Ymul-19.27)*Y1mul**3)
!      cgvfl=facb*((clgffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
!     .clgffnsl(y,eps)*0.5*exp(1-1/eps**2))-clgvfsub(ymul,eps)/xcmax)
!      cqvfl=facb*((clqffnsh(y,eps)*(1-0.5*exp(1-1/eps**2))+
!     .clqffnsl(y,eps)*0.5*exp(1-1/eps**2)))
!      fflb=fflb+0.5/xbmax*(xbmax-x)*wi(i)*al**2.*(fbxy*FNS2LQmul+
!     xfbxy*FS2LQmul)*1.25*(1/(1+4.*eps)-0.2)
!c$$$      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*singxy)
!C--   G.W. 25/07/2007 Add charm and bottom contributions to singlet.
!      fflb=fflb+0.5*(xbmax-x)*wi(i)*al**2.*(cgvfl*gluxy+cqvfl*
!     &     ( singxy + 9./4.*ftemp(3) + 9.*ftemp(4) ))
!
!  625 CONTINUE
!      endif
!      END IF
!
!  421 CONTINUE
!
!      if(ffb.lt.0.) ffb=0.
!
!      F2 = ffx+ffc+ffb
!      F2C = ffc
!      F2B = ffb
!      FL = fflx+fflc+fflb
!      FLC = fflc
!      FLB = fflb

      RETURN
      END
C--   End of MSTWNC subroutine.

