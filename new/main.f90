program main
implicit none
double precision  X,q,f2,f2c,f2b,fl,flc,flb,q2
double precision mc,mt,mb,asmur,mur,fr2
integer ipn,i,iordmstw,j,isetjam
integer dim_f2
parameter (dim_f2=52)
double precision xdata(dim_f2),Q2data(dim_f2)
double precision rs(dim_f2),sigr(dim_f2),norm(dim_f2)
double precision stat(dim_f2),cor(dim_f2),ucor(dim_f2)
double precision f2cjam,flcjam,thy,y,yp,m2,s
COMMON/iordCommon/iordmstw
COMMON/setjam/isetjam
INTEGER alphaSorder,alphaSnfmax
DOUBLE PRECISION distance,tolerance, mCharm,mBottom,alphaSQ0,alphaSMZ
COMMON/mstwCommon/distance,tolerance,mc,mb,ASMUR,alphaSMZ,alphaSorder,alphaSnfmax
CHARACTER prefix*50,prefix1*50,cl*4,namejam*60
integer n
double precision alphas

iordmstw = 1 !--0=LO, 1=NLO, 2=NNLO
ipn=1        !--1 = proton, 2 = neutron, 3 = nucleus DSSZ, 4 = nucleus EPS09 
isetjam = 0  !--replica 0

FR2 = 1.D0                !--ratio of mu_f^2 to mu_r^
MUR = 1.D0                !--input mu_r in G
ASMUR = 0.49128D0         !--input value of alpha_s at mu_r
MC = 1.4D0                !--charm quark mass
MB = 4.75D0               !--bottom quark mass
MT = 175d10               !--top quark mas          

CALL WATE96
CALL INITALPHAS(IORDmstw, FR2, mur, ASMUR, MC, MB, MT)

namejam = 'JAM19PDF_proton_nlo'
call initpdfsetbyname(namejam)
!call numberpdf(n) !number of replicas
call initpdf(isetjam)

x  = 0.5d0
q2 = 10d0
q  = dsqrt(q2)
call MSTWNC(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
write(*,'(F20.10)') ALPHAS(q)
write(*,'(F20.10)') f2
write(*,'(F20.10)') f2c
write(*,'(F20.10)') f2b
write(*,'(F20.10)') fl
write(*,'(F20.10)') flc
write(*,'(F20.10)') flb
end program






























    
