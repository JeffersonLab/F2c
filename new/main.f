      IMPLICIT NONE

      double precision  X,q,f2,f2c,f2b,fl,flc,flb,q2
      double precision mc,mt,mb,asmur,mur,fr2
      integer ipn,i,iordmstw,j,isetjam
      integer dim_f2
      parameter (dim_f2=52)
      double precision xdata(dim_f2),Q2data(dim_f2),
     1  rs(dim_f2),sigr(dim_f2),norm(dim_f2),
     1  stat(dim_f2),cor(dim_f2),ucor(dim_f2)
      double precision f2cjam,flcjam,thy,y,yp,m2,s
      COMMON/iordCommon/iordmstw
      COMMON/setjam/isetjam
      INTEGER alphaSorder,alphaSnfmax
      DOUBLE PRECISION distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ
      COMMON/mstwCommon/distance,tolerance,
     &     mc,mb,ASMUR,alphaSMZ,alphaSorder,alphaSnfmax

            FR2 = 1.D0                ! ratio of mu_f^2 to mu_r^
            MUR = 1.D0                ! input mu_r in G
            ASMUR = 0.49128D0         ! input value of alpha_s at mu_r
            MC = 1.4D0                ! charm quark mass
            MB = 4.75D0               ! bottom quark mass
            MT = 175d10               ! top quark mas          

        iordmstw = 1 !0=LO, 1=NLO, 2=NNLO
        CALL WATE96
        CALL INITALPHAS(IORDmstw, FR2, mur, ASMUR, MC, MB, MT)
c      ipn=1 ! 1 = proton, 2 = neutron, 3 = nucleus DSSZ, 4 = nucleus EPS09 
      

c everything
      open(20, file='data_10037.txt',status='unknown')
      do i=1,dim_f2,1
       read(20,*) Xdata(i),Q2data(i),rs(i),sigr(i),stat(i),ucor(i)
     1 ,cor(i),norm(i)
      enddo
      close(20)   
      
      m2=0.8815688 !GeV**2 proton mass**2   
c     F2c proton (JAM19)     
      ipn=1
    
      isetjam = 0 !replica 0
      open(20,file= 
     1 'test-mstwcode-jam19-sigmaredc.dat'
     2 ,status='new')
      do i=1,dim_f2,1
         x = Xdata(i)
         q2 = Q2data(i)
         s = rs(i)**2.d0
         q = dsqrt(q2)
         y = q2 /x / (s-M2)
         yp = 1.d0+(1.d0-y)**2.d0
         call MSTWNC(x,q,ipn,f2,f2c,f2b,fl,flc,flb)
         f2cjam = f2c
         flcjam = flc
         thy = f2cjam - y**2.d0/yp*flcjam
c         print*, x, q2, f2cjam, flcjam, thy
         write(20,*) x, q2, f2cjam, flcjam, thy
      enddo
      close(20)
      end program

     
