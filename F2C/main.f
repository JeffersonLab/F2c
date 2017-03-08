      PROGRAM SFDAT
      IMPLICIT DOUBLE PRECISION (A - H, O - Z)
      COMMON /PDFI/ IPDF
      integer ipdf
      CHARACTER(len=255) PATH,root
      common /root/root
      root='./'
      X    = 0.01d0
      Q2   = 5d0
      IPDF = 1
      CALL SFH(F2H, X, Q2,1)
      print*,F2H
      END
