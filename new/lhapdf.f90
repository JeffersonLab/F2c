subroutine setup(fname,iset)
integer iset
character(len=255) fname
call initpdfsetbyname(fname)
call initpdf(iset)
!call numberpdf(n) !number of replicas
end subroutine



