subroutine setup(fname,iset)
integer iset
character(len=255) fname
call initpdfsetbyname(fname)
call initpdf(iset)
!call numberpdf(n) !number of replicas
end subroutine


subroutine setup_pdf(fname)
character(len=255) fname
call initpdfsetbyname(fname)
end subroutine


subroutine setup_pdfset(iset)
integer iset
call initpdf(iset)
end subroutine


