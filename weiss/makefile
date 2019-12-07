compile:
	gfortran -o pro \
	grv98.f  \
	EPS09.f  \
	slmach.f \
	dqags.f  \
	hq_usr_pdf.f \
	hq_core_sfh.f \
	main.f

pycode:
	f2py -c --fcompiler=gnu95 -m fF2C  \
	grv98.f  \
	EPS09.f  \
	slmach.f \
	dqags.f  \
	hq_usr_pdf.f \
	hq_core_sfh.f 


