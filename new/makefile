#lhapdf=-L/home/nsato/apps/lhapdf/lib  -lLHAPDF
lhapdf=-L/work/JAM/apps/lhapdf/lib  -lLHAPDF

compile:
	gfortran -o pro main.f90 alphaS.f hqlib.f  $(lhapdf) 

pylib:
	rm -f *.so
	f2py -c \
	-m hqlib \
	--fcompiler=gnu95 \
	--f90flags="-Wtabs" \
	$(lhapdf) \
	lhapdf.f90 \
	alphaS.f \
	hqlib.f
	#mv hqlib*so hqlib.so


clean:
	rm *.so


