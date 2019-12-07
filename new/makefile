lhapdf=-L/home/nsato/apps/lhapdf/lib  -l LHAPDF

compile:
	gfortran -o pro main.f90 alphaS.f f2clib-bkp.f $(lhapdf) 

pylib:
	rm -f *.so
	f2py -c \
	-m f2clib \
	--fcompiler=gnu95 \
	--f90flags="-Wtabs" \
	alphaS.f \
	f2clib.f
	mv f2clib*so f2clib.so



clean:
	rm *.so


