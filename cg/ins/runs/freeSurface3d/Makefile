ogen = ${Overture}/bin/ogen

grids:
	$(ogen) -noplot freeSurfaceGrid3d -interp=e -factor=.25 -bc=p
	$(ogen) -noplot freeSurfaceGrid3d -interp=e -factor=1 -bc=p
	$(ogen) -noplot freeSurfaceGrid3d -interp=e -factor=2 -bc=p
	$(ogen) -noplot freeSurfaceGrid3d -interp=e -factor=4 -bc=p
	$(ogen) -noplot freeSurfaceGrid3d -interp=e -factor=8 -bc=p


ins:
	cd ../../; \
	make
