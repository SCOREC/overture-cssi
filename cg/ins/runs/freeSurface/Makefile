ogen = ${Overture}/bin/ogen
conv = ${CG}/common/src/conv.p

doc:
	pdflatex freeSurface

grids:
	$(ogen) -noplot freeSurfaceGrid2d -interp=e -periodic=p -factor=1 -amp=0
	$(ogen) -noplot freeSurfaceGrid2d -interp=e -periodic=p -factor=2 -amp=0
	$(ogen) -noplot freeSurfaceGrid2d -interp=e -periodic=p -factor=4 -amp=0
	$(ogen) -noplot freeSurfaceGrid2d -interp=e -periodic=p -factor=8 -amp=0
	$(ogen) -noplot freeSurfaceGrid2d -interp=e -periodic=p -factor=16 -amp=0

ins:
	cd ../../; \
	make

cfConv:
	# $(conv) ../../conv/surfaceTension.conv -option=cf -ts=pc -tf=.2 -tp=.01 -psolver=yale -cg=squarep -ng=4 -caseNumber=1 -amp=.00001 -cfl=.00001
	# $(conv) ../../conv/surfaceTension.conv -option=cf -ts=pc -tf=.2 -tp=.01 -psolver=yale -cg=squarep -ng=4 -caseNumber=2 -amp=.00001 -cfl=.00001
	# $(conv) ../../conv/surfaceTension.conv -option=cf -ts=pc -tf=.2 -tp=.01 -psolver=yale -cg=squarep -ng=4 -caseNumber=3 -amp=.00001 -cfl=.00001
	$(conv) ../../conv/surfaceTension.conv -option=cf -ts=pc -tf=.2 -tp=.01 -psolver=yale -cg=squarep -ng=4 -caseNumber=5 -amp=.00001 -cfl=.1
	$(conv) ../../conv/surfaceTension.conv -option=cf -ts=pc -tf=.2 -tp=.01 -psolver=yale -cg=squarep -ng=4 -caseNumber=6 -amp=.00001 -cfl=.1
	$(conv) ../../conv/surfaceTension.conv -option=cf -ts=pc -tf=.2 -tp=.01 -psolver=yale -cg=squarep -ng=4 -caseNumber=7 -amp=.00001 -cfl=.1
