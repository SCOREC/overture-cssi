ogen = ${Overture}/bin/ogen
conv = ${CG}/common/src/conv.p


grids:
	$(ogen) freeSurfaceCylinder.cmd -factor=1
	$(ogen) freeSurfaceCylinder.cmd -factor=2
	$(ogen) freeSurfaceCylinder.cmd -factor=4
	$(ogen) freeSurfaceCylinder.cmd -factor=8



csConv: # need more time steps for accuracy
	# $(conv) ../../conv/surfaceTension.conv -option=cs -ts=pc -tf=0.2 -psolver=best -ng=4 -caseNumber=1 -amp=0.000001 -ng=3 -cfl=0.1
	# $(conv) ../../conv/surfaceTension.conv -option=cs -ts=pc -tf=0.2 -psolver=best -ng=4 -caseNumber=2 -amp=0.000001 -ng=3 -cfl=0.1
	# $(conv) ../../conv/surfaceTension.conv -option=cs -ts=pc -tf=0.2 -psolver=best -ng=4 -caseNumber=3 -amp=0.000001 -ng=3 -cfl=0.1
	$(conv) ../../conv/surfaceTension.conv -option=cs -ts=pc -tf=0.2 -psolver=best -ng=4 -caseNumber=4 -amp=0.000001 -ng=3 -cfl=0.1
	$(conv) ../../conv/surfaceTension.conv -option=cs -ts=pc -tf=0.2 -psolver=best -ng=4 -caseNumber=5 -amp=0.000001 -ng=3 -cfl=0.1
	$(conv) ../../conv/surfaceTension.conv -option=cs -ts=pc -tf=0.2 -psolver=best -ng=4 -caseNumber=6 -amp=0.000001 -ng=3 -cfl=0.1

ins:
	cd ../../; \
	make
