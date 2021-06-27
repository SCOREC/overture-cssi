# CgMp examples -- curved pipe CHT+INS+AD

ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp


#  Explicit time-stepping + coupled interface iteration
curvedPipe: curvedPipeGrid90e1.order2.hdf
	$(cgmp) -noplot curvedPipe.cmd -g=curvedPipeGrid90e1.order2.hdf -solver=best -tf=1 -tp=.1 -show=curvedPipe.show -go=go

#  Partitioned implicit-time-stepping + interface iteration
curvedPipeImplicit: curvedPipeGrid90e2.order2.hdf
	$(cgmp) -noplot curvedPipe.cmd -g=curvedPipeGrid90e2.order2.hdf -tf=1 -tp=.1 -ts=im -nu=.1 -kappa=.025 -solver=best -coupled=0 -nc=15 -show=curvedPipeImp.show -go=go

# grids
curvedPipeGrid90e1.order2.hdf:
	$(ogen) -noplot curvedPipeGrid -angle=90 -interp=e -factor=1

curvedPipeGrid90e2.order2.hdf:
	$(ogen) -noplot curvedPipeGrid -angle=90 -interp=e -factor=2

