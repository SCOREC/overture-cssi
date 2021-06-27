# CgMp examples -- cocentric cylinders CHT+INS+AD

ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp


#  Partitioned implicit-time-stepping + interface iteration
cocentricCylinders: cocentricCylindersGridi2p.order2.hdf
	$(cgmp) cgmp -noplot cocentricCylinders.cmd -g=cocentricCylindersGridi2p.order2.hdf -tf=1 -tp=.1 -ts=im -nu=.05 -kSolid=1. -kappaSolid=1. -kFluid=.2 -kappaFluid=.2 -solver=best -coupled=0 -nc=10 -iTol=1.e-3 -iOmega=.67 -debug=3 -show=cocentricCyl.show -go=go > cocentricCyl.out

# grid 
cocentricCylindersGridi2p.order2.hdf:
	$(ogen) -noplot cocentricCylindersGrid -zPeriodic=1 -factor=2

