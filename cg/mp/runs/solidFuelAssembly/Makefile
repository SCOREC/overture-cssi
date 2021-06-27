# CgMp examples -- solid fuel assembly

ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp

# ------- Fuel assembly with heated fuel rods ------
#  Partitioned implicit time-stepping (predictor corrector):
#  Here is short run with one rod
fuelAssemblyOneRod: heatExchangerGride2.order2.hdf
	$(cgmp) -noplot solidFuelAssembly.cmd -g=sfaGridOneRodHeight1e2.order2.hdf -solver=best -ts=imp -coupled=0 -nc=8 -dtMax=.05 -tp=.1 -tf=1 -show=sfaOneRod.show -go=go > sfaOneRod.out 
# grid for one rod
sfaGridOneRodHeight1e2.order2.hdf:
	$(ogen) -noplot solidFuelAssemblyGrid.cmd -prefix=sfaOneRodHeight1Grid -npins=1 -zb=1 -interp=e -factor=2


#  Here is short run with five rods
fuelAssembly5Rods: sfa5RodsGride1.order2.hdf
	$(cgmp) -noplot solidFuelAssembly.cmd -g=sfa5RodsGride1.order2.hdf -solver=best -ts=imp -coupled=0 -nc=8 -dtMax=.05 -tp=.1 -tf=1 -show=sfa5Rods.show -go=go > sfa5Rods.out 
# grid for 5 rods
sfa5RodsGride1.order2.hdf:
	$(ogen) -noplot solidFuelAssemblyGrid.cmd -prefix=sfa5RodsGrid -npins=3 -zb=2 -interp=e -factor=1
