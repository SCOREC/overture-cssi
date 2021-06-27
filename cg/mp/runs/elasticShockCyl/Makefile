all : linkFiles diskDeform

linkFiles:
	ln -sf $(Overture)/sampleGrids/diskDeform.cmd .

ogen = $(Overture)/bin/ogen
ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp

diskDeform:
	$(ogen) -noplot diskDeform -factor=1 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="diskDeformBig1e.hdf"
	$(ogen) -noplot diskDeform -factor=2 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="diskDeformBig2e.hdf"
	$(ogen) -noplot diskDeform -factor=4 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="diskDeformBig4e.hdf"
	$(ogen) -noplot diskDeform -factor=8 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="diskDeformBig8e.hdf"
	$(ogen) -noplot diskDeform -factor=16 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="diskDeformBig16e.hdf"


# Shock hitting two elastic disks:
# coarse grid run:
shockTwoElasticDisks: multiDiskDeformGridi2.hdf
	$(cgmp) shockMultiDisk -method=cns -cnsVariation=godunov -smVariation=g -g=multiDiskDeformGridi2 -scf=5.0 -tp=.1 -tf=2.0 -nc=1 -fic=shock -pOffset=.71440 -lambdaSolid=4.0 -muSolid=2.0 -diss=5.0 -cnsGodunovOrder=1 -debug=1 -show=shockTwoElasticDisks.show -go=go > shockTwoElasticDisks.out

# fine grid:
shockTwoElasticDisksFine: multiDiskDeformGride8.hdf
	$(cgmp) -noplot shockMultiDisk -method=cns -cnsVariation=godunov -smVariation=g -g=multiDiskDeformGride8 -scf=5.0 -tp=.1 -tf=2.0 -nc=1 -fic=shock -pOffset=.71440 -lambdaSolid=4.0 -muSolid=2.0 -diss=5.0 -cnsGodunovOrder=1 -debug=1 -show=shockTwoElasticDisksG8.show -go=go > shockTwoElasticDisksG8.out
#
# grids
multiDiskDeformGridi2.hdf:
		$(ogen) -noplot multiDiskDeformGrid -factor=2 -width=.4

multiDiskDeformGride8.hdf:		
		$(ogen) -noplot multiDiskDeformGrid -interp=e -factor=8