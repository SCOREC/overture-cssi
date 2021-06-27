#
# CgMP: rotating solid disk (non-linear solid) in an annular compressible fluid
#

ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp


# Rotating SVK (St. Venant Kirchoff nonlinear solid) disk
rotatingDiskSVK: diskDeformFixedBig4e.hdf
	$(cgmp) -noplot rotatingDisk.cmd -g=diskDeformFixedBig4e -method=cns -cfl=.75 -tp=.1 -tf=1.0 -tz=none -problem=1 -cnsVariation=godunov -cnsGodunovOrder=2 -smVariation=g -pi=1 -en=max -adCns=10.0 -bcOption=4 -rotationRate=.1 -piGhostOption=3 -piOption=0 -godunovType=2 -stressRelaxation=1 -d1=outerDomain -d2=innerDomain -rho0=14 -pOffset=10 -lambdaSolid=1.0 -muSolid=1.0 -diss=1.0 -multiDomainAlgorithm=1 -show=rotatingDiskSVKG4.show -go=go >! rotatingDiskSVKG4.show


# grids
diskDeformFixedBig4e.hdf:
	$(ogen) -noplot diskDeform -factor=4 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name=diskDeformBig4e.hdf

