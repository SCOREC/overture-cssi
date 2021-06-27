# Make deforming ellipse grids


ogen = $(Overture)/bin/ogen
cgmp = $(CGBUILDPREFIX)/mp/bin/cgmp


grids: deformingEllipseGrid

deformingEllipseGrid:
	$(ogen) -noplot deformingEllipseGrid -factor=2 -interp=e -width=.25
	$(ogen) -noplot deformingEllipseGrid -factor=4 -interp=e -width=.25
	$(ogen) -noplot deformingEllipseGrid -factor=8 -interp=e -width=.25
	$(ogen) -noplot deformingEllipseGrid -factor=16 -interp=e -width=.25

# Shock in a compressible fluid hitting a solid disk (Neo-Hookean solid)
# Coarse grid run:
neoHookean: deformingEllipseGride4.hdf
	$(cgmp) -noplot shockCyl -method=cns -cnsVariation=godunov -smVariation=g  -g=deformingEllipseGride4 -scf=1.0 -tp=.1 -tf=1.0 -nc=1 -d1=outerDomain -d2=innerDomain -fic=shock -a1=1.183215957 -Mshock=1.3 -pOffset=1.0 -lambdaSolid=4.0 -muSolid=4. -rhoSolid=.1  -diss=2.0 -tangentialDissipationSolid=1.0 -cfl=.9 -cnsGodunovOrder=2 -debug=1 -xShock=-2.2500001 -a1= -flushFrequency=2 -godunovType=4 -pi=1 -piGhostOption=3 -multiDomainAlgorithm=1 -show=shockCylNeoHookean.show -go=go > shockCylNeoHookean.out 

deformingEllipseGride4.hdf:
	$(ogen) -noplot deformingEllipseGrid -factor=4 -interp=e -width=.25

# Shock in a compressible fluid hitting a solid disk (Neo-Hookean solid)
# Finer grid run:
neoHookeanG8: deformingEllipseGride8.hdf
	$(cgmp) -noplot shockCyl -method=cns -cnsVariation=godunov -smVariation=g  -g=deformingEllipseGride8 -scf=1.0 -tp=.1 -tf=1.0 -nc=1 -d1=outerDomain -d2=innerDomain -fic=shock -a1=1.183215957 -Mshock=1.3 -pOffset=1.0 -lambdaSolid=4.0 -muSolid=4. -rhoSolid=.1  -diss=2.0 -tangentialDissipationSolid=1.0 -cfl=.9 -cnsGodunovOrder=2 -debug=1 -xShock=-2.2500001 -a1= -flushFrequency=2 -godunovType=4 -pi=1 -piGhostOption=3 -multiDomainAlgorithm=1 -show=shockCylNeoHookeanG8.show -go=go > shockCylNeoHookeanG8.out 

# grids 
deformingEllipseGride4.hdf:
	$(ogen) -noplot deformingEllipseGrid -factor=4 -interp=e -width=.25

deformingEllipseGride8.hdf:
	$(ogen) -noplot deformingEllipseGrid -factor=8 -interp=e -width=.25