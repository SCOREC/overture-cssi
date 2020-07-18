# ogen = /home/henshw/Overture.g/bin/ogen
ogen = $(Overture)/bin/ogen

cgmx = $(CGBUILDPREFIX)/mx/bin/cgmx

slabsGrids: slabsGrid1 slabsGrid1p slabsGrid3d1 slabsGrid2 slabsGrid3 slabsGrid2p

# 1 slab (3 domains)
slabsGrid1:
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -order=2 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -order=2 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -order=2 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -order=2 -interp=e -factor=16
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -order=4 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -order=4 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -order=4 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -order=4 -interp=e -factor=16

# 1 slab - periodic 
slabsGrid1p:
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -periodic=p -order=2 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -periodic=p -order=2 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -periodic=p -order=2 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -periodic=p -order=2 -interp=e -factor=16
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -periodic=p -order=4 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -periodic=p -order=4 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -periodic=p -order=4 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=1 -x1=-.5 -x2=.5 -periodic=p -order=4 -interp=e -factor=16

# 1 slab in 3D 
slabsGrid3d1:
	$(ogen) -noplot slabsGrid3d -numSlabs=1 -x1=-.5 -x2=.5 -order=2 -interp=e -factor=2
	$(ogen) -noplot slabsGrid3d -numSlabs=1 -x1=-.5 -x2=.5 -order=2 -interp=e -factor=4
	$(ogen) -noplot slabsGrid3d -numSlabs=1 -x1=-.5 -x2=.5 -order=2 -interp=e -factor=8
	$(ogen) -noplot slabsGrid3d -numSlabs=1 -x1=-.5 -x2=.5 -order=2 -interp=e -factor=16
	$(ogen) -noplot slabsGrid3d -numSlabs=1 -x1=-.5 -x2=.5 -order=4 -interp=e -factor=2
	$(ogen) -noplot slabsGrid3d -numSlabs=1 -x1=-.5 -x2=.5 -order=4 -interp=e -factor=4
	$(ogen) -noplot slabsGrid3d -numSlabs=1 -x1=-.5 -x2=.5 -order=4 -interp=e -factor=8
	$(ogen) -noplot slabsGrid3d -numSlabs=1 -x1=-.5 -x2=.5 -order=4 -interp=e -factor=16



# 2 slabs (4 domains) 
slabsGrid2:
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -order=2 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -order=2 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -order=2 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -order=2 -interp=e -factor=16
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -order=4 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -order=4 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -order=4 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -order=4 -interp=e -factor=16

# 2 slabs - periodic 
slabsGrid2p:
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -periodic=p -order=2 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -periodic=p -order=2 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -periodic=p -order=2 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -periodic=p -order=2 -interp=e -factor=16
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -periodic=p -order=4 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -periodic=p -order=4 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -periodic=p -order=4 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=2 -x1=-.5 -x2=.0 -x3=.5 -periodic=p -order=4 -interp=e -factor=16


# 3 slabs (5 domains) 
slabsGrid3:
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -x1=-.5 -x2=.0 -x3=.5 -x4=1 -order=2 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -x1=-.5 -x2=.0 -x3=.5 -x4=1 -order=2 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -x1=-.5 -x2=.0 -x3=.5 -x4=1 -order=2 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -x1=-.5 -x2=.0 -x3=.5 -x4=1 -order=2 -interp=e -factor=16
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -x1=-.5 -x2=.0 -x3=.5 -x4=1 -order=4 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -x1=-.5 -x2=.0 -x3=.5 -x4=1 -order=4 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -x1=-.5 -x2=.0 -x3=.5 -x4=1 -order=4 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -x1=-.5 -x2=.0 -x3=.5 -x4=1 -order=4 -interp=e -factor=16

# 3 slabs (5 domains) periodic  total width=1
sixth = .1666666666666666666667
slabsGrid3p:
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -ya=-.1 -yb=.1 -x1=-.5 -x2=-$(sixth) -x3=$(sixth) -x4=.5 -periodic=p -order=2 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -ya=-.1 -yb=.1 -x1=-.5 -x2=-$(sixth) -x3=$(sixth) -x4=.5 -periodic=p -order=2 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -ya=-.1 -yb=.1 -x1=-.5 -x2=-$(sixth) -x3=$(sixth) -x4=.5 -periodic=p -order=2 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -ya=-.1 -yb=.1 -x1=-.5 -x2=-$(sixth) -x3=$(sixth) -x4=.5 -periodic=p -order=2 -interp=e -factor=16
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -ya=-.1 -yb=.1 -x1=-.5 -x2=-$(sixth) -x3=$(sixth) -x4=.5 -periodic=p -order=4 -interp=e -factor=2
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -ya=-.1 -yb=.1 -x1=-.5 -x2=-$(sixth) -x3=$(sixth) -x4=.5 -periodic=p -order=4 -interp=e -factor=4
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -ya=-.1 -yb=.1 -x1=-.5 -x2=-$(sixth) -x3=$(sixth) -x4=.5 -periodic=p -order=4 -interp=e -factor=8
	$(ogen) -noplot slabsGrid2d -numSlabs=3 -ya=-.1 -yb=.1 -x1=-.5 -x2=-$(sixth) -x3=$(sixth) -x4=.5 -periodic=p -order=4 -interp=e -factor=16



# new block grid, width=1
blockGrid:
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -ya=-1 -yb=1 -interp=e -order=2 -factor=2
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -ya=-1 -yb=1 -interp=e -order=2 -factor=4
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -ya=-1 -yb=1 -interp=e -order=2 -factor=8
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -ya=-1 -yb=1 -interp=e -order=4 -factor=2 
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -ya=-1 -yb=1 -interp=e -order=4 -factor=4 
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -ya=-1 -yb=1 -interp=e -order=4 -factor=8 
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -ya=-1 -yb=1 -interp=e -order=4 -factor=2 -numGhost=3
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -ya=-1 -yb=1 -interp=e -order=4 -factor=4 -numGhost=3
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -ya=-1 -yb=1 -interp=e -order=4 -factor=8 -numGhost=3
