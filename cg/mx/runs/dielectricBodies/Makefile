#
#   Makefile for dielectric bodies
#

# ogen = /home/henshw/Overture.g/bin/ogen
ogen = $(Overture)/bin/ogen

cgmx = $(CGBUILDPREFIX)/mx/bin/cgmx

# Make the documentation:
doc:
	pdflatex dielectricBodies
	bibtex dielectricBodies

# generate some grids
grids : grids2d grids3d

# 2D grids: 
grids2d: disk ellipse rod cross rectangle dielectricBlockGrid2d multiblock 

# 3D grids: 
grids3d: sphere pillbox


# one disk for MLA self convergence
oneDiskGrids:
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=2 -periodic=np -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=2 -periodic=np -factor=4
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=2 -periodic=np -factor=8
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=2 -periodic=np -factor=16
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=4 -numGhost=3 -periodic=np -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=4 -numGhost=3 -periodic=np -factor=4
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=4 -numGhost=3 -periodic=np -factor=8
	$(ogen) -noplot diskArrayGrid -prefix=oneDiskMLAGrid -deltaRadius=.1 -xa=-1.0 -xb=1.0 -ya=-1.0 -yb=1.0 -rad=.5 -nCylx=1 -nCyly=1 -deltaX=1.25 -deltaY=1.25 -order=4 -numGhost=3 -periodic=np -factor=16

# make -j6 oneDiskRuns 
oneDiskRuns: oneDiskRunsG2O2 oneDiskRunsG4O2 oneDiskRunsG8O2 oneDiskRunsG2O4 oneDiskRunsG4O4 oneDiskRunsG8O4

oneDiskRunsG2O2:
	$(cgmx) -noplot mxScript -g=oneDiskMLAGridFixednpe2.order4.ng3 -dm=gdm -nm=mla -numBodies=1 -materialFile=disk4MlaMaterials.h -tf=2 -tp=.5 -useSosupDissipation=1 -ic=gp -kx=1 -ay=1 -az=0 -x0=-1.5 -rbc=rbcNonLocal -background=backGround -debug=0 -show=oneDiskMLAO4G2.show -go=go > oneDiskMLAO4G2.out
oneDiskRunsG4O2:
	$(cgmx) -noplot mxScript -g=oneDiskMLAGridFixednpe4.order4.ng3 -dm=gdm -nm=mla -numBodies=1 -materialFile=disk4MlaMaterials.h -tf=2 -tp=.5 -useSosupDissipation=1 -ic=gp -kx=1 -ay=1 -az=0 -x0=-1.5 -rbc=rbcNonLocal -background=backGround -debug=0 -show=oneDiskMLAO4G4.show -go=go > oneDiskMLAO4G4.out
oneDiskRunsG8O2:
	$(cgmx) -noplot mxScript -g=oneDiskMLAGridFixednpe8.order4.ng3 -dm=gdm -nm=mla -numBodies=1 -materialFile=disk4MlaMaterials.h -tf=2 -tp=.5 -useSosupDissipation=1 -ic=gp -kx=1 -ay=1 -az=0 -x0=-1.5 -rbc=rbcNonLocal -background=backGround -debug=0 -show=oneDiskMLAO4G8.show -go=go > oneDiskMLAO4G8.out

oneDiskRunsG2O4:
	$(cgmx) -noplot mxScript -g=oneDiskMLAGridFixednpe2.order4.ng3 -dm=gdm -nm=mla -numBodies=1 -materialFile=disk4MlaMaterials.h -tf=2 -tp=.5 -useSosupDissipation=1 -ic=gp -kx=1 -ay=1 -az=0 -x0=-1.5 -rbc=rbcNonLocal -background=backGround -debug=0 -show=oneDiskMLAO4G2.show -go=go > oneDiskMLAO4G2.out
oneDiskRunsG4O4:
	$(cgmx) -noplot mxScript -g=oneDiskMLAGridFixednpe4.order4.ng3 -dm=gdm -nm=mla -numBodies=1 -materialFile=disk4MlaMaterials.h -tf=2 -tp=.5 -useSosupDissipation=1 -ic=gp -kx=1 -ay=1 -az=0 -x0=-1.5 -rbc=rbcNonLocal -background=backGround -debug=0 -show=oneDiskMLAO4G4.show -go=go > oneDiskMLAO4G4.out
oneDiskRunsG8O4:
	$(cgmx) -noplot mxScript -g=oneDiskMLAGridFixednpe8.order4.ng3 -dm=gdm -nm=mla -numBodies=1 -materialFile=disk4MlaMaterials.h -tf=2 -tp=.5 -useSosupDissipation=1 -ic=gp -kx=1 -ay=1 -az=0 -x0=-1.5 -rbc=rbcNonLocal -background=backGround -debug=0 -show=oneDiskMLAO4G8.show -go=go > oneDiskMLAO4G8.out


# Make a 2D rectangular region with no body (extra ghost needed for sosup)
rectangle:
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-2 -xb=2 -ya=-.75 -yb=.75 -periodic=np -factor=4
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-2 -xb=2 -ya=-.75 -yb=.75 -periodic=np -factor=8
	 $(ogen) -noplot rectangleArg -prefix=rect -order=4 -xa=-2 -xb=2 -ya=-.75 -yb=.75 -periodic=np -factor=4
	 $(ogen) -noplot rectangleArg -prefix=rect -order=4 -xa=-2 -xb=2 -ya=-.75 -yb=.75 -periodic=np -factor=8
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-8 -xb=2 -ya=-.75 -yb=.75 -periodic=np -factor=4
	 $(ogen) -noplot rectangleArg -prefix=airSlabGrid -order=4 -xa=-4 -xb=4 -ya=-2 -yb=2 -periodic=np -factor=2
	 $(ogen) -noplot rectangleArg -prefix=airSlabGrid -order=4 -xa=-4 -xb=4 -ya=-2 -yb=2 -periodic=np -factor=4
	 $(ogen) -noplot rectangleArg -prefix=airSlabGrid -order=4 -xa=-4 -xb=4 -ya=-2 -yb=2 -periodic=np -factor=8
	 $(ogen) -noplot rectangleArg -prefix=channel -order=2 -xa=-3 -xb=3 -ya=-.5 -yb=.5 -periodic=np -factor=2
	 $(ogen) -noplot rectangleArg -prefix=channel -order=2 -xa=-3 -xb=3 -ya=-.5 -yb=.5 -periodic=np -factor=4
	 $(ogen) -noplot rectangleArg -prefix=channel -order=2 -xa=-3 -xb=3 -ya=-.5 -yb=.5 -periodic=np -factor=4
	 $(ogen) -noplot rectangleArg -prefix=noPecDisk -order=4 -xa=-3 -xb=3 -ya=-1.25 -yb=1.25 -periodic=np -factor=2
	 $(ogen) -noplot rectangleArg -prefix=noPecDisk -order=4 -xa=-3 -xb=3 -ya=-1.25 -yb=1.25 -periodic=np -factor=4
	 $(ogen) -noplot rectangleArg -prefix=noPecDisk -order=4 -xa=-3 -xb=3 -ya=-1.25 -yb=1.25 -periodic=np -factor=8


# block grid 3d width=1
blockGrid3d:
	 $(ogen) -noplot dielectricBlockGrid3d -prefix=blockGrid3d -width=1 -interp=e -order=2 -factor=2
	 $(ogen) -noplot dielectricBlockGrid3d -prefix=blockGrid3d -width=1 -interp=e -order=2 -factor=4
	 $(ogen) -noplot dielectricBlockGrid3d -prefix=blockGrid3d -width=1 -interp=e -order=2 -factor=8
	 $(ogen) -noplot dielectricBlockGrid3d -prefix=blockGrid3d -width=1 -interp=e -order=4 -factor=2 
	 $(ogen) -noplot dielectricBlockGrid3d -prefix=blockGrid3d -width=1 -interp=e -order=4 -factor=4 
	 $(ogen) -noplot dielectricBlockGrid3d -prefix=blockGrid3d -width=1 -interp=e -order=4 -factor=8 
	 $(ogen) -noplot dielectricBlockGrid3d -prefix=blockGrid3d -width=1 -interp=e -order=4 -factor=2 -numGhost=3
	 $(ogen) -noplot dielectricBlockGrid3d -prefix=blockGrid3d -width=1 -interp=e -order=4 -factor=4 -numGhost=3
	 $(ogen) -noplot dielectricBlockGrid3d -prefix=blockGrid3d -width=1 -interp=e -order=4 -factor=8 -numGhost=3


# new block grid, width=1
blockGrid:
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -interp=e -order=2 -periodic=np -factor=2
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -interp=e -order=2 -periodic=np -factor=4
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -interp=e -order=2 -periodic=np -factor=8
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -interp=e -order=4 -periodic=np -factor=2 
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -interp=e -order=4 -periodic=np -factor=4 
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -interp=e -order=4 -periodic=np -factor=8 
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -interp=e -order=4 -periodic=np -factor=2 -numGhost=3
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -interp=e -order=4 -periodic=np -factor=4 -numGhost=3
	 $(ogen) -noplot dielectricBlockGrid2d -prefix=blockGrid -width=1 -interp=e -order=4 -periodic=np -factor=8 -numGhost=3

# Old dielectric block grid:
dielectricBlockGrid2d:
	 $(ogen) -noplot dielectricBlockGrid2d -interp=e -order=2 -periodic=np -factor=2
	 $(ogen) -noplot dielectricBlockGrid2d -interp=e -order=2 -periodic=np -factor=4
	 $(ogen) -noplot dielectricBlockGrid2d -interp=e -order=2 -periodic=np -factor=8
	 $(ogen) -noplot dielectricBlockGrid2d -interp=e -order=4 -periodic=np -factor=2 
	 $(ogen) -noplot dielectricBlockGrid2d -interp=e -order=4 -periodic=np -factor=4 
	 $(ogen) -noplot dielectricBlockGrid2d -interp=e -order=4 -periodic=np -factor=8 
	 $(ogen) -noplot dielectricBlockGrid2d -interp=e -order=4 -periodic=np -factor=2 -numGhost=3
	 $(ogen) -noplot dielectricBlockGrid2d -interp=e -order=4 -periodic=np -factor=4 -numGhost=3
	 $(ogen) -noplot dielectricBlockGrid2d -interp=e -order=4 -periodic=np -factor=8 -numGhost=3

# Multiple dielectric blocs
multiblock:
	$(ogen) -noplot dielectricMultiBlockGrid2d -numBlocks=2 -x1=-.25 -x2=.0 -x3=.25 -order=4 -interp=e -factor=2
	$(ogen) -noplot dielectricMultiBlockGrid2d -numBlocks=2 -x1=-.25 -x2=.0 -x3=.25 -order=4 -interp=e -factor=4
	$(ogen) -noplot dielectricMultiBlockGrid2d -numBlocks=2 -x1=-.25 -x2=.0 -x3=.25 -order=4 -interp=e -factor=8
	$(ogen) -noplot dielectricMultiBlockGrid2d -numBlocks=2 -x1=-.25 -x2=.0 -x3=.25 -order=4 -interp=e -factor=16


# curved dielectric block 
curvedBlock:
	$(ogen) -noplot curvedBlockGrid2d -prefix=convexBlockGrid -order=2 -interp=e -width=.5 -deltaLeft=-.075 -deltaRight=.075 -factor=2
	$(ogen) -noplot curvedBlockGrid2d -prefix=convexBlockGrid -order=2 -interp=e -width=.5 -deltaLeft=-.075 -deltaRight=.075 -factor=4
	$(ogen) -noplot curvedBlockGrid2d -prefix=convexBlockGrid -order=2 -interp=e -width=.5 -deltaLeft=-.075 -deltaRight=.075 -factor=8
	$(ogen) -noplot curvedBlockGrid2d -prefix=convexBlockGrid -order=4 -interp=e -width=.5 -deltaLeft=-.075 -deltaRight=.075 -factor=2
	$(ogen) -noplot curvedBlockGrid2d -prefix=convexBlockGrid -order=4 -interp=e -width=.5 -deltaLeft=-.075 -deltaRight=.075 -factor=4
	$(ogen) -noplot curvedBlockGrid2d -prefix=convexBlockGrid -order=4 -interp=e -width=.5 -deltaLeft=-.075 -deltaRight=.075 -factor=8

pecDisk:
	$(ogen) -noplot diskInAChannelGrid -xa=-3 -xb=3. -interp=e -order=2 -factor=2
	$(ogen) -noplot diskInAChannelGrid -xa=-3 -xb=3. -interp=e -order=2 -factor=4
	$(ogen) -noplot diskInAChannelGrid -xa=-3 -xb=3. -interp=e -order=2 -factor=8
	$(ogen) -noplot diskInAChannelGrid -xa=-3 -xb=3. -interp=e -order=4 -factor=2
	$(ogen) -noplot diskInAChannelGrid -xa=-3 -xb=3. -interp=e -order=4 -factor=4
	$(ogen) -noplot diskInAChannelGrid -xa=-3 -xb=3. -interp=e -order=4 -factor=8


disk:
	$(ogen) -noplot io -prefix=diskInBox -order=2 -interp=e -xa=-5.5 -xb=5.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=2
	$(ogen) -noplot io -prefix=diskInBox -order=2 -interp=e -xa=-5.5 -xb=5.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=4
	$(ogen) -noplot io -prefix=diskInBox -order=2 -interp=e -xa=-5.5 -xb=5.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=8
	$(ogen) -noplot io -prefix=diskInBox -order=4 -interp=e -xa=-5.5 -xb=5.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=2 -numGhost=3
	$(ogen) -noplot io -prefix=diskInBox -order=4 -interp=e -xa=-5.5 -xb=5.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=4 -numGhost=3
	$(ogen) -noplot io -prefix=diskInBox -order=4 -interp=e -xa=-5.5 -xb=5.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=8 -numGhost=3

# For BA Homogenization 
solidDisk:
	$(ogen) -noplot io -prefix=solidDisk -order=2 -interp=e -xa=-2.5 -xb=2.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=1
	$(ogen) -noplot io -prefix=solidDisk -order=2 -interp=e -xa=-2.5 -xb=2.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=2
	$(ogen) -noplot io -prefix=solidDisk -order=2 -interp=e -xa=-2.5 -xb=2.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=4
	$(ogen) -noplot io -prefix=solidDisk -order=2 -interp=e -xa=-2.5 -xb=2.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=8
	$(ogen) -noplot io -prefix=solidDisk -order=4 -interp=e -xa=-2.5 -xb=2.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=2 -numGhost=3
	$(ogen) -noplot io -prefix=solidDisk -order=4 -interp=e -xa=-2.5 -xb=2.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=4 -numGhost=3
	$(ogen) -noplot io -prefix=solidDisk -order=4 -interp=e -xa=-2.5 -xb=2.5 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=8 -numGhost=3

# short domain 
shortSolidDiskGrid:
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=2 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=1
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=2 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=2
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=2 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=4
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=2 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=8
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=2 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=16
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=2 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=32
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=2 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=64
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=4 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=2 -numGhost=3
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=4 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=4 -numGhost=3
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=4 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=8 -numGhost=3
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=4 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=16 -numGhost=3
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=4 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=32 -numGhost=3
	$(ogen) -noplot io -prefix=shortSolidDiskGrid -order=4 -interp=e -xa=-0.75 -xb=0.75 -ya=-.75 -yb=.75 -outerBC=yPeriodic -factor=64 -numGhost=3

# Two disks in a vertical configuration
twoDiskGrid:
	$(ogen) -noplot diskArrayGrid -prefix=twoDiskGrid -xa=-.75 -xb=.75 -ya=-1.5 -yb=1.5 -nCylx=1 -nCyly=2 -dist=1.5 -order=2 -periodic=np -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=twoDiskGrid -xa=-.75 -xb=.75 -ya=-1.5 -yb=1.5 -nCylx=1 -nCyly=2 -dist=1.5 -order=2 -periodic=np -factor=4
	$(ogen) -noplot diskArrayGrid -prefix=twoDiskGrid -xa=-.75 -xb=.75 -ya=-1.5 -yb=1.5 -nCylx=1 -nCyly=2 -dist=1.5 -order=4 -periodic=np -numGhost=3 -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=twoDiskGrid -xa=-.75 -xb=.75 -ya=-1.5 -yb=1.5 -nCylx=1 -nCyly=2 -dist=1.5 -order=4 -periodic=np -numGhost=3 -factor=4

# Four disks in a 2 x 2 configuration
fourDiskGrid:
	$(ogen) -noplot diskArrayGrid -prefix=fourDiskGrid -xa=-1.75 -xb=1.75 -ya=-1.5 -yb=1.5 -nCylx=2 -nCyly=2 -deltaX=1.5 -deltaY=1.5 -order=2 -periodic=np -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=fourDiskGrid -xa=-1.75 -xb=1.75 -ya=-1.5 -yb=1.5 -nCylx=2 -nCyly=2 -deltaX=1.5 -deltaY=1.5 -order=2 -periodic=np -factor=4
	$(ogen) -noplot diskArrayGrid -prefix=fourDiskGrid -xa=-1.75 -xb=1.75 -ya=-1.5 -yb=1.5 -nCylx=2 -nCyly=2 -deltaX=1.5 -deltaY=1.5 -order=4 -periodic=np -numGhost=3 -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=fourDiskGrid -xa=-1.75 -xb=1.75 -ya=-1.5 -yb=1.5 -nCylx=2 -nCyly=2 -deltaX=1.5 -deltaY=1.5 -order=4 -periodic=np -numGhost=3 -factor=4

# 16 disks:
16DiskGrid:
	$(ogen) -noplot diskArrayGrid -prefix=16DiskGrid -xa=-3.25 -xb=3.25 -ya=-3.0 -yb=3.0 -nCylx=4 -nCyly=4 -deltaX=1.5 -deltaY=1.5 -order=2 -periodic=np -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=16DiskGrid -xa=-3.25 -xb=3.25 -ya=-3.0 -yb=3.0 -nCylx=4 -nCyly=4 -deltaX=1.5 -deltaY=1.5 -order=2 -periodic=np -factor=4
	$(ogen) -noplot diskArrayGrid -prefix=16DiskGrid -xa=-3.25 -xb=3.25 -ya=-3.0 -yb=3.0 -nCylx=4 -nCyly=4 -deltaX=1.5 -deltaY=1.5 -order=4 -periodic=np -numGhost=3 -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=16DiskGrid -xa=-3.25 -xb=3.25 -ya=-3.0 -yb=3.0 -nCylx=4 -nCyly=4 -deltaX=1.5 -deltaY=1.5 -order=4 -periodic=np -numGhost=3 -factor=4

# 64 disks:
64DiskGrid:
	$(ogen) -noplot diskArrayGrid -prefix=64DiskGrid -xa=-6.25 -xb=6.25 -ya=-6.0 -yb=6.0 -nCylx=8 -nCyly=8 -deltaX=1.5 -deltaY=1.5 -order=2 -periodic=np -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=64DiskGrid -xa=-6.25 -xb=6.25 -ya=-6.0 -yb=6.0 -nCylx=8 -nCyly=8 -deltaX=1.5 -deltaY=1.5 -order=2 -periodic=np -factor=4
	$(ogen) -noplot diskArrayGrid -prefix=64DiskGrid -xa=-6.25 -xb=6.25 -ya=-6.0 -yb=6.0 -nCylx=8 -nCyly=8 -deltaX=1.5 -deltaY=1.5 -order=4 -periodic=np -numGhost=3 -factor=2
	$(ogen) -noplot diskArrayGrid -prefix=64DiskGrid -xa=-6.25 -xb=6.25 -ya=-6.0 -yb=6.0 -nCylx=8 -nCyly=8 -deltaX=1.5 -deltaY=1.5 -order=4 -periodic=np -numGhost=3 -factor=4


ellipse:
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -xa=-5 -xb=5 -interp=e -factor=2
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -xa=-5 -xb=5 -interp=e -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -xa=-5 -xb=5 -interp=e -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -xa=-5 -xb=5 -interp=e -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -xa=-5 -xb=5 -interp=e -order=4 -factor=16


ellipseSmallDomain:
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalSmallGrid -radX=.5 -radY=.45 -xa=-1.5 -xb=1.5 -ya=-.75 -yb=.75 -interp=e -nr0=6 -factor=2
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalSmallGrid -radX=.5 -radY=.45 -xa=-1.5 -xb=1.5 -ya=-.75 -yb=.75 -interp=e -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalSmallGrid -radX=.5 -radY=.45 -xa=-1.5 -xb=1.5 -ya=-.75 -yb=.75 -interp=e -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalSmallGrid -radX=.5 -radY=.45 -xa=-1.5 -xb=1.5 -ya=-.75 -yb=.75 -interp=e -order=4 -factor=4  -numGhost=3
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalSmallGrid -radX=.5 -radY=.45 -xa=-1.5 -xb=1.5 -ya=-.75 -yb=.75 -interp=e -order=4 -factor=8  -numGhost=3
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalSmallGrid -radX=.5 -radY=.45 -xa=-1.5 -xb=1.5 -ya=-.75 -yb=.75 -interp=e -order=4 -factor=16 -numGhost=3

ellipseAngle30:
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -angle=30 -xa=-5 -xb=5 -interp=e -order=2 -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -angle=30 -xa=-5 -xb=5 -interp=e -order=2 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -angle=30 -xa=-5 -xb=5 -interp=e -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -angle=30 -xa=-5 -xb=5 -interp=e -order=4 -factor=16
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -angle=30 -xa=-5 -xb=5 -interp=e -numGhost=3 -order=4 -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -angle=30 -xa=-5 -xb=5 -interp=e -numGhost=3 -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGrid -angle=30 -xa=-5 -xb=5 -interp=e -numGhost=3 -order=4 -factor=16
 
# Long domain L=20 : [-10,10]
ellipseL20Angle30:
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL20 -angle=30 -xa=-10. -xb=10. -interp=e -order=2 -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL20 -angle=30 -xa=-10. -xb=10. -interp=e -order=2 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL20 -angle=30 -xa=-10. -xb=10. -interp=e -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL20 -angle=30 -xa=-10. -xb=10. -interp=e -order=4 -factor=16
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL20 -angle=30 -xa=-10. -xb=10. -interp=e -numGhost=3 -order=4 -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL20 -angle=30 -xa=-10. -xb=10. -interp=e -numGhost=3 -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL20 -angle=30 -xa=-10. -xb=10. -interp=e -numGhost=3 -order=4 -factor=16

# Mediumdomain L=10 : [-5,5]
ellipseL10Angle30:
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL10 -angle=30 -xa=-5 -xb=5 -interp=e -order=2 -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL10 -angle=30 -xa=-5 -xb=5 -interp=e -order=2 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL10 -angle=30 -xa=-5 -xb=5 -interp=e -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL10 -angle=30 -xa=-5 -xb=5 -interp=e -order=4 -factor=16
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL10 -angle=30 -xa=-5 -xb=5 -interp=e -numGhost=3 -order=4 -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL10 -angle=30 -xa=-5 -xb=5 -interp=e -numGhost=3 -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL10 -angle=30 -xa=-5 -xb=5 -interp=e -numGhost=3 -order=4 -factor=16

# Short domain L=4 : [-2,2]
ellipseL4Angle30:
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL4 -angle=30 -xa=-2. -xb=2. -interp=e -order=2 -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL4 -angle=30 -xa=-2. -xb=2. -interp=e -order=2 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL4 -angle=30 -xa=-2. -xb=2. -interp=e -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL4 -angle=30 -xa=-2. -xb=2. -interp=e -order=4 -factor=16
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL4 -angle=30 -xa=-2. -xb=2. -interp=e -numGhost=3 -order=4 -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL4 -angle=30 -xa=-2. -xb=2. -interp=e -numGhost=3 -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL4 -angle=30 -xa=-2. -xb=2. -interp=e -numGhost=3 -order=4 -factor=16

# Short domain L=4 : [-1,1]
ellipseL2Angle30:
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL2 -angle=30 -xa=-1 -xb=1 -interp=e -order=2 -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL2 -angle=30 -xa=-1 -xb=1 -interp=e -order=2 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL2 -angle=30 -xa=-1 -xb=1 -interp=e -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL2 -angle=30 -xa=-1 -xb=1 -interp=e -order=4 -factor=16
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL2 -angle=30 -xa=-1 -xb=1 -interp=e -numGhost=3 -order=4 -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL2 -angle=30 -xa=-1 -xb=1 -interp=e -numGhost=3 -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridL2p -angle=30 -xa=-1 -xb=1 -interp=e -numGhost=3 -order=4 -factor=16
 
# Non periodic grid for oblique incidence test 
ellipsenp:
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridn -per=0 -xa=-2. -xb=2. -interp=e -factor=4
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridn -per=0  -xa=-2. -xb=2. -interp=e -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridn -per=0  -xa=-2. -xb=2. -interp=e -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=ellipse -prefix=ellipticalGridn -per=0  -xa=-2. -xb=2. -interp=e -order=4 -factor=16



rod:
	$(ogen) -noplot rodGrid2d -order=2 -xa=-5 -xb=5 -interp=e -factor=4
	$(ogen) -noplot rodGrid2d -order=2 -xa=-5 -xb=5 -interp=e -factor=8
	$(ogen) -noplot rodGrid2d -order=2 -xa=-5 -xb=5 -interp=e -factor=16
	$(ogen) -noplot rodGrid2d -order=4 -xa=-5 -xb=5 -interp=e -factor=8
	$(ogen) -noplot rodGrid2d -order=4 -xa=-5 -xb=5 -interp=e -factor=16

cross:
	$(ogen) -noplot shapedDiskGrid -shape=cross -prefix=crossGrid -xa=-5 -xb=5 -interp=e -order=2 -factor=16
	$(ogen) -noplot shapedDiskGrid -shape=cross -prefix=crossGrid -xa=-5 -xb=5 -interp=i -order=4 -factor=16
	$(ogen) -noplot shapedDiskGrid -shape=cross -prefix=crossGrid -xa=-5 -xb=5 -interp=e -order=4 -factor=16

# IBeam -- shorter domain Nov 26, 2020
IBeam:
	$(ogen) -noplot shapedDiskGrid -shape=IBeam -prefix=IBeamGrid -xa=-1 -xb=1 -interp=e -order=2 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=IBeam -prefix=IBeamGrid -xa=-1 -xb=1 -interp=e -order=2 -factor=16
	$(ogen) -noplot shapedDiskGrid -shape=IBeam -prefix=IBeamGrid -xa=-1 -xb=1 -interp=i -order=4 -factor=8
	$(ogen) -noplot shapedDiskGrid -shape=IBeam -prefix=IBeamGrid -xa=-1 -xb=1 -interp=e -order=4 -factor=16
	$(ogen) -noplot shapedDiskGrid -shape=IBeam -prefix=IBeamGrid -xa=-1 -xb=1 -interp=e -order=4 -numGhost=3 -factor=16

sphere:
	$(ogen) -noplot solidSphereInABox -prefix=dielectricSphereGrid -xa=-3 -xb=3 -ya=-1.5 -yb=1.5 -za=-1.5 -zb=1.5 -outerBC=yzPeriodic -order=2 -interp=e -factor=2
	$(ogen) -noplot solidSphereInABox -prefix=dielectricSphereGrid -xa=-3 -xb=3 -ya=-1.5 -yb=1.5 -za=-1.5 -zb=1.5 -outerBC=yzPeriodic -order=2 -interp=e -factor=4
	$(ogen) -noplot solidSphereInABox -prefix=dielectricSphereGrid -xa=-3 -xb=3 -ya=-1.5 -yb=1.5 -za=-1.5 -zb=1.5 -outerBC=yzPeriodic -order=2 -interp=e -factor=8

pillbox:
	$(ogen) -noplot solidPillBoxGrid.cmd -interp=e -factor=4 

# NOTE: USE -overlapOption=maximize to avoid trouble with matching at interface: 
# ellipsoid543: a=.5, b=.4, c=.3
ellipsoid543Grid:
	$(ogen) -noplot solidEllipsoidInABoxGrid -prefix=solidEllipsoid543Grid -overlapOption=maximize -order=2 -per=npp -a=.5 -b=.4 -c=.3 -factor=4
	$(ogen) -noplot solidEllipsoidInABoxGrid -prefix=solidEllipsoid543Grid -overlapOption=maximize -order=2 -per=npp -a=.5 -b=.4 -c=.3 -factor=8

# to avoid negative volumes at order=4, start at factor=8 *fix me*
rpiGrids:
	$(ogen) -noplot rpiGrid -direction=vertical -interp=e -order=2 -factor=4
	$(ogen) -noplot rpiGrid -direction=vertical -interp=e -numGhost=3 -order=4 -factor=8


# Probe runs
#   Single dielectric block, eps=1, 2, 3, 4
probeBlock:
	$(cgmx) -noplot dielectricBodies -g=blockGride8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps4G8 -eps1=4 -go=go > blockEps4G8.out
	$(cgmx) -noplot dielectricBodies -g=blockGride4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps4G4 -eps1=4 -go=go > blockEps4G4.out
	$(cgmx) -noplot dielectricBodies -g=blockGride2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps4G2 -eps1=4 -go=go > blockEps4G2.out
	$(cgmx) -noplot dielectricBodies -g=blockGride8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps3G8 -eps1=3 -go=go > blockEps3G8.out
	$(cgmx) -noplot dielectricBodies -g=blockGride4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps3G4 -eps1=3 -go=go > blockEps4G3.out
	$(cgmx) -noplot dielectricBodies -g=blockGride2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps3G2 -eps1=3 -go=go > blockEps4G3.out
	$(cgmx) -noplot dielectricBodies -g=blockGride8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps2G8 -eps1=2 -go=go > blockEps2G8.out
	$(cgmx) -noplot dielectricBodies -g=blockGride4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps2G4 -eps1=2 -go=go > blockEps2G4.out
	$(cgmx) -noplot dielectricBodies -g=blockGride2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps2G2 -eps1=2 -go=go > blockEps2G2.out
	$(cgmx) -noplot dielectricBodies -g=blockGride8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps1G8 -eps1=1 -go=go > blockEps1G8.out
	$(cgmx) -noplot dielectricBodies -g=blockGride4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps1G4 -eps1=1 -go=go > blockEps1G4.out
	$(cgmx) -noplot dielectricBodies -g=blockGride2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=-2 -rbc=rbcNonLocal -background=OuterSquare -xLeftProbe=-2 -xRightProbe=2 -dm=none -probeFrequency=1 -probeFileName=blockEps1G2 -eps1=1 -go=go > blockEps1G2.out

# Probe runs in 3d 
#   Single dielectric block, eps=1, 2, 3, 4
#   Polarization:
#     P1 or PEy :  Ex-Ey-Hz  -->  x-direction: (0,Ey,Hz)
#     P2 or PEz :  Hx-Hy-Ez  -->  x-direction: (0,Hy,Ez) 
# Notes:
#   using abcEM2 + probeLocation=-2 : looks 2nd order
probeBlock3dEps1:
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps1G2 -eps1=1 -go=go > block3dP1Eps1G2.out
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps1G4 -eps1=1 -go=go > block3dP1Eps1G4.out
#	$(cgmx) -noplot dielectricBodies -g=blockGrid3de8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps1G8 -eps1=1 -go=go > block3dP1Eps1G8.out

probeBlock3dEps2:
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps2G2 -eps1=2 -go=go > block3dP1Eps2G2.out
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps2G4 -eps1=2 -go=go > block3dP1Eps2G4.out
#	$(cgmx) -noplot dielectricBodies -g=blockGrid3de8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps2G8 -eps1=2 -go=go > block3dP1Eps2G8.out


probeBlock3dEps4:
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps4G2 -eps1=4 -go=go > block3dP1Eps4G2.out
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps4G4 -eps1=4 -go=go > block3dP1Eps4G4.out
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps4G8 -eps1=4 -go=go > block3dP1Eps4G8.out


# PML VERSIONS
probeBlock3dEps1PML:
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcPML -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps1G2 -eps1=1 -go=go > block3dP1Eps1G2.out
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcPML -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps1G4 -eps1=1 -go=go > block3dP1Eps1G4.out
#	$(cgmx) -noplot dielectricBodies -g=blockGrid3de8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcPML -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps1G8 -eps1=1 -go=go > block3dP1Eps1G8.out

probeBlock3dEps2PML:
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcPML -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps2G2 -eps1=2 -go=go > block3dP1Eps2G2.out
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcPML -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps2G4 -eps1=2 -go=go > block3dP1Eps2G4.out
#	$(cgmx) -noplot dielectricBodies -g=blockGrid3de8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcPML -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps2G8 -eps1=2 -go=go > block3dP1Eps2G8.out


probeBlock3dEps4PML:
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcPML -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps4G2 -eps1=4 -go=go > block3dP1Eps4G2.out
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcPML -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps4G4 -eps1=4 -go=go > block3dP1Eps4G4.out
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -xb=1 -rbc=abcPML -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP1Eps4G8 -eps1=4 -go=go > block3dP1Eps4G8.out



# **** POLARIZATION 2 : Ez 
probeBlock3dP2Eps1:
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -az=1 -xb=-1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP2Eps1G2 -eps1=1 -go=go > block3dP2Eps1G2.out
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -az=1 -xb=1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP2Eps1G4 -eps1=1 -go=go > block3dP2Eps1G4.out
#	$(cgmx) -noplot dielectricBodies -g=blockGrid3de8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -az=1 -xb=-1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP2Eps1G8 -eps1=1 -go=go > block3dP2Eps1G8.out

probeBlock3dP2Eps2:
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de2.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -az=1 -xb=-1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP2Eps2G2 -eps1=2 -go=go > block3dP2Eps2G2.out
	$(cgmx) -noplot dielectricBodies -g=blockGrid3de4.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -az=1 -xb=-1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP2Eps2G4 -eps1=2 -go=go > block3dP2Eps2G4.out
#	$(cgmx) -noplot dielectricBodies -g=blockGrid3de8.order4.ng3.hdf -tf=30 -tp=1 -useSosupDissipation=1 -ic=pw -kx=1 -az=1 -xb=-1 -rbc=abcEM2 -background=OuterSquare -xLeftProbe=-1 -xRightProbe=1 -dm=none -probeFrequency=1 -probeFileName=block3dP2Eps2G8 -eps1=2 -go=go > block3dP2Eps2G8.out

