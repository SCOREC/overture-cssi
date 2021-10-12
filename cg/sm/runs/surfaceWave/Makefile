# ogen = /home/henshw/Overture.g/bin/ogen
ogen = $(Overture)/bin/ogen

cgsm= $(CGBUILDPREFIX)/sm/bin/cgsm

grids: case1Grid case1GridCurvilinear stripTractionTractionGrid

# strip - case 1 grid (traction-displacement BC)
case1Grid:
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=2
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=4
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=8
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=16
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=32
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=64
	$(ogen) -noplot rectangleArg -prefix=case1Grid40x40 -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=4 -numGhost=3 -Nx=41 -Ny=41
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=4 -numGhost=3 -factor=2
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=4 -numGhost=3 -factor=4
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=4 -numGhost=3 -factor=8
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=4 -numGhost=3 -factor=16
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=4 -numGhost=3 -factor=32	
	$(ogen) -noplot rectangleArg -prefix=case1Grid -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=4 -numGhost=3 -factor=64


# case 1 grid
case1GridCurvilinear:
	$(ogen) -noplot nonRectangle -prefix=case1GridCurvilinear -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=2
	$(ogen) -noplot nonRectangle -prefix=case1GridCurvilinear -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=4
	$(ogen) -noplot nonRectangle -prefix=case1GridCurvilinear -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=8
	$(ogen) -noplot nonRectangle -prefix=case1GridCurvilinear -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=16
	$(ogen) -noplot nonRectangle -prefix=case1GridCurvilinear -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=32	
	$(ogen) -noplot nonRectangle -prefix=case1GridCurvilinear -xa=0 -xb=1 -ya=0 -yb=1.2710078982637464465 -periodic=np -order=2 -factor=64	



# strip case 2 grid (traction-traction BC)
stripTractionTractionGrid:
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=2 -factor=2
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=2 -factor=4
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=2 -factor=8
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=2 -factor=16
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=2 -factor=32
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=2 -factor=64
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=4 -numGhost=3 -factor=2
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=4 -numGhost=3 -factor=4
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=4 -numGhost=3 -factor=8
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=4 -numGhost=3 -factor=16
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=4 -numGhost=3 -factor=32	
	$(ogen) -noplot rectangleArg -prefix=stripTractionTractionGrid -xa=0 -xb=1 -ya=0 -yb=1.0 -periodic=np -order=4 -numGhost=3 -factor=64
