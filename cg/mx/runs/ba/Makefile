# Makefile for BA tests


# ogen = /home/henshw/Overture.g/bin/ogen
ogen = $(Overture)/bin/ogen
cgmx = ${CGBUILDPREFIX}/mx/bin/cgmx

bigSmall: bigSquare bigSmallWidth0p05 bigSmallWidth0p1 bigSmallWidth0p2

# BA GDM  Material WIDTH=.2 
bigSmallWidth0p2:
	$(cgmx) -noplot baScat -g=cSquareSmalle16.order4.ng3 -method=bamx -dm=gdm -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=10 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -compareToShowFile=cSquareBaMatGDMGenNdDRef16.show -show=baMatGDMGenNDSGWidth0p2.show -go=go > baMatGDMGenNDSGWidth0p2.out


# BA Material WIDTH=.1
bigSmallWidth0p1:
	$(cgmx) -noplot baScat -g=cSquareSmalle16.order4.ng3 -method=bamx -dm=gdm -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=10 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.1 -compareToShowFile=cSquareBaMatGDMGenNdDRef16.show -show=baMatGDMGenNDSGWidth0p1.show -go=go > baMatGDMGenNDSGWidth0p1.out

# BA Material WIDTH=.05
bigSmallWidth0p05:
	$(cgmx) -noplot baScat -g=cSquareSmalle16.order4.ng3 -method=bamx -dm=gdm -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=10 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.05 -compareToShowFile=cSquareBaMatGDMGenNdDRef16.show -show=baMatGDMGenNDSGWidth0p05.show -go=go > baMatGDMGenNDSGWidth0p05.out

# Big domain for reference 
bigSquare:
	$(cgmx) -noplot baScat -g=cSquareBige16.order4.ng3 -method=bamx  -dm=gdm -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=10 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.1 -show=cSquareBaMatGDMGenNdDRef16.show -go=go >! cSquareBaMatGenRef16.out


# 2D baMatGDMgenND for self convergence 
gsgdm2d:
	$(cgmx) -noplot baScat -g=cSquareSmalle4.order4.ng3 -method=bamx -dm=none -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatGDMgenNDG4.show -go=go > baMatGDMgenNDG4.out
	$(cgmx) -noplot baScat -g=cSquareSmalle8.order4.ng3 -method=bamx -dm=none -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatGDMgenNDG8.show -go=go > baMatGDMgenNDG8.out
	$(cgmx) -noplot baScat -g=cSquareSmalle16.order4.ng3 -method=bamx -dm=none -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatGDMgenNDG16.show -go=go > baMatGDMgenNDG16.out
	$(cgmx) -noplot baScat -g=cSquareSmalle32.order4.ng3 -method=bamx -dm=none -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatGDMgenNDG32.show -go=go > baMatGDMgenNDG32.out
	$(cgmx) -noplot baScat -g=cSquareSmalle64.order4.ng3 -method=bamx -dm=none -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatGDMgenNDG64.show -go=go > baMatGDMgenNDG64.out


# 2D Isotropic for self convergence 
gs2dIso:
	$(cgmx) -noplot baScat -g=cSquareSmalle4.order4.ng3 -method=bamx -dm=none -matFile=baIsoEps1.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatIsoEps1G4.show -go=go > baMatIsoEps1G4.out
	$(cgmx) -noplot baScat -g=cSquareSmalle8.order4.ng3 -method=bamx -dm=none -matFile=baIsoEps1.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatIsoEps1G8.show -go=go > baMatIsoEps1G8.out
	$(cgmx) -noplot baScat -g=cSquareSmalle16.order4.ng3 -method=bamx -dm=none -matFile=baIsoEps1.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatIsoEps1G16.show -go=go > baMatIsoEps1G16.out
	$(cgmx) -noplot baScat -g=cSquareSmalle32.order4.ng3 -method=bamx -dm=none -matFile=baIsoEps1.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatIsoEps1G32.show -go=go > baMatIsoEps1G32.out
#	$(cgmx) -noplot baScat -g=cSquareSmalle64.order4.ng3 -method=bamx -dm=none -matFile=baIsoEps1.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatIsoEps1G64.show -go=go > baMatIsoEps1G64.out

# 3D Isotropic for self convergence 
gs3dIso:
	$(cgmx) -noplot baScat -g=cboxSmalle4.order4.ng3 -method=bamx -dm=none -matFile=baIsoEps1.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatIsoEps1BoxG4.show -go=go > baMatIsoEps1BoxG4.out
	$(cgmx) -noplot baScat -g=cboxSmalle8.order4.ng3 -method=bamx -dm=none -matFile=baIsoEps1.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatIsoEps1BoxG8.show -go=go > baMatIsoEps1BoxG8.out
	$(cgmx) -noplot baScat -g=cboxSmalle16.order4.ng3 -method=bamx -dm=none -matFile=baIsoEps1.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatIsoEps1BoxG16.show -go=go > baMatIsoEps1BoxG16.out
	$(cgmx) -noplot baScat -g=cboxSmalle32.order4.ng3 -method=bamx -dm=none -matFile=baIsoEps1.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatIsoEps1BoxG32.show -go=go > baMatIsoEps1BoxG32.out
#	$(cgmx) -noplot baScat -g=cboxSmalle64.order4.ng3 -method=bamx -dm=none -matFile=baIsoEps1.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatIsoEps1BoxG64.show -go=go > baMatIsoEps1BoxG64.out


# 3D baMatGDMgenND for self convergence 
gsgdm3d:
	$(cgmx) -noplot baScat -g=cboxSmalle4.order4.ng3 -method=bamx -dm=none -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatGDMgenNDBoxG4.show -go=go > baMatGDMgenNDBoxG4.out
	$(cgmx) -noplot baScat -g=cboxSmalle8.order4.ng3 -method=bamx -dm=none -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatGDMgenNDBoxG8.show -go=go > baMatGDMgenNDBoxG8.out
	$(cgmx) -noplot baScat -g=cboxSmalle16.order4.ng3 -method=bamx -dm=none -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatGDMgenNDBoxG16.show -go=go > baMatGDMgenNDBoxG16.out
	$(cgmx) -noplot baScat -g=cboxSmalle32.order4.ng3 -method=bamx -dm=none -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatGDMgenNDBoxG32.show -go=go > baMatGDMgenNDBoxG32.out
#	$(cgmx) -noplot baScat -g=cboxSmalle64.order4.ng3 -method=bamx -dm=none -matFile=baMatGDMgenND.txt -ts=rk4 -eps1=1 -eps2=1 -ic=gs -x0=.0 -y0=.0 -beta=200 -solveForAllFields=1 -tp=.1 -tf=2 -diss=1 -rbc=absorbing  -useSuperGrid=1 -superGridWidth=.2 -show=baMatGDMgenNDBoxG64.show -go=go > baMatGDMgenNDBoxG64.out



# Make a 2D rectangular region 
rectangle:
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -factor=2
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -factor=4
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -factor=8
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -factor=16
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -factor=32
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -factor=2
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -factor=4
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -factor=8
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -factor=16
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -factor=32

rectanglep:
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=p -factor=2
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=p -factor=4
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=p -factor=8
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=p -factor=16
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=p -factor=32
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=p -factor=2
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=p -factor=4
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=p -factor=8
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=p -factor=16
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=p -factor=32

rectanglenp:
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=2
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=4
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=8
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=16
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=32
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=64
	 $(ogen) -noplot rectangleArg -prefix=rect -order=2 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=128
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=2
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=4
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=8
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=16
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=32
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=64
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=128 

bigRect:
	 $(ogen) -noplot rectangleArg -prefix=bigRect -order=2 -xa=-2 -xb=2 -ya=-1  -yb=1  -factor=2
	 $(ogen) -noplot rectangleArg -prefix=bigRect -order=2 -xa=-2 -xb=2 -ya=-1  -yb=1  -factor=4
	 $(ogen) -noplot rectangleArg -prefix=bigRect -order=2 -xa=-2 -xb=2 -ya=-1  -yb=1  -factor=8
	 $(ogen) -noplot rectangleArg -prefix=bigRect -order=2 -xa=-2 -xb=2 -ya=-1  -yb=1  -factor=16
	 $(ogen) -noplot rectangleArg -prefix=bigRect -order=2 -xa=-2 -xb=2 -ya=-1  -yb=1  -factor=32
	 $(ogen) -noplot rectangleArg -prefix=bigRect -numGhost=3 -order=4 -xa=-2 -xb=2 -ya=-1  -yb=1  -factor=2
	 $(ogen) -noplot rectangleArg -prefix=bigRect -numGhost=3 -order=4 -xa=-2 -xb=2 -ya=-1  -yb=1  -factor=4
	 $(ogen) -noplot rectangleArg -prefix=bigRect -numGhost=3 -order=4 -xa=-2 -xb=2 -ya=-1  -yb=1  -factor=8
	 $(ogen) -noplot rectangleArg -prefix=bigRect -numGhost=3 -order=4 -xa=-2 -xb=2 -ya=-1  -yb=1  -factor=16
	 $(ogen) -noplot rectangleArg -prefix=bigRect -numGhost=3 -order=4 -xa=-2 -xb=2 -ya=-1  -yb=1  -factor=32

# adjust grid so interface points at x=-.5 and x=.5 are on cell centers
blockGrid:
	 $(ogen) -noplot rectangleArg -prefix=blockGrid -order=2 -xa=-3 -xb=3 -ya=-.1  -yb=.1 -periodic=np -adjustCC=1 -factor=2
	 $(ogen) -noplot rectangleArg -prefix=blockGrid -order=2 -xa=-3 -xb=3 -ya=-.1  -yb=.1 -periodic=np -adjustCC=1 -factor=4
	 $(ogen) -noplot rectangleArg -prefix=blockGrid -order=2 -xa=-3 -xb=3 -ya=-.1  -yb=.1 -periodic=np -adjustCC=1 -factor=8
	 $(ogen) -noplot rectangleArg -prefix=blockGrid -order=2 -xa=-3 -xb=3 -ya=-.1  -yb=.1 -periodic=np -adjustCC=1 -factor=16
	 $(ogen) -noplot rectangleArg -prefix=blockGrid -order=2 -xa=-3 -xb=3 -ya=-.1  -yb=.1 -periodic=np -adjustCC=1 -factor=32
	 $(ogen) -noplot rectangleArg -prefix=blockGrid -numGhost=3 -order=4 -xa=-3 -xb=3 -ya=-.1  -yb=.1 -periodic=np -adjustCC=1 -factor=2
	 $(ogen) -noplot rectangleArg -prefix=blockGrid -numGhost=3 -order=4 -xa=-3 -xb=3 -ya=-.1  -yb=.1 -periodic=np -adjustCC=1 -factor=4
	 $(ogen) -noplot rectangleArg -prefix=blockGrid -numGhost=3 -order=4 -xa=-3 -xb=3 -ya=-.1  -yb=.1 -periodic=np -adjustCC=1 -factor=8
	 $(ogen) -noplot rectangleArg -prefix=blockGrid -numGhost=3 -order=4 -xa=-3 -xb=3 -ya=-.1  -yb=.1 -periodic=np -adjustCC=1 -factor=16
	 $(ogen) -noplot rectangleArg -prefix=blockGrid -numGhost=3 -order=4 -xa=-3 -xb=3 -ya=-.1  -yb=.1 -periodic=np -adjustCC=1 -factor=32


# (centered) square for the domain [-1,1]x[-1,1]
cSquare:
	 $(ogen) -noplot rectangleArg -prefix=cSquaree -order=2 -xa=-1 -xb=1 -ya=-1 -yb=1 -factor=2
	 $(ogen) -noplot rectangleArg -prefix=cSquaree -order=2 -xa=-1 -xb=1 -ya=-1 -yb=1 -factor=4
	 $(ogen) -noplot rectangleArg -prefix=cSquaree -order=2 -xa=-1 -xb=1 -ya=-1 -yb=1 -factor=8
	 $(ogen) -noplot rectangleArg -prefix=cSquaree -order=2 -xa=-1 -xb=1 -ya=-1 -yb=1 -factor=16
	 $(ogen) -noplot rectangleArg -prefix=cSquaree -order=2 -xa=-1 -xb=1 -ya=-1 -yb=1 -factor=32
	 $(ogen) -noplot rectangleArg -prefix=cSquaree -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-1 -yb=1 -factor=2
	 $(ogen) -noplot rectangleArg -prefix=cSquaree -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-1 -yb=1 -factor=4
	 $(ogen) -noplot rectangleArg -prefix=cSquaree -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-1 -yb=1 -factor=8
	 $(ogen) -noplot rectangleArg -prefix=cSquaree -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-1 -yb=1 -factor=16
	 $(ogen) -noplot rectangleArg -prefix=cSquaree -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-1 -yb=1 -factor=32

# small square for the domain [-.5,.5]^2
cSquareSmall:
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -order=2 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=2
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -order=2 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=4
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -order=2 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=8
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -order=2 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=16
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -order=2 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=32
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -order=2 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=64
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -numGhost=3 -order=4 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=2
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -numGhost=3 -order=4 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=4
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -numGhost=3 -order=4 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=8
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -numGhost=3 -order=4 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=16
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -numGhost=3 -order=4 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=32
	 $(ogen) -noplot rectangleArg -prefix=cSquareSmalle -numGhost=3 -order=4 -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -factor=64
# big square for the domain [4,4]^2
cSquareBig:
	 $(ogen) -noplot rectangleArg -prefix=cSquareBige -order=2 -xa=-4 -xb=4 -ya=-4 -yb=4 -factor=2
	 $(ogen) -noplot rectangleArg -prefix=cSquareBige -order=2 -xa=-4 -xb=4 -ya=-4 -yb=4 -factor=4
	 $(ogen) -noplot rectangleArg -prefix=cSquareBige -order=2 -xa=-4 -xb=4 -ya=-4 -yb=4 -factor=8
	 $(ogen) -noplot rectangleArg -prefix=cSquareBige -order=2 -xa=-4 -xb=4 -ya=-4 -yb=4 -factor=16
	 $(ogen) -noplot rectangleArg -prefix=cSquareBige -order=2 -xa=-4 -xb=4 -ya=-4 -yb=4 -factor=32
	 $(ogen) -noplot rectangleArg -prefix=cSquareBige -numGhost=3 -order=4 -xa=-4 -xb=4 -ya=-4 -yb=4 -factor=2
	 $(ogen) -noplot rectangleArg -prefix=cSquareBige -numGhost=3 -order=4 -xa=-4 -xb=4 -ya=-4 -yb=4 -factor=4
	 $(ogen) -noplot rectangleArg -prefix=cSquareBige -numGhost=3 -order=4 -xa=-4 -xb=4 -ya=-4 -yb=4 -factor=8
	 $(ogen) -noplot rectangleArg -prefix=cSquareBige -numGhost=3 -order=4 -xa=-4 -xb=4 -ya=-4 -yb=4 -factor=16
	 $(ogen) -noplot rectangleArg -prefix=cSquareBige -numGhost=3 -order=4 -xa=-4 -xb=4 -ya=-4 -yb=4 -factor=32


# Make a 3D rectangular region for an interface problem
box3d:
	$(ogen) -noplot boxArg -prefix=box3de -xa=-1 -xb=1 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -order=2 -factor=1 
	$(ogen) -noplot boxArg -prefix=box3de -xa=-1 -xb=1 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -order=2 -factor=2 
	$(ogen) -noplot boxArg -prefix=box3de -xa=-1 -xb=1 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -order=2 -factor=4 
	$(ogen) -noplot boxArg -prefix=box3de -xa=-1 -xb=1 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -order=2 -factor=8 
	$(ogen) -noplot boxArg -prefix=box3de -xa=-1 -xb=1 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -order=2 -factor=16
	$(ogen) -noplot boxArg -prefix=box3de -xa=-1 -xb=1 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -numGhost=3 -order=4 -factor=1 
	$(ogen) -noplot boxArg -prefix=box3de -xa=-1 -xb=1 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -numGhost=3 -order=4 -factor=2 
	$(ogen) -noplot boxArg -prefix=box3de -xa=-1 -xb=1 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -numGhost=3 -order=4 -factor=4 
	$(ogen) -noplot boxArg -prefix=box3de -xa=-1 -xb=1 -ya=-.5 -yb=.5 -za=-.5 -zb=.5  -numGhost=3-order=4 -factor=8 
	$(ogen) -noplot boxArg -prefix=box3de -xa=-1 -xb=1 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -numGhost=3 -order=4 -factor=16


# for BA scattering from a sphere 
cbox:
	$(ogen) -noplot boxArg -prefix=cboxe -xa=-2 -xb=2 -ya=-2 -yb=2 -za=-2 -zb=2 -order=2 -factor=1 
	$(ogen) -noplot boxArg -prefix=cboxe -xa=-2 -xb=2 -ya=-2 -yb=2 -za=-2 -zb=2 -order=2 -factor=2 
	$(ogen) -noplot boxArg -prefix=cboxe -xa=-2 -xb=2 -ya=-2 -yb=2 -za=-2 -zb=2 -order=2 -factor=4 
	$(ogen) -noplot boxArg -prefix=cboxe -xa=-2 -xb=2 -ya=-2 -yb=2 -za=-2 -zb=2 -order=2 -factor=8 
	$(ogen) -noplot boxArg -prefix=cboxe -xa=-2 -xb=2 -ya=-2 -yb=2 -za=-2 -zb=2 -order=2 -factor=16
	$(ogen) -noplot boxArg -prefix=cboxe -xa=-2 -xb=2 -ya=-2 -yb=2 -za=-2 -zb=2 -numGhost=3 -order=4 -factor=1 
	$(ogen) -noplot boxArg -prefix=cboxe -xa=-2 -xb=2 -ya=-2 -yb=2 -za=-2 -zb=2 -numGhost=3 -order=4 -factor=2 
	$(ogen) -noplot boxArg -prefix=cboxe -xa=-2 -xb=2 -ya=-2 -yb=2 -za=-2 -zb=2 -numGhost=3 -order=4 -factor=4 
	$(ogen) -noplot boxArg -prefix=cboxe -xa=-2 -xb=2 -ya=-2 -yb=2 -za=-2 -zb=2 -numGhost=3 -order=4 -factor=8 
	$(ogen) -noplot boxArg -prefix=cboxe -xa=-2 -xb=2 -ya=-2 -yb=2 -za=-2 -zb=2 -numGhost=3 -order=4 -factor=16


# for BA scattering from a sphere 
cboxSmall:
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -order=2 -factor=1 
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -order=2 -factor=2 
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -order=2 -factor=4 
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -order=2 -factor=8 
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -order=2 -factor=16
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -order=2 -factor=32
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -numGhost=3 -order=4 -factor=1 
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -numGhost=3 -order=4 -factor=2 
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -numGhost=3 -order=4 -factor=4 
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -numGhost=3 -order=4 -factor=8 
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -numGhost=3 -order=4 -factor=16
	$(ogen) -noplot boxArg -prefix=cboxSmalle -xa=-.5 -xb=.5 -ya=-.5 -yb=.5 -za=-.5 -zb=.5 -numGhost=3 -order=4 -factor=32

# rectangle with refinements
rectRefine:
	$(ogen) -noplot rectRefine -order=2 -factor=1
	$(ogen) -noplot rectRefine -order=2 -factor=2
	$(ogen) -noplot rectRefine -order=2 -factor=4
	$(ogen) -noplot rectRefine -order=2 -factor=8
	$(ogen) -noplot rectRefine -order=4 -factor=1
	$(ogen) -noplot rectRefine -order=4 -factor=2
	$(ogen) -noplot rectRefine -order=4 -factor=4
	$(ogen) -noplot rectRefine -order=4 -factor=8
