# Makefile for BA tests


# ogen = /home/henshw/Overture.g/bin/ogen
ogen = $(Overture)/bin/ogen


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
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=2
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=4
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=8
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=16
	 $(ogen) -noplot rectangleArg -prefix=rect -numGhost=3 -order=4 -xa=-1 -xb=1 -ya=-.5 -yb=.5 -periodic=np -factor=32


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
