#
#   Makefile for the soliton examples
#


ogen = $(Overture)/bin/ogen
plotStuff = $(Overture)/bin/plotStuff
cgmx = $(CGBUILDPREFIX)/mx/bin/cgmx

# Make a 2D rectangular grids
grids:
	$(ogen) -noplot rectangleArg -prefix=thinRectangle -order=2 -xa=0 -xb=400 -ya=0.0 -yb=.5 -periodic=np -factor=1
	$(ogen) -noplot rectangleArg -prefix=thinRectangle -order=2 -xa=0 -xb=400 -ya=0.0 -yb=.5 -periodic=np -factor=2
	$(ogen) -noplot rectangleArg -prefix=thinRectangle -order=2 -xa=0 -xb=400 -ya=0.0 -yb=.5 -periodic=np -factor=4
	$(ogen) -noplot rectangleArg -prefix=thinRectangle -order=4 -xa=0 -xb=400 -ya=0.0 -yb=.5 -periodic=np -numGhost=3 -factor=1
	$(ogen) -noplot rectangleArg -prefix=thinRectangle -order=4 -xa=0 -xb=400 -ya=0.0 -yb=.5 -periodic=np -numGhost=3 -factor=2
	$(ogen) -noplot rectangleArg -prefix=thinRectangle -order=4 -xa=0 -xb=400 -ya=0.0 -yb=.5 -periodic=np -numGhost=3 -factor=4

# Order=2 runs: 
runs2:
	$(cgmx) -noplot soliton -g=thinRectangle1np.order2 -tf=400 -tp=10 -bc=d -diss=0 -method=fd -dm=gdm -npv=1 -nm=mla -matFile1=mlaSolitonMat3.txt -ic=pulse -show=solitonO2G1.show -go=go > solitonO2G1.out 


# Order=4 runs: 
runs4:
	$(cgmx) -noplot soliton -g=thinRectangle1np.order4.ng3 -tf=400 -tp=10 -bc=d -diss=0 -method=fd -dm=gdm -npv=1 -nm=mla -matFile1=mlaSolitonMat3.txt -ic=pulse -show=solitonO4G1.show -go=go > solitonO4G1.out

# plot results from runs 
plots:
	$(plotStuff) plotSoliton.cmd -show=solitonO4 -name=solitonO2 -solution=21
	$(plotStuff) plotSoliton.cmd -show=solitonO4 -name=solitonO4 -solution=21
