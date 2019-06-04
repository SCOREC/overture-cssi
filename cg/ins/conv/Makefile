common = ../../common

# The convergence rate script is called conv.p
all:  
	ln -sf $(common)/src/conv.p .

cgins:
	@cd ..;      $(MAKE) -$(MAKEFLAGS) cgins

ogen = $(Overture)/bin/ogen


taylorGreen:
	conv.p exact.conv -option=tg -ts=pc   -nu=.05 -kp=2 -tf=.2 -tp=.2 -psolver=yale -cg=squarep -order=4 -orderInTime=4 -ng=5
	conv.p exact.conv -option=tg -ts=imex -nu=.05 -kp=2 -tf=.2 -tp=.2 -psolver=yale -cg=squarep -order=2 -orderInTime=2 -ng=5
	conv.p exact.conv -option=tg -ts=imex -nu=.05 -kp=2 -tf=.2 -tp=.2 -psolver=yale -cg=squarep -order=4 -orderInTime=4 -ng=5


circularCouette:
	conv.p exact.conv -option=rotating -ts=imex -amp1=1 -nFreq1=1 -nu=.05 -tf=.3 -tp=.1 -psolver=yale -cg=annulus -order=2 -orderInTime=2 -omegaInner=1 -omegaOuter=2  -ng=5
	conv.p exact.conv -option=rotating -ts=imex -amp1=1 -nFreq1=1 -nu=.05 -tf=.3 -tp=.1 -psolver=yale -cg=annulus -order=4 -orderInTime=2 -omegaInner=1 -omegaOuter=2  -ng=5
	conv.p exact.conv -option=rotating -ts=imex -amp1=1 -nFreq1=1 -nu=.05 -tf=.3 -tp=.1 -psolver=yale -cg=annulus -order=4 -orderInTime=4 -omegaInner=1 -omegaOuter=2  -ng=5

circularCouetteDouble:
	conv.p exact.conv -option=rotating -ts=imex -amp1=1 -nFreq1=1 -nu=.05 -tf=.3 -tp=.1 -psolver=yale -cg=splitAnnulus -order=2 -orderInTime=2 -omegaInner=1 -omegaOuter=2  -ng=4
	conv.p exact.conv -option=rotating -ts=imex -amp1=1 -nFreq1=1 -nu=.05 -tf=.3 -tp=.1 -psolver=yale -cg=splitAnnulus -order=4 -orderInTime=2 -omegaInner=1 -omegaOuter=2  -ng=4
	conv.p exact.conv -option=rotating -ts=imex -amp1=1 -nFreq1=1 -nu=.05 -tf=.3 -tp=.1 -psolver=yale -cg=splitAnnulus -order=4 -orderInTime=4 -omegaInner=1 -omegaOuter=2  -ng=4


diskCouette:
	conv.p exact.conv -option=rotating -ts=imex -amp1=.5 -nFreq1=1 -nu=.05 -tf=.3 -tp=.1 -psolver=yale -cg=sic -rInner=0 -order=2 -orderInTime=2 -ng=4
	conv.p exact.conv -option=rotating -ts=imex -amp1=.5 -nFreq1=1 -nu=.05 -tf=.3 -tp=.1 -psolver=yale -cg=sic -rInner=0 -order=4 -orderInTime=2 -ng=4
	conv.p exact.conv -option=rotating -ts=imex -amp1=.5 -nFreq1=1 -nu=.05 -tf=.3 -tp=.1 -psolver=yale -cg=sic -rInner=0 -order=4 -orderInTime=4 -ng=4

taylorCouette3d:
	conv.p exact.conv -option=rotating -ts=imex -order=2 -orderInTime=2 -nu=.05 -tf=.1 -tp=.1 -psolver=best -cg=cylinderp -amp1=1. -nFreq1=1 -ng=3
	conv.p exact.conv -option=rotating -ts=imex -order=4 -orderInTime=2 -nu=.05 -tf=.1 -tp=.1 -psolver=best -cg=cylinderp -amp1=1. -nFreq1=1 -ng=3
	conv.p exact.conv -option=rotating -ts=imex -order=4 -orderInTime=4 -nu=.05 -tf=.1 -tp=.1 -psolver=best -cg=cylinderp -amp1=1. -nFreq1=1 -ng=3

pipeFlow:
	conv.p exact.conv -option=pipe -ts=imex -order=2 -orderInTime=2 -nu=.1 -tf=.1 -tp=.1 -psolver=best -solver=best -cg=shortPipep -pGrad=.1 -length=.5 -amp1=1 -nFreq1=1 -axialAxis=0 -omega=1. -amp2=1 -nFreq2=1 -ng=4
	conv.p exact.conv -option=pipe -ts=imex -order=4 -orderInTime=2 -nu=.1 -tf=.1 -tp=.1 -psolver=best -solver=best -cg=shortPipep -pGrad=.1 -length=.5 -amp1=1 -nFreq1=1 -axialAxis=0 -omega=1. -amp2=1 -nFreq2=1 -ng=4
	conv.p exact.conv -option=pipe -ts=imex -order=4 -orderInTime=4 -nu=.1 -tf=.1 -tp=.1 -psolver=best -solver=best -cg=shortPipep -pGrad=.1 -length=.5 -amp1=1 -nFreq1=1 -axialAxis=0 -omega=1. -amp2=1 -nFreq2=1 -ng=4

grids: channelGrids

channelGrids:
	$(ogen) -noplot channelArg -order=2 -length=2. -per=0 -factor=1
	$(ogen) -noplot channelArg -order=2 -length=2. -per=0 -factor=2
	$(ogen) -noplot channelArg -order=2 -length=2. -per=0 -factor=4
	$(ogen) -noplot channelArg -order=2 -length=2. -per=0 -factor=8
	$(ogen) -noplot channelArg -order=2 -length=2. -per=0 -factor=16
	$(ogen) -noplot channelArg -order=4 -length=2. -per=0 -factor=1
	$(ogen) -noplot channelArg -order=4 -length=2. -per=0 -factor=2
	$(ogen) -noplot channelArg -order=4 -length=2. -per=0 -factor=4
	$(ogen) -noplot channelArg -order=4 -length=2. -per=0 -factor=8
	$(ogen) -noplot channelArg -order=4 -length=2. -per=0 -factor=16
	$(ogen) -noplot channelArg -order=2 -length=2. -per=1 -factor=1
	$(ogen) -noplot channelArg -prefix=channelGridp -order=2 -length=2. -per=1 -factor=2
	$(ogen) -noplot channelArg -prefix=channelGridp -order=2 -length=2. -per=1 -factor=4
	$(ogen) -noplot channelArg -prefix=channelGridp -order=2 -length=2. -per=1 -factor=8
	$(ogen) -noplot channelArg -prefix=channelGridp -order=2 -length=2. -per=1 -factor=16
	$(ogen) -noplot channelArg -prefix=channelGridp -order=4 -length=2. -per=1 -factor=1
	$(ogen) -noplot channelArg -prefix=channelGridp -order=4 -length=2. -per=1 -factor=2
	$(ogen) -noplot channelArg -prefix=channelGridp -order=4 -length=2. -per=1 -factor=4
	$(ogen) -noplot channelArg -prefix=channelGridp -order=4 -length=2. -per=1 -factor=8
	$(ogen) -noplot channelArg -prefix=channelGridp -order=4 -length=2. -per=1 -factor=16


annulusGrids:
	$(ogen) -noplot annulusArg -order=2 -factor=2 
	$(ogen) -noplot annulusArg -order=2 -factor=4 
	$(ogen) -noplot annulusArg -order=2 -factor=8 
	$(ogen) -noplot annulusArg -order=2 -factor=16 
	$(ogen) -noplot annulusArg -order=2 -factor=32
	$(ogen) -noplot annulusArg -order=4 -factor=2 
	$(ogen) -noplot annulusArg -order=4 -factor=4 
	$(ogen) -noplot annulusArg -order=4 -factor=8 
	$(ogen) -noplot annulusArg -order=4 -factor=16 
	$(ogen) -noplot annulusArg -order=4 -factor=32 

cylGrids:
	$(ogen) -noplot cylinder -periodic=p -order=2 -factor=1
	$(ogen) -noplot cylinder -periodic=p -order=2 -factor=2
	$(ogen) -noplot cylinder -periodic=p -order=2 -factor=4
	$(ogen) -noplot cylinder -periodic=p -order=2 -factor=8
	$(ogen) -noplot cylinder -periodic=p -order=4 -factor=1
	$(ogen) -noplot cylinder -periodic=p -order=4 -factor=2
	$(ogen) -noplot cylinder -periodic=p -order=4 -factor=4
	$(ogen) -noplot cylinder -periodic=p -order=4 -factor=8

pipeGrids: pipeGridsOrder2 pipeGridsOrder4

pipeGridsOrder2:
	$(ogen) -noplot pipe -rgd=fixed -order=2 -interp=e -factor=1 
	$(ogen) -noplot pipe -rgd=fixed -order=2 -interp=e -factor=2 
	$(ogen) -noplot pipe -rgd=fixed -order=2 -interp=e -factor=4
	$(ogen) -noplot pipe -rgd=fixed -order=2 -interp=e -factor=8
	$(ogen) -noplot pipe -prefix=pipep -rgd=fixed -order=2 -interp=e -per=1 -factor=1 
	$(ogen) -noplot pipe -prefix=pipep -rgd=fixed -order=2 -interp=e -per=1 -factor=2 
	$(ogen) -noplot pipe -prefix=pipep -rgd=fixed -order=2 -interp=e -per=1 -factor=4 
	$(ogen) -noplot pipe -prefix=pipep -rgd=fixed -order=2 -interp=e -per=1 -factor=8 

pipeGridsOrder4: 
	$(ogen) -noplot pipe -rgd=fixed -order=4 -interp=e -factor=1 
	$(ogen) -noplot pipe -rgd=fixed -order=4 -interp=e -factor=2
	$(ogen) -noplot pipe -rgd=fixed -order=4 -interp=e -factor=4
	$(ogen) -noplot pipe -rgd=fixed -order=4 -interp=e -factor=8
	$(ogen) -noplot pipe -prefix=pipep -rgd=fixed -order=4 -interp=e -per=1 -factor=1 
	$(ogen) -noplot pipe -prefix=pipep -rgd=fixed -order=4 -interp=e -per=1 -factor=2 
	$(ogen) -noplot pipe -prefix=pipep -rgd=fixed -order=4 -interp=e -per=1 -factor=4 
	$(ogen) -noplot pipe -prefix=pipep -rgd=fixed -order=4 -interp=e -per=1 -factor=8 

shortPipeGrids:
	$(ogen) -noplot pipe -prefix=shortPipep -sb=.5 -rgd=fixed -order=2 -interp=e -per=1 -factor=1 
	$(ogen) -noplot pipe -prefix=shortPipep -sb=.5 -rgd=fixed -order=2 -interp=e -per=1 -factor=2 
	$(ogen) -noplot pipe -prefix=shortPipep -sb=.5 -rgd=fixed -order=2 -interp=e -per=1 -factor=4 
	$(ogen) -noplot pipe -prefix=shortPipep -sb=.5 -rgd=fixed -order=2 -interp=e -per=1 -factor=8 
	$(ogen) -noplot pipe -prefix=shortPipep -sb=.5 -rgd=fixed -order=4 -interp=e -per=1 -factor=1 
	$(ogen) -noplot pipe -prefix=shortPipep -sb=.5 -rgd=fixed -order=4 -interp=e -per=1 -factor=2 
	$(ogen) -noplot pipe -prefix=shortPipep -sb=.5 -rgd=fixed -order=4 -interp=e -per=1 -factor=4 
	$(ogen) -noplot pipe -prefix=shortPipep -sb=.5 -rgd=fixed -order=4 -interp=e -per=1 -factor=8 


# Square with an extra ghost for singular problems 
squareGrids: 
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=2
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=3
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=5
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=10
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=20
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=40
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=80
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=160
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=4
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=8
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=16
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=32
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=64
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=128
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=256
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=512
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=1024
	$(ogen) -noplot squareArg -addGhost=1 -order=2 -nx=2048
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=10
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=20
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=40
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=80
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=160
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=8
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=16
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=32
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=64
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=128
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=256
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=512
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=1024
	$(ogen) -noplot squareArg -addGhost=1 -order=4 -nx=2048
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=2
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=3
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=5
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=10
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=20
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=40
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=80
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=160
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=320
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=4
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=8
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=16
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=32
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=64
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=128
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=256
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=512
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=1024
	$(ogen) -noplot squareArg -periodic=p -order=2 -nx=2048
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=10
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=20
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=40
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=80
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=160
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=320
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=8
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=16
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=32
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=64
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=128
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=256
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=512
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=1024
	$(ogen) -noplot squareArg -periodic=p -order=4 -nx=2048
