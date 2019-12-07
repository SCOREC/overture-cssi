SHELL := /usr/bin/tcsh

ogen = ${Overture}/bin/ogen
plotStuff = ${Overture}/bin/plotStuff
cgins = ${CGBUILDPREFIX}/ins/bin/cgins


grids:
	$(ogen) -noplot twoCylindersGrid -order=2 -interp=e -factor=1
	$(ogen) -noplot twoCylindersGrid -order=2 -interp=e -factor=2
	$(ogen) -noplot twoCylindersGrid -order=2 -interp=e -factor=4
	$(ogen) -noplot twoCylindersGrid -order=2 -interp=e -factor=8  
	$(ogen) -noplot twoCylindersGrid -order=2 -interp=e -factor=16  
	$(ogen) -noplot twoCylindersGrid -order=2 -interp=e -factor=32  
	$(ogen) -noplot twoCylindersGrid -order=4 -interp=e -factor=1
	$(ogen) -noplot twoCylindersGrid -order=4 -interp=e -factor=2
	$(ogen) -noplot twoCylindersGrid -order=4 -interp=e -factor=4
	$(ogen) -noplot twoCylindersGrid -order=4 -interp=e -factor=8
	$(ogen) -noplot twoCylindersGrid -order=4 -interp=e -factor=16  
	$(ogen) -noplot twoCylindersGrid -order=4 -interp=e -factor=32

# fixed radius for grids 
gridsFixed:
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=2 -interp=e -factor=1 -deltaRadius0=.4
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=2 -interp=e -factor=2
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=2 -interp=e -factor=4
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=2 -interp=e -factor=8  
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=2 -interp=e -factor=16  
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=2 -interp=e -factor=32  
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=4 -interp=e -factor=1 -deltaRadius0=.6
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=4 -interp=e -factor=2
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=4 -interp=e -factor=4
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=4 -interp=e -factor=8
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=4 -interp=e -factor=16  
	$(ogen) -noplot twoCylindersGrid -rgd=fixed -order=4 -interp=e -factor=32

# plot for paper:
plotsIMEX44:
	$(plotStuff) plotTwoCyls.cmd -show=twoCyl2IMEX44nu0001.show -name=twoCyl2IMEX44 -vorMin=-70 -vorMax=50
	$(plotStuff) plotTwoCyls.cmd -show=twoCyl4IMEX44nu0001.show -name=twoCyl4IMEX44 -vorMin=-70 -vorMax=50
	$(plotStuff) plotTwoCyls.cmd -show=twoCyl8IMEX44nu0001.show -name=twoCyl8IMEX44 -vorMin=-70 -vorMax=50
	$(plotStuff) plotTwoCyls.cmd -show=twoCyl16IMEX44nu0001.show -name=twoCyl16IMEX44 -vorMin=-70 -vorMax=50
#	$(plotStuff) plotTwoCyls.cmd -show=twoCyl32IMEX44nu0001.show -name=twoCyl32IMEX44 -vorMin=-70 -vorMax=50

plotsIMEX22:
	$(plotStuff) plotTwoCyls.cmd -show=twoCyl2IMEX22nu0001.show -name=twoCyl2IMEX22 -vorMin=-70 -vorMax=50
	$(plotStuff) plotTwoCyls.cmd -show=twoCyl4IMEX22nu0001.show -name=twoCyl4IMEX22 -vorMin=-70 -vorMax=50
	$(plotStuff) plotTwoCyls.cmd -show=twoCyl8IMEX22nu0001.show -name=twoCyl8IMEX22 -vorMin=-70 -vorMax=50
	$(plotStuff) plotTwoCyls.cmd -show=twoCyl16IMEX22nu0001.show -name=twoCyl16IMEX22 -vorMin=-70 -vorMax=50
	$(plotStuff) plotTwoCyls.cmd -show=twoCyl32IMEX22nu0001.show -name=twoCyl32IMEX22 -vorMin=-70 -vorMax=50


# one cylinder in a channel
oneCylinderGrids:
	$(ogen) -noplot twoCylindersGrid -prefix=oneCyl -numCyl=1 -x1=0 -y1=0 -xa=-3 -xb=6 -ya=-2 -yb=2 -bl=1 -rgd=fixed -order=2 -interp=e -factor=2
	$(ogen) -noplot twoCylindersGrid -prefix=oneCyl -numCyl=1 -x1=0 -y1=0 -xa=-3 -xb=6 -ya=-2 -yb=2 -bl=1 -rgd=fixed -order=2 -interp=e -factor=4
	$(ogen) -noplot twoCylindersGrid -prefix=oneCyl -numCyl=1 -x1=0 -y1=0 -xa=-3 -xb=6 -ya=-2 -yb=2 -bl=1 -rgd=fixed -order=2 -interp=e -factor=8
	$(ogen) -noplot twoCylindersGrid -prefix=oneCyl -numCyl=1 -x1=0 -y1=0 -xa=-3 -xb=6 -ya=-2 -yb=2 -bl=1 -rgd=fixed -order=2 -interp=e -factor=16
	$(ogen) -noplot twoCylindersGrid -prefix=oneCyl -numCyl=1 -x1=0 -y1=0 -xa=-3 -xb=6 -ya=-2 -yb=2 -bl=1 -rgd=fixed -order=4 -interp=e -factor=2
	$(ogen) -noplot twoCylindersGrid -prefix=oneCyl -numCyl=1 -x1=0 -y1=0 -xa=-3 -xb=6 -ya=-2 -yb=2 -bl=1 -rgd=fixed -order=4 -interp=e -factor=4
	$(ogen) -noplot twoCylindersGrid -prefix=oneCyl -numCyl=1 -x1=0 -y1=0 -xa=-3 -xb=6 -ya=-2 -yb=2 -bl=1 -rgd=fixed -order=4 -interp=e -factor=8
	$(ogen) -noplot twoCylindersGrid -prefix=oneCyl -numCyl=1 -x1=0 -y1=0 -xa=-3 -xb=6 -ya=-2 -yb=2 -bl=1 -rgd=fixed -order=4 -interp=e -factor=16

# one cylinder in a square
oneCylGrids:
	$(ogen) -noplot oneCylGrid -order=2 -interp=e -factor=1
	$(ogen) -noplot oneCylGrid -order=2 -interp=e -factor=2
	$(ogen) -noplot oneCylGrid -order=2 -interp=e -factor=4
	$(ogen) -noplot oneCylGrid -order=2 -interp=e -factor=8
	$(ogen) -noplot oneCylGrid -order=2 -interp=e -factor=16
	$(ogen) -noplot oneCylGrid -order=4 -interp=e -factor=2
	$(ogen) -noplot oneCylGrid -order=4 -interp=e -factor=4
	$(ogen) -noplot oneCylGrid -order=4 -interp=e -factor=8
	$(ogen) -noplot oneCylGrid -order=4 -interp=e -factor=16
	$(ogen) -noplot oneCylGrid -prefix=oneCylStretchedGrid -blf=3 -order=4 -interp=e -factor=2
	$(ogen) -noplot oneCylGrid -prefix=oneCylStretchedGrid -blf=3 -order=4 -interp=e -factor=4
	$(ogen) -noplot oneCylGrid -prefix=oneCylStretchedGrid -blf=3 -order=4 -interp=e -factor=8

oneCylGridsFixed:
	$(ogen) -noplot oneCylGrid -rgd=fixed -order=2 -interp=e -factor=1
	$(ogen) -noplot oneCylGrid -rgd=fixed -order=2 -interp=e -factor=2
	$(ogen) -noplot oneCylGrid -rgd=fixed -order=2 -interp=e -factor=4
	$(ogen) -noplot oneCylGrid -rgd=fixed -order=2 -interp=e -factor=8
	$(ogen) -noplot oneCylGrid -rgd=fixed -order=2 -interp=e -factor=16
	$(ogen) -noplot oneCylGrid -rgd=fixed -order=4 -interp=e -factor=2
	$(ogen) -noplot oneCylGrid -rgd=fixed -order=4 -interp=e -factor=4
	$(ogen) -noplot oneCylGrid -rgd=fixed -order=4 -interp=e -factor=8
	$(ogen) -noplot oneCylGrid -rgd=fixed -order=4 -interp=e -factor=16



# stretched grid 
oneCylMoveIMEX44s:
	$(cgins) -noplot oneCylMove -g=oneCylStretchedGride2.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.01   -tolFactor=1e-10 -tp=.1 -tf=.5 -rampOrder=5 -show=oneCylStretched2IMEX44s.show -go=go > oneCylStretched2IMEX44s.out 
	$(cgins) -noplot oneCylMove -g=oneCylStretchedGride4.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.005  -tolFactor=1e-10 -tp=.1 -tf=.5 -rampOrder=5 -show=oneCylStretched4IMEX44s.show -go=go > oneCylStretched4IMEX44s.out 
	$(cgins) -noplot oneCylMove -g=oneCylStretchedGride8.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.0025 -tolFactor=1e-10 -tp=.1 -tf=.5 -rampOrder=5 -show=oneCylStretched8IMEX44s.show -go=go > oneCylStretched8IMEX44s.out 



# --------- PAPER ------
# IMEX44, nu=.05, dt=[.05,.025,.0125,.00625]
oneCylMoveIMEX44:
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede2.order4.hdf -ts=imex -orderInTime=4 -nu=.05 -dtMax=.05   -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl2IMEX44.show -go=go > oneCyl2IMEX44.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede4.order4.hdf -ts=imex -orderInTime=4 -nu=.05 -dtMax=.025  -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl4IMEX44.show -go=go > oneCyl4IMEX44.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede8.order4.hdf -ts=imex -orderInTime=4 -nu=.05 -dtMax=.0125 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl8IMEX44.show -go=go > oneCyl8IMEX44.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede16.order4.hdf -ts=imex -orderInTime=4 -nu=.05 -dtMax=.00625 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl16IMEX44.show -go=go > oneCyl16IMEX44.out 

# IMEX24, nu=.05, dt=[.05,.025,.0125,.00625]
oneCylMoveIMEX24:
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede2.order4.hdf -ts=imex -orderInTime=2 -nu=.05 -dtMax=.05   -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl2IMEX24.show -go=go > oneCyl2IMEX24.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede4.order4.hdf -ts=imex -orderInTime=2 -nu=.05 -dtMax=.025  -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl4IMEX24.show -go=go > oneCyl4IMEX24.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede8.order4.hdf -ts=imex -orderInTime=2 -nu=.05 -dtMax=.0125 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl8IMEX24.show -go=go > oneCyl8IMEX24.out
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede16.order4.hdf -ts=imex -orderInTime=2 -nu=.05 -dtMax=.00625 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl16IMEX24.show -go=go > oneCyl16IMEX24.out

# IMEX22, nu=.05, dt=[.05,.025,.0125,.00625]
oneCylMoveIMEX22:
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede2.order2.hdf -ts=imex -orderInTime=2 -nu=.05 -dtMax=.05   -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl2IMEX22.show -go=go > oneCyl2IMEX22.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede4.order2.hdf -ts=imex -orderInTime=2 -nu=.05 -dtMax=.025  -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl4IMEX22.show -go=go > oneCyl4IMEX22.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede8.order2.hdf -ts=imex -orderInTime=2 -nu=.05 -dtMax=.0125 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl8IMEX22.show -go=go > oneCyl8IMEX22.out
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede16.order2.hdf -ts=imex -orderInTime=2 -nu=.05 -dtMax=.00625 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl16IMEX22.show -go=go > oneCyl16IMEX22.out

# PC44, nu=.05, dt=[.05,.025,.0125,.00625]
oneCylMovePC44:
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede2.order4.hdf -ts=pc4 -orderInTime=4 -nu=.05 -dtMax=.05   -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl2PC44.show -go=go > oneCyl2PC44.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede4.order4.hdf -ts=pc4 -orderInTime=4 -nu=.05 -dtMax=.025  -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl4PC44.show -go=go > oneCyl4PC44.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede8.order4.hdf -ts=pc4 -orderInTime=4 -nu=.05 -dtMax=.0125 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl8PC44.show -go=go > oneCyl8PC44.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede16.order4.hdf -ts=pc4 -orderInTime=4 -nu=.05 -dtMax=.00625 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl16PC44.show -go=go > oneCyl16PC44.out 


# PC22, nu=.05, dt=[.05,.025,.0125,.00625]
oneCylMovePC22:
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede2.order2.hdf -ts=pc -orderInTime=2 -nu=.05 -dtMax=.05   -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl2PC22.show -go=go > oneCyl2PC22.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede4.order2.hdf -ts=pc -orderInTime=2 -nu=.05 -dtMax=.025  -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl4PC22.show -go=go > oneCyl4PC22.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede8.order2.hdf -ts=pc -orderInTime=2 -nu=.05 -dtMax=.0125 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl8PC22.show -go=go > oneCyl8PC22.out
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede16.order2.hdf -ts=pc -orderInTime=2 -nu=.05 -dtMax=.00625 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl16PC22.show -go=go > oneCyl16PC22.out

# run G16 
oneCylMoveFine:
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede16.order4.hdf -ts=imex -orderInTime=2 -nu=.05 -dtMax=.00625 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl16IMEX24.show -go=go > oneCyl16IMEX24.out
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede16.order4.hdf -ts=imex -orderInTime=4 -nu=.05 -dtMax=.00625 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl16IMEX44.show -go=go > oneCyl16IMEX44.out 
	$(cgins) -noplot oneCylMove -g=oneCylGridFixede16.order4.hdf -ts=pc4 -orderInTime=4 -nu=.05 -dtMax=.00625 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl16PC44.show -go=go > oneCyl16PC44.out 


# rampOrder=3
oneCylMoveIMEX44d:
	$(cgins) -noplot oneCylMove -g=oneCylGride2.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.01   -tolFactor=1e-10 -tp=.1 -tf=.5 -rampOrder=3 -show=oneCyl2IMEX44d.show -go=go > oneCyl2IMEX44d.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride4.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.005  -tolFactor=1e-10 -tp=.1 -tf=.5 -rampOrder=3 -show=oneCyl4IMEX44d.show -go=go > oneCyl4IMEX44d.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride8.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.0025 -tolFactor=1e-10 -tp=.1 -tf=.5 -rampOrder=3 -show=oneCyl8IMEX44d.show -go=go > oneCyl8IMEX44d.out 


# rampOrder=4 
oneCylMoveIMEX44b:
	$(cgins) -noplot oneCylMove -g=oneCylGride2.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.01   -tolFactor=1e-10 -tp=.1 -tf=.5 -rampOrder=4 -show=oneCyl2IMEX44b.show -go=go > oneCyl2IMEX44b.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride4.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.005  -tolFactor=1e-10 -tp=.1 -tf=.5 -rampOrder=4 -show=oneCyl4IMEX44b.show -go=go > oneCyl4IMEX44b.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride8.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.0025 -tolFactor=1e-10 -tp=.1 -tf=.5 -rampOrder=4 -show=oneCyl8IMEX44b.show -go=go > oneCyl8IMEX44b.out 

# rampOrder=5 dt=[.04,.02,.05] Bigger dt 
# rampOrder=5 dt=[.05,.025,.0125] Bigger dt OK
oneCylMoveIMEX44ab:
	$(cgins) -noplot oneCylMove -g=oneCylGride2.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.05   -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl2IMEX44ab.show -go=go > oneCyl2IMEX44ab.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride4.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.025  -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl4IMEX44ab.show -go=go > oneCyl4IMEX44ab.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride8.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.0125 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl8IMEX44ab.show -go=go > oneCyl8IMEX44ab.out 

# rampOrder=5 dt=[.02,.01,.005] Bigger dt 
oneCylMoveIMEX44aa:
	$(cgins) -noplot oneCylMove -g=oneCylGride2.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.02  -tolFactor=1e-10 -tp=.1 -tf=1 -rampOrder=5 -show=oneCyl2IMEX44aa.show -go=go > oneCyl2IMEX44aa.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride4.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.01  -tolFactor=1e-10 -tp=.1 -tf=1 -rampOrder=5 -show=oneCyl4IMEX44aa.show -go=go > oneCyl4IMEX44aa.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride8.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.005 -tolFactor=1e-10 -tp=.1 -tf=1 -rampOrder=5 -show=oneCyl8IMEX44aa.show -go=go > oneCyl8IMEX44aa.out 

# rampOrder=5 dt=[.01,.005,.0025]
oneCylMoveIMEX44a:
	$(cgins) -noplot oneCylMove -g=oneCylGride2.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.01   -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl2IMEX44a.show -go=go > oneCyl2IMEX44a.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride4.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.005  -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl4IMEX44a.show -go=go > oneCyl4IMEX44a.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride8.order4.hdf -ts=imex -orderInTime=4 -nu=.1 -dtMax=.0025 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl8IMEX44a.show -go=go > oneCyl8IMEX44a.out 

 



# IMEX24a bigger dt  
oneCylMoveIMEX24a:
	$(cgins) -noplot oneCylMove -g=oneCylGride2.order4.hdf -ts=imex -orderInTime=2 -nu=.1 -dtMax=.02  -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl2IMEX24a.show -go=go > oneCyl2IMEX24a.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride4.order4.hdf -ts=imex -orderInTime=2 -nu=.1 -dtMax=.01  -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl4IMEX24a.show -go=go > oneCyl4IMEX24a.out 
	$(cgins) -noplot oneCylMove -g=oneCylGride8.order4.hdf -ts=imex -orderInTime=2 -nu=.1 -dtMax=.005 -tolFactor=1e-10 -tp=.1 -tf=2 -rampOrder=5 -show=oneCyl8IMEX24a.show -go=go > oneCyl8IMEX24a.out 



# One cylinder, PC44
oneCylPC44n0p1:
	$(cgins) -noplot tcilc -g=oneCylFixede2.order4 -nu=.1 -inflow=rampedPressure -inflowPresssure=1 -cpn=0 -dtMax=.01 -tp=.5 -tf=1.0 -project=0 -solver=yale -psolver=yale -ts=pc4 -orderInTime=4 -show=oneCyl2PC44nu0p1.show -go=go > oneCyl2PC44nu0p1.out 
	$(cgins) -noplot tcilc -g=oneCylFixede4.order4 -nu=.1 -inflow=rampedPressure -inflowPresssure=1 -cpn=0 -dtMax=.005 -tp=.5 -tf=1.0 -project=0 -solver=yale -psolver=yale -ts=pc4 -orderInTime=4 -show=oneCyl4PC44nu0p1.show -go=go > oneCyl4PC44nu0p1.out 
	$(cgins) -noplot tcilc -g=oneCylFixede8.order4 -nu=.1 -inflow=rampedPressure -inflowPresssure=1 -cpn=0 -dtMax=.0025 -tp=.5 -tf=1.0 -project=0 -solver=yale -psolver=yale -ts=pc4 -orderInTime=4 -show=oneCyl8PC44nu0p1.show -go=go > oneCyl8PC44nu0p1.out 



# One cylinder, IMEX44 **CHANGE RAMP INTERVAL 
oneCylIMEX44n0p1:
	$(cgins) -noplot tcilc -g=oneCylFixede2.order4 -nu=.1 -inflow=rampedPressure -inflowPresssure=1 -cpn=0 -dtMax=.01 -tp=.5 -tf=2.0 -rampOrder=4 -tb=4 -project=0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=oneCyl2IMEX44nu0p1.show -go=go > oneCyl2IMEX44nu0p1.out 
	$(cgins) -noplot tcilc -g=oneCylFixede4.order4 -nu=.1 -inflow=rampedPressure -inflowPresssure=1 -cpn=0 -dtMax=.005 -tp=.5 -tf=2.0 -rampOrder=4 -tb=4 -project=0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=oneCyl4IMEX44nu0p1.show -go=go > oneCyl4IMEX44nu0p1.out 
	$(cgins) -noplot tcilc -g=oneCylFixede8.order4 -nu=.1 -inflow=rampedPressure -inflowPresssure=1 -cpn=0 -dtMax=.0025 -tp=.5 -tf=2.0 -rampOrder=4 -tb=4 -project=0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=oneCyl8IMEX44nu0p1.show -go=go > oneCyl8IMEX44nu0p1.out 
	$(cgins) -noplot tcilc -g=oneCylFixede16.order4 -nu=.1 -inflow=rampedPressure -inflowPresssure=1 -cpn=0 -dtMax=.00125 -tp=.5 -tf=2.0 -rampOrder=4 -tb=4 -project=0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=oneCyl16IMEX44nu0p1.show -go=go > oneCyl16IMEX44nu0p1.out 


# One cylinder, IMEX22 
oneCylIMEX22n0p1:
	$(cgins) -noplot tcilc -g=oneCylFixede2.order2 -nu=.1 -inflow=rampedPressure -inflowPresssure=1 -dtMax=.01 -tp=.5 -tf=1.0 -project=0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=oneCyl2IMEX2nu0p1.show -go=go > oneCyl2IMEX2nu0p1.out 
	$(cgins) -noplot tcilc -g=oneCylFixede4.order2 -nu=.1 -inflow=rampedPressure -inflowPresssure=1 -dtMax=.005 -tp=.5 -tf=1.0 -project=0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=oneCyl4IMEX2nu0p1.show -go=go > oneCyl4IMEX2nu0p1.out 
	$(cgins) -noplot tcilc -g=oneCylFixede8.order2 -nu=.1 -inflow=rampedPressure -inflowPresssure=1 -dtMax=.0025 -tp=.5 -tf=1.0 -project=0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=oneCyl8IMEX2nu0p1.show -go=go > oneCyl8IMEX2nu0p1.out 
	$(cgins) -noplot tcilc -g=oneCylFixede16.order2 -nu=.1 -inflow=rampedPressure -inflowPresssure=1 -dtMax=.00125 -tp=.5 -tf=1.0 -project=0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=oneCyl16IMEX2nu0p1.show -go=go > oneCyl16IMEX2nu0p1.out 


#  IMEX44 Bewno, nu=1e-4  LONGER RUN FOR THE PAPER
# Note add AD for coarse grid to stabilize
# NOTE: Use iterative solvers for G32
caseIMEX44nu0001:
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede2.order4 -nu=.0001 -inflow=ramp -dtMax=.01 -tp=.5 -tf=10 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -ad4=1 -ad41=.2 -ad42=.2 -ad21=.5 -ad22=.5 -show=twoCyl2IMEX44nu0001.show -go=go >! twoCyl2IMEX44nu0001.out   
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede4.order4 -nu=.0001 -inflow=ramp -dtMax=.005 -tp=.5 -tf=10 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl4IMEX44nu0001.show -go=go >! twoCyl4IMEX44nu0001.out   
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede8.order4 -nu=.0001 -inflow=ramp -dtMax=.0025 -tp=.5 -tf=10 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl8IMEX44nu0001.show -go=go >! twoCyl8IMEX44nu0001.out   
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede16.order4 -nu=.0001 -inflow=ramp -dtMax=.00125 -tp=.5 -tf=10 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl16IMEX44nu0001.show -go=go >! twoCyl16IMEX44nu0001.out 
#	$(cgins) -noplot tcilc -g=twoCylindersGridFixede32.order4 -nu=.0001 -inflow=ramp -dtMax=.000625 -tp=0.5 -tf=8 -tolFactor=1e-8  -iluLevels=10 -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl32IMEX44nu0001.show -debug=3 -go=go >! twoCyl32IMEX44nu0001.out


#  IMEX22 + AD2, nu=1e-4 - LONGER RUN FOR THE PAPER
caseIMEX22nu0001:
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede2.order2 -nu=.0001 -inflow=ramp -dtMax=.01 -tp=.5 -tf=10 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=.25 -ad22=.25 -show=twoCyl2IMEX22nu0001.show -go=go >! twoCyl2IMEX22nu0001.out   
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede4.order2 -nu=.0001 -inflow=ramp -dtMax=.005 -tp=.5 -tf=10 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=.25 -ad22=.25 -show=twoCyl4IMEX22nu0001.show -go=go >! twoCyl4IMEX22nu0001.out   
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede8.order2 -nu=.0001 -inflow=ramp -dtMax=.0025 -tp=.5 -tf=10 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=.25 -ad22=.25 -show=twoCyl8IMEX22nu0001.show -go=go >! twoCyl8IMEX22nu0001.out 
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede16.order2 -nu=.0001 -inflow=ramp -dtMax=.00125 -tp=.5 -tf=10 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=.25 -ad22=.25 -show=twoCyl16IMEX22nu0001.show -go=go >! twoCyl16IMEX22nu0001.out 
#	$(cgins) -noplot tcilc -g=twoCylindersGridFixede32.order2 -nu=.0001 -inflow=ramp -dtMax=.000625 -tp=.5 -tf=10 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=.25 -ad22=.25 -show=twoCyl32IMEX22nu0001.show -go=go >! twoCyl32IMEX22nu0001.out 



# Test convergence: IMEX44 ramped p TURN OFF BEWNO 
caseIMEX44nu0p1nb:
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede1.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.02 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=twoCyl1IMEX44nu0p1nb.show -go=go >! twoCyl1IMEX44nu0p1nb.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede2.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.01 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=twoCyl2IMEX44nu0p1nb.show -go=go >! twoCyl2IMEX44nu0p1nb.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede4.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.005 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=twoCyl4IMEX44nu0p1nb.show -go=go >! twoCyl4IMEX44nu0p1nb.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede8.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.0025 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=twoCyl8IMEX44nu0p1nb.show -go=go >! twoCyl8IMEX44nu0p1nb.out



# Test convergence: IMEX44 ramped p (rampOrder=5 for pressure gives almost the same results as rampOrder=4)
# Reduce dtMax -- Question : how often is dt recmputed? -- almost the same 
caseIMEX44nu0p1a:
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede1.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.01 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl1IMEX44nu0p1a.show -go=go >! twoCyl1IMEX44nu0p1a.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede2.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.005 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl2IMEX44nu0p1a.show -go=go >! twoCyl2IMEX44nu0p1a.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede4.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.0025 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl4IMEX44nu0p1a.show -go=go >! twoCyl4IMEX44nu0p1a.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede8.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.00125 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl8IMEX44nu0p1a.show -go=go >! twoCyl8IMEX44nu0p1a.out

caseIMEX44nu0p1:
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede1.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.02 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl1IMEX44nu0p1.show -go=go >! twoCyl1IMEX44nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede2.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.01 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl2IMEX44nu0p1.show -go=go >! twoCyl2IMEX44nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede4.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.005 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl4IMEX44nu0p1.show -go=go >! twoCyl4IMEX44nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede8.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.0025 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl8IMEX44nu0p1.show -go=go >! twoCyl8IMEX44nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede16.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.00125 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl16IMEX44nu0p1.show -go=go >! twoCyl16IMEX44nu0p1.out

# Test convergence: PC44, ramped p
casePC44nu0p1:
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede1.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.02 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=pc4 -orderInTime=4 -ao=bweno -show=twoCyl1PC44nu0p1.show -go=go >! twoCyl1PC44nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede2.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.01 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=pc4 -orderInTime=4 -ao=bweno -show=twoCyl2PC44nu0p1.show -go=go >! twoCyl2PC44nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede4.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.005 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=pc4 -orderInTime=4 -ao=bweno -show=twoCyl4PC44nu0p1.show -go=go >! twoCyl4PC44nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede8.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.0025 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=pc4 -orderInTime=4 -ao=bweno -show=twoCyl8PC44nu0p1.show -go=go >! twoCyl8PC44nu0p1.out


# Test convergence: IMEX24 ramped p (rampOrder=5 for pressure gives almost the same results as rampOrder=4)
caseIMEX24nu0p1:
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede1.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.02 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=2 -ao=bweno -show=twoCyl1IMEX24nu0p1.show -go=go >! twoCyl1IMEX24nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede2.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.01 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=2 -ao=bweno -show=twoCyl2IMEX24nu0p1.show -go=go >! twoCyl2IMEX24nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede4.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.005 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=2 -ao=bweno -show=twoCyl4IMEX24nu0p1.show -go=go >! twoCyl4IMEX24nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede8.order4 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -rampOrder=4 -dtMax=.0025 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -orderInTime=2 -ao=bweno -show=twoCyl8IMEX24nu0p1.show -go=go >! twoCyl8IMEX24nu0p1.out

# Test convergence: IMEX44 ramped v 
caseIMEX44nu0p1v:
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede1.order4 -nu=.1 -inflow=ramp -dtMax=.02 -tp=.5 -tf=0.5 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl1IMEX44nu0p1v.show -go=go >! twoCyl1IMEX44nu0p1v.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede2.order4 -nu=.1 -inflow=ramp -dtMax=.01 -tp=.5 -tf=0.5 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl2IMEX44nu0p1v.show -go=go >! twoCyl2IMEX44nu0p1v.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede4.order4 -nu=.1 -inflow=ramp -dtMax=.005 -tp=.5 -tf=0.5 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl4IMEX44nu0p1v.show -go=go >! twoCyl4IMEX44nu0p1v.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede8.order4 -nu=.1 -inflow=ramp -dtMax=.0025 -tp=.5 -tf=0.5 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -ao=bweno -show=twoCyl8IMEX44nu0p1v.show -go=go >! twoCyl8IMEX44nu0p1v.out


# Test convergence: IMEX22 + AD2 ramped p
caseIMEX22nu0p1:
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede1.order2 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -dtMax=.02 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=twoCyl1IMEX2nu0p1.show -go=go >! twoCyl1IMEX2nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede2.order2 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -dtMax=.01 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=twoCyl2IMEX2nu0p1.show -go=go >! twoCyl2IMEX2nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede4.order2 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -dtMax=.005 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=twoCyl4IMEX2nu0p1.show -go=go >! twoCyl4IMEX2nu0p1.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede8.order2 -nu=.1 -inflow=rampedPressure -inflowPressure=2 -dtMax=.0025 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=twoCyl8IMEX2nu0p1.show -go=go >! twoCyl8IMEX2nu0p1.out

# Test convergence: IMEX22 + AD2 ramped v
caseO2p1v:
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede1.order2 -nu=.1 -inflow=ramp -dtMax=.02 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=twoCyl1IMEX2nu0p1v.show -go=go >! twoCyl1IMEX2nu0p1v.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede2.order2 -nu=.1 -inflow=ramp -dtMax=.01 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=twoCyl2IMEX2nu0p1v.show -go=go >! twoCyl2IMEX2nu0p1v.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede4.order2 -nu=.1 -inflow=ramp -dtMax=.005 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=twoCyl4IMEX2nu0p1v.show -go=go >! twoCyl4IMEX2nu0p1v.out
	$(cgins) -noplot tcilc -g=twoCylindersGridFixede8.order2 -nu=.1 -inflow=ramp -dtMax=.0025 -tp=.5 -tf=1.0 -solver=yale -psolver=yale -ts=imex -ad2=1 -ad21=0.2 -ad22=0.2 -show=twoCyl8IMEX2nu0p1v.show -go=go >! twoCyl8IMEX2nu0p1v.out


# nu=.01  - Bweno
case01b :

# nu=.005, No bweno 
case005 :
	$(cgins) -noplot tcilc -g=tcilce2.order4 -nu=.005 -inflow=rampedPressure -inflowPressure=1 -tp=.5 -tf=2 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=tcilc2IMEX4nu005.show -go=go >! tcilc2IMEX4nu005.out 
	$(cgins) -noplot tcilc -g=tcilce4.order4 -nu=.005 -inflow=rampedPressure -inflowPressure=1 -tp=.5 -tf=2 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=tcilc4IMEX4nu005.show -go=go >! tcilc4IMEX4nu005.out 
	$(cgins) -noplot tcilc -g=tcilce8.order4 -nu=.005 -inflow=rampedPressure -inflowPressure=1 -tp=.5 -tf=2 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=tcilc8IMEX4nu005.show -go=go >! tcilc8IMEX4nu005.out 


# nu=.001, No bweno 
case001 :
	$(cgins) -noplot tcilc -g=tcilce2.order4 -nu=.001 -inflow=rampedPressure -inflowPressure=1 -tp=.5 -tf=.5 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=tcilc2IMEX4nu001.show -go=go >! tcilc2IMEX4nu001.out 
	$(cgins) -noplot tcilc -g=tcilce4.order4 -nu=.001 -inflow=rampedPressure -inflowPressure=1 -tp=.5 -tf=.5 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=tcilc4IMEX4nu001.show -go=go >! tcilc4IMEX4nu001.out 
	$(cgins) -noplot tcilc -g=tcilce8.order4 -nu=.001 -inflow=rampedPressure -inflowPressure=1 -tp=.5 -tf=.5 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=tcilc8IMEX4nu001.show -go=go >! tcilc8IMEX8nu001.out 
	$(cgins) -noplot tcilc -g=tcilce16.order4 -nu=.001 -inflow=rampedPressure -inflowPressure=1 -tp=.5 -tf=.5 -solver=yale -psolver=yale -ts=imex -orderInTime=4 -show=tcilc16IMEX4nu001.show -go=go >! tcilc16IMEX4nu001.out 
