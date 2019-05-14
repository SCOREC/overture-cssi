ogen = ${Overture}/bin/ogen
conv = ${CG}/common/src/conv.p


bubble:
	$(ogen) -noplot freeSurfaceDropGrid -interp=e -radx=1. -rady=1.25 -factor=2 
	$(ogen) -noplot freeSurfaceDropGrid -interp=e -radx=1. -rady=1.25 -factor=4
	$(ogen) -noplot freeSurfaceDropGrid -interp=e -radx=1. -rady=1. -prefix=oscillatingBubbleGrid -factor=1 
	$(ogen) -noplot freeSurfaceDropGrid -interp=e -radx=1. -rady=1. -prefix=oscillatingBubbleGrid -factor=2 
	$(ogen) -noplot freeSurfaceDropGrid -interp=e -radx=1. -rady=1. -prefix=oscillatingBubbleGrid -factor=4
	$(ogen) -noplot freeSurfaceDropGrid -interp=e -radx=1. -rady=1. -prefix=oscillatingBubbleGrid -factor=8
	$(ogen) -noplot freeSurfaceDropGrid -interp=e -radx=1. -rady=1. -prefix=oscillatingBubbleGrid -factor=16


obConv:
	$(conv) ../../conv/surfaceTension.conv -option=ob -ts=pc -tf=0.2 -psolver=best -ng=4 -caseNumber=1 -amp=0.00001 
	$(conv) ../../conv/surfaceTension.conv -option=ob -ts=pc -tf=0.2 -psolver=best -ng=4 -caseNumber=2 -amp=0.00001 
	$(conv) ../../conv/surfaceTension.conv -option=ob -ts=pc -tf=0.2 -psolver=best -ng=4 -caseNumber=3 -amp=0.00001 

ins:
	cd ../../; \
	make
