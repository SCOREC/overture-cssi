#
#
# For paper: 
#   plotStuff plotOneCyl.cmd -show=oneCyl4IMEX44.show -name=oneCyl4IMEX44 
#
$show="oneCyl4IMEX44.show"; $vMax=""; $vorMin=-.5; $vorMax=.50001; $solution=-1; $name="fiveBodiesG4";  $tp=""; $nDomains=5; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"vorMin=f"=>\$vorMin,"vorMax=f"=>\$vorMax,"vMax=f"=>\$vMax,\
            "solution=i"=>\$solution,"tp=s"=>\$tp,"nDomains=i"=>\$nDomains );
#
$show
#
previous
# 
derived types
  vorticity
exit
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 3
contour
  plot:vorticity
  # plot contour lines (toggle)
  vertical scale factor 0.
  DISPLAY AXES:0 0
exit
x-:0
solution: 11
hardcopy file name:0 oneCylVort1p0.ps
hardcopy save:0
pause
solution: 21
hardcopy file name:0 oneCylVort2p0.ps
hardcopy save:0
pause
erase
#
#
stream lines
  arrow size 0.05
  arrow size 3.000000e-02
exit
hardcopy file name:0 oneCylSL2p0.ps
hardcopy save:0
pause
solution: 11
hardcopy file name:0 oneCylSL1p0.ps
hardcopy save:0







  vertical scale factor 0.
  coarsening factor 1 (<0 : adaptive)
exit



pause
DISPLAY AXES:0 0
set view:0 -0.0583647 -0.0115004 0 1.43487 1 0 0 0 1 0 0 0 1
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
#
solution: 7
$plotName= $name . "t3p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
solution: 11
$plotName= $name . "t5p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0



pause
#
DISPLAY COLOUR BAR:0 0
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
set view:0 0 0.00564972 0 1.30467 1 0 0 0 1 0 0 0 1
pause
#
movie file name: $name
save movie files 1
solution: 1

show movie
