#
#    plotStuff plotSpiralWire -show=spiralWire.show -name=spiralWire
#    plotStuff plotSpiralWire -show=sprialWire.show -name=spiralWire
#
$show="spiralWire.show"; $name=""; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name );
# -------------------------------------------------------------------------------------------------
$show
previous
grid
  toggle grid 0 0
  plot grid lines 0
  plot block boundaries 0
  grid colour 1 BRASS
  grid colour 2 BRASS
  grid colour 3 BRASS
exit this menu
# 
derived types
E field norm
exit
plot:eFieldNorm
contour
  contour lines 0
exit


y+r:0
y+r:0
x+r:0
x+r:0
x+r:0
contour
  pick to delete contour planes
  delete contour plane 2
  delete contour plane 0
  plot the grid
    exit this menu
  pick to add contour plane y
  add contour plane  0.00000e+00  1.00000e+00  0.00000e+00 -2.96339e-01 -1.93972e+00  1.66602e+00 
  pick to add contour plane z
  add contour plane  0.00000e+00  0.00000e+00  1.00000e+00 -2.48087e-01  5.45270e-01 -2.50156e+00 
  reset:0
  y+r:0
  y+r:0
  y+r:0
  y+r:0
  x+r:0
  x+r:0
  x+r:0
  x-r:0
  exit
previous
previous
contour
  base min/max on contour plane values
  exit
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
DISPLAY COLOUR BAR:0 0
hardcopy file name:0 spiralWireEfieldNorm.ps
hardcopy save:0
y+r:0
y+r:0
reset:0
y+r:0
y+r:0
y+r:0
y+r:0
y+r:0
y+r:0
y+r:0
y+r:0
y+r:0
y-r:0
y-r:0
x+r:0
x+r:0
hardcopy file name:0 spiralWireEfieldNormB.ps
hardcopy save:0
reset:0
y-r:0
x+r:0
x+r:0
y-r:0
y-:0
hardcopy file name:0 spiralWireEfieldNormC.ps
y-r:0
x-r:0
x-r:0
x-r:0
x-r:0
x+r:0
y-r:0
x+r:0
z-r:0
x-r:0
y+:0
y+:0
hardcopy file name:0 spiralWireEfieldNormD.ps
hardcopy save:0
