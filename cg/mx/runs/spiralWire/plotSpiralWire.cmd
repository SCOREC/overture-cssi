#
#    plotStuff plotSpiralWire -show=spiralWire.show -name=spiralWire
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

