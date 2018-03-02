#
#  Plot results from Cgmx
#
#    plotStuff plotSphere.cmd -show=sphereInABoxG2.show
#    plotStuff plotSphere.cmd -show=sphereInABoxG4.show
#
$show="sphereInABoxG2.show"; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show );
# -------------------------------------------------------------------------------------------------
$show
#
DISPLAY AXES:0 0
# DISPLAY COLOUR BAR:0 0
contour
  plot the grid
    toggle grid 0 0
    grid colour 1 BRASS
    grid colour 2 BRASS
    grid colour 3 BRASS
    plot block boundaries 0
    exit this menu
  exit
previous
hardcopy vertical resolution:0 1024
hardcopy horizontal resolution:0 1024
set view:0 -0.0909366 -0.0530211 0 1.0216 0.866025 0.17101 -0.469846 0 0.939693 0.34202 0.5 -0.296198 0.813798
#
plot:Ex
hardcopy file name:0 sphereInABoxExt0p5.ps
hardcopy save:0
#
plot:Ey
hardcopy file name:0 sphereInABoxEyt0p5.ps
hardcopy save:0