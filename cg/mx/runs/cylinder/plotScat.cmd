#
#  Plot results from Cgmx scattering from the crew reetnery vehichle
#
#    plotStuff plotScat.cmd -show=cicG4Order2.show
#    plotStuff plotScat.cmd -show=cicG4Order4.show -name=cicG4Order4
#
$show="cicG4Order2.show"; $name="cicG4Order2";
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show, "name=s"=>\$name );
# -------------------------------------------------------------------------------------------------
$show
#
DISPLAY AXES:0 0
# DISPLAY COLOUR BAR:0 0
previous
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  line width scale factor:0 3
#
contour
  plot:Ex
  $plotName=$name . "Ex.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
  plot:Ey
  $plotName=$name . "Ey.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
  plot:Ey error
  $plotName=$name . "EyError.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0
  plot:Ex error
  $plotName=$name . "ExError.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0


set view:0 -0.01251 0.000152542 0 1.09579 0.966822 -0.105335 0.232722 0.0244387 0.944987 0.326194 -0.254279 -0.309684 0.916208
#
grid
  toggle grid 0 0
  plot block boundaries 0
  coarsening factor 2
  grid colour 1 BRASS
  grid colour 2 BRASS
  grid colour 3 BRASS
  plot grid lines 0
exit this menu
# 
contour
  pick to delete contour planes
  delete contour plane 0
  contour lines 0
exit