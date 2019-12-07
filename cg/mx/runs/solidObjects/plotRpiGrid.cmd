#
#   plotStuff plotRpiGrid.cmd -show=rpiGride4.order2.hdf -name=rpiGrid
# 
#
$show="solidObjectsGride4.order2.hdf";
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name );
#
$show
#
  colour grid lines from chosen name
  grid colour 0 BLUE
  grid colour 1 BLUE
  grid colour 2 BLUE
  grid colour 3 BLUE
#
#
  grid colour 4 RED
  grid colour 5 RED
# 
  grid colour 6 LIMEGREEN
  grid colour 7 LIMEGREEN
#
  grid colour 8 YELLOW
  grid colour 9 YELLOW
#
  bigger 
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
#
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048  
  line width scale factor:0 3
#
  plot
  $plotName = $name . ".ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
  pause
  # 
  # ZOOM
  set view:0 -0.0277331 -0.250927 0 4.28964 1 0 0 0 1 0 0 0 1
  $plotName = $name . "Zoom.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0










#
  grid colour 4 RED
  grid colour 5 RED
# 
  grid colour 6 MEDIUMORCHID
  grid colour 7 MEDIUMORCHID
#
  grid colour 8 SPRINGGREEN
  grid colour 9 SPRINGGREEN


  set view:0 -0.108471 -0.0058006 0 1.79677 1 0 0 0 1 0 0 0 1
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
#
  colour grid lines from chosen name
  grid colour 0 BLUE
  grid colour 1 BLUE
  grid colour 2 BLUE
  grid colour 3 BLUE
  grid colour 4 BLUE
  grid colour 5 BLUE
  grid colour 6 BLUE
#
#
  grid colour 7 ORANGERED
  grid colour 8 ORANGERED
# 
  grid colour 9 AQUAMARINE
  grid colour 10 AQUAMARINE
#
  grid colour 11 ORANGE
  grid colour 12 ORANGE
#
  grid colour 13 YELLOW
  grid colour 14 YELLOW
#
  grid colour 15 MEDIUMGOLDENROD
  grid colour 16 MEDIUMGOLDENROD
#
  grid colour 17 SPRINGGREEN
  grid colour 18 SPRINGGREEN
#
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048  
  line width scale factor:0 3
#
  plot
  $plotName = $name . ".ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
  # 
  # ZOOM
  set view:0 -0.362271 0.273521 0 3.53788 1 0 0 0 1 0 0 0 1
  $plotName = $name . "Zoom.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0

