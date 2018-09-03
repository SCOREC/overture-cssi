#
#  plotStuff plotGrid -show=objectsGride2.order2.hdf -name=objectsGride2
#  plotStuff plotGrid -show=objectsGride4.order2.hdf -name=objectsGride4
#
$show="diskInAChannelGride2.order2.hdf";
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name  );
#
$show
#
  # bigger
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  colour boundaries by chosen name
  colour grid lines from chosen name
  set view:0 -0.0602637 0.0489642 0 3.35864 1 0 0 0 1 0 0 0 1
#  grid colour 1 WHEAT
#  grid colour 1 BLACK
  grid colour 0 GREEN
  grid colour 1 BLUE
  grid colour 7 MEDIUMTURQUOISE
  grid colour 13 ORANGERED
  grid colour 9 MEDIUMORCHID
  grid colour 10 SEAGREEN
#
  set view:0 0.0133301 0.0132588 0 2.25128 1 0 0 0 1 0 0 0 1
#
  coarsening factor 4
#
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  line width scale factor:0 2
  plot interpolation points 0
  # colour interpolation points 1
#
  hardcopy file name:0 $name.ps
  hardcopy save:0
pause
#
#  Zoom
#
  coarsening factor 2
  plot
  set view:0 0.178801 0.00147887 0 9.40291 1 0 0 0 1 0 0 0 1
  $name = $name . "Zoom"; 
  hardcopy file name:0 $name.ps
  hardcopy save:0
