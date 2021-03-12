# 
#   plotStuff plotGrid.cmd -show=twoSquaresInterfaceGride2.order4 -name=twoSquaresInterfaceGrid
#
#
$show="ellipticalDiskGride8.order2.hdf"; $opt="#"; 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "opt=s"=>\$opt );
#
$show
#
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048  
  line width scale factor:0 3
  DISPLAY SQUARES:0 0
  DISPLAY AXES:0 0
  bigger:0
  colour boundaries by grid number
  #
  $opt
  #
  coarsening factor 1
  plot
  $plotName = $name . ".ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0


pause
#
  set view:0 0.187103 -0.191008 0 12.3143 1 0 0 0 1 0 0 0 1
  $plotName = $name . "Zoom.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
pause
