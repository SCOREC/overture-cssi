#
#   plotStuff plotGrid.cmd -show=multiDiskGride2.order2.hdf -name=multiDiskGrid2
# 
#
$show="multiDiskGride2.order2.hdf";
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name );
#
$show
#
  DISPLAY SQUARES:0 0
  DISPLAY AXES:0 0
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048  
  bigger:0
  plot block boundaries 0
  colour grid lines by domain
  colour boundaries by domain
#
  hardcopy file name:0 multiDiskGrid2.ps
  hardcopy save:0
pause
#
  line width scale factor:0 3
  plot block boundaries 0
  set view:0 -0.173957 0.138729 0 6.40384 1 0 0 0 1 0 0 0 1
  hardcopy file name:0 multiDiskGrid2Zoom.ps
  hardcopy save:0
pause
  reset:0
  bigger:0
  bigger:0
  bigger:0
  bigger:0
  plot grid lines 0
  plot block boundaries 1
  hardcopy file name:0 multiDiskGridBlockBoundaries.ps
  hardcopy save:0  


  plot
  $plotName = $name . ".ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0

  set view:0 0.233167 -0.00650583 0 11.4097 1 0 0 0 1 0 0 0 1
  bigger:0
  smaller:0
  DISPLAY SQUARES:0 0
  $plotName = $name . "Zoom.ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0
