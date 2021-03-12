#
#  plotStuff plotGrid.cmd -show=curvedChannel4.order2.hdf -name=curvedChannelGrid
#
echo to terminal 1
$show="curvedChannel4.order2.hdf"; $opt="#"; 
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
  #
  $opt
  #
  coarsening factor 1
  plot
  $plotName = $name . ".ps"; 
  hardcopy file name:0 $plotName
  hardcopy save:0

