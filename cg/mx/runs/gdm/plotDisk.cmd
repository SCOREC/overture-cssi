#
#  plotStuff plotDisk.cmd -show=diskG8.show -name=diskBroadBand
#  
GetOptions( "show=s"=>\$show,"name=s"=>\$name );
#
$show
#
DISPLAY AXES:0 0
x-:0
#
contour
  plot:Ey
exit
#
solution: 1
$plotName = $name . "Eyt0p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
solution: 3
$plotName = $name . "Eyt0p2.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
solution: 5
$plotName = $name . "Eyt0p4.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
solution: 101
$plotName = $name . "Eyt10p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
