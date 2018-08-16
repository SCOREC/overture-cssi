#
#  plotStuff plotInterface.cmd -show=interfaceG8.show -name=interfaceG8
#  plotStuff plotInterface.cmd -show=twoSquareInterfaceG8Rotated45.show -name=interfaceG8Angle45
#  
GetOptions( "show=s"=>\$show,"name=s"=>\$name );
#
$show
#
solution: 11
contour
  plot:Ex
exit
#
DISPLAY AXES:0 0
#
plot:Hz
$plotName = $name . "Hzt1p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
plot:Ex
$plotName = $name . "Ext1p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
#
plot:Ey
$plotName = $name . "Eyt1p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
plot:Py
$plotName = $name . "Pyt1p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0