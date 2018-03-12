#
#  plotStuff plotEigen.cmd -show=diskEigen16.show -name=diskEig
#  plotStuff plotEigen.cmd -show=diskEigen16Np2.show -name=diskEig
#  
GetOptions( "show=s"=>\$show,"name=s"=>\$name );
#
$show
#
DISPLAY AXES:0 0
x-:0
solution: 11
contour
  plot:Ex
exit
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
plot:Ex error
$plotName = $name . "ExErrt1p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
# plot:Py
# $plotName = $name . "Pyt1p0.ps"; 
# hardcopy file name:0 $plotName
# hardcopy save:0
# 

