#
#  plotStuff plotEigen.cmd -show=squareEigen128.show -name=sqEig
#  
GetOptions( "show=s"=>\$show,"name=s"=>\$name );
#
$show
#
DISPLAY AXES:0 0
solution: 11
contour
  plot:Ex
exit
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
plot:Py
$plotName = $name . "Pyt1p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0