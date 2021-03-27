#
#  plotStuff plotInterface.cmd -show=interfaceG8.show -name=interfaceG8
#  plotStuff plotInterface.cmd -show=twoSquareInterfaceG8Rotated45.show -name=interfaceG8Angle45
#
#  plotStuff plotInterface.cmd -show=twoSquaresInterfaceG8.show -name=twoSquaresInterfaceG8Mla2Mla3 -solution=11
#
#  plotStuff plotInterface.cmd -show=oneDiskMLA.show -name=oneDiskMLA3 -solution=4 -time=t0p3 
#  plotStuff plotInterface.cmd -show=oneDiskMLA.show -name=oneDiskMLA3 -solution=11 -time=t1p0 
#
$solution=11; $time="t1p0"; 
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"solution=i"=>\$solution,"time=s"=>\$time );
#
$show
#
solution: $solution
contour
  plot:Ex
  vertical scale factor 0.3
exit
DISPLAY AXES:0 0
set view:0 -0.05 -0.0060423 0 1.04416 0.642788 -0.383022 0.663414 0.766044 0.321394 -0.55667 -1.0631e-16 0.866025 0.5
# ondeDiskMLA:
set view:0 -0.098565 -0.0113293 0 0.832704 0.939693 0.219846 -0.262003 -0.34202 0.604023 -0.719846 4.74648e-17 0.766044 0.642788
pause
#
#
plot:Ey
$plotName = $name . "Ey$time.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
plot:Ey error
$plotName = $name . "EyErr$time.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
plot:Py
$plotName = $name . "Py$time.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0



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
