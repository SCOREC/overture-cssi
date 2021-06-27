#
# plot results from the 3D heatExchanger
#
#   plotStuff plotSolution3d.cmd -show=hexExplicit.show -name=hexExplicitTt0p1 -Tmin=0 -Tmax=10 -solution=2
#   plotStuff plotSolution3d.cmd -show=hexExplicit.show -name=hexExplicitTt0p2 -Tmin=0 -Tmax=10 -solution=3
# 
#
$show="he.show"; $solution="-1"; $name="plot"; $field="T"; $Tmin=0; $Tmax=10;
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field,"Tmin=f"=>\$Tmin,"Tmax=f"=>\$Tmax );
#
$show
#
solution: $solution
frame series:solidDomain
contour
 contour lines 0
 min max $Tmin $Tmax
exit
#
frame series:tubeDomain
contour
  contour lines 0
  plot:T
  min max $Tmin $Tmax
exit
grid
  plot block boundaries 0
  plot shaded surfaces (3D) 0
  coarsening factor 2
exit this menu
DISPLAY AXES:0 0
set view:0 -0.0422716 -0.0595288 0 1.10114 0.766044 0.321394 -0.55667 0 0.866025 0.5 0.642788 -0.383022 0.663414
pause
# 
$plotName = "$name.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0



contour



  plot:eFieldNorm
  delete contour plane 0
  contour lines 0
  min max 0 $emax
exit
