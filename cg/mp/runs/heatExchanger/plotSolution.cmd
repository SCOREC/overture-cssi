#
# plot results from the 2D heatExchanger
#
#   plotStuff plotSolution.cmd -show=hex2d.show -name=hex2dExplicitTt1p0 -Tmin=0 -Tmax=10 -solution=11
#   plotStuff plotSolution.cmd -show=hex2d.show -name=hex2dExplicitTt5p0 -Tmin=0 -Tmax=10 -solution=51
#
# plotStuff plotSolution.cmd -show=hex2dG16.show -name=hex2dG16ExplicitTt1p0 -Tmin=0 -Tmax=10 -solution=3
#
$show="he.show"; $solution="-1"; $name="plot"; $field="T"; $Tmin=0; $Tmax=10;
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field,"Tmin=f"=>\$Tmin,"Tmax=f"=>\$Tmax );
#
$show
#
#
solution: $solution
frame series:solidDomain
contour
 plot contour lines (toggle)
 min max $Tmin $Tmax
exit
#
frame series:tubeDomain
derived types
 vorticity
  speed
exit
contour
  plot contour lines (toggle)
  plot:T
  min max $Tmin $Tmax
exit
#
DISPLAY AXES:0 0
# set view:0 -0.0422716 -0.0595288 0 1.10114 0.766044 0.321394 -0.55667 0 0.866025 0.5 0.642788 -0.383022 0.663414
pause
# 
$plotName = "$name.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0