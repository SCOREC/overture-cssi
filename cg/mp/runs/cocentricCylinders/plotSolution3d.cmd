#
# plot results from the 3D cocentric cylinders
#
#   plotStuff plotSolution3d.cmd -show=cocentricCyl.show -name=cocentricCylTt1p0 -Tmin=1 -Tmax=5 -solution=11
#    
#
$show="he.show"; $solution="-1"; $name="plot"; $field="T"; $Tmin=0; $Tmax=-1;
$tSave=1; $numPerTime=2; $numToSave=5; # save solution at these time intervals
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "solution=i"=>\$solution,"tSave=f"=>\$tSave,\
      "numPerTime=i"=>\$numPerTime, "numToSave=i"=>\$numToSave,"field=s"=>\$field,"Tmin=f"=>\$Tmin,"Tmax=f"=>\$Tmax );
#
$show
#
solution: $solution
frame series:innerDomain
contour
 contour lines 0
 if( $Tmax > $Tmin ){ $cmd="min max $Tmin $Tmax"; }else{ $cmd="#"; }
 $cmd
exit
#
frame series:outerDomain
contour
  contour lines 0
  plot:T
  if( $Tmax > $Tmin ){ $cmd="min max $Tmin $Tmax"; }else{ $cmd="#"; }
  $cmd
exit
grid
  plot block boundaries 0
  plot shaded surfaces (3D) 0
  coarsening factor 2
exit this menu
DISPLAY AXES:0 0
set view:0 -0.0422716 -0.0595288 0 1.10114 0.766044 0.321394 -0.55667 0 0.866025 0.5 0.642788 -0.383022 0.663414
x-
pause
# 
$plotName = "$name.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0


