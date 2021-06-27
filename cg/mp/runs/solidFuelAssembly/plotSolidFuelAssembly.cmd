#
# plotStuff plotSolidFuelAssembly.cmd -show=sfaOneRod.show -name=sfaOneRod -solution=11 -Tmin=0 -Tmax=.9
# plotStuff plotSolidFuelAssembly.cmd -show=sfa5Rods.show -name=sfa5Rods -solution=11 -Tmin=0 -Tmax=.9
# plotStuff plotSolidFuelAssembly.cmd -show=sfa19Rods.show -name=sfa19Rods -solution=11 -Tmin=0 -Tmax=.9
# plotStuff plotSolidFuelAssembly.cmd -show=sfa37Rods.show -name=sfa37Rods -solution=11 -Tmin=0 -Tmax=.9
# plotStuff plotSolidFuelAssembly.cmd -show=sfa37RodsH10.show -name=sfa37RodsH10 -solution=11 -Tmin=0 -Tmax=.9 -view=1 
#
$Tmin=0.; $Tmax=.2; 
# $show="sfa.show";
# $show="sfa1.show";
# $show="sfa3.show";
# $show="sfa3lp5.show";
# $show="sfa3l2.show";
# $show="sfa3l5.show";
# 
# $show="sfa3l1f2.show"; $Tmax=.2;  $name="solidFuelAssembly3pinsl1"; 
# $show="sfa3l4f2.show"; $Tmax=.48; $name="solidFuelAssembly3pinsl4"; 
# $show="sfa3l4f2t0p5.show"; $Tmax=.48; $name="sfa3l4f2t0p5"; 
# $show="sfal2f1.show";  $Tmax=-1.; $name="sfal2f1";
$show="sfa3l4f1.show";  $Tmax=-1.; $name="sfa3l4f1"; $view=0; 
# get command line arguments
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "opt=s"=>\$opt, "Tmin=f"=>\$Tmin, "Tmax=f"=>\$Tmax, "solution=f"=>\$solution, "view=i"=>\$view );
#
$show
solution: $solution
# 
set view:0 0.00483384 -0.00483384 0 1.30521 0.866025 0.17101 -0.469846 -0.5 0.296198 -0.813798 5.55112e-17 0.939693 0.34202
if( $view == 1 ){ $cmd="set view:0 -0.036939 -0.0646901 0 0.915549 0.866025 0.17101 -0.469846 -0.5 0.296198 -0.813798 5.55112e-17 0.939693 0.34202"; }else{ $cmd="#"; }
$cmd
#
frame series:fluidDomain 
contour 
  plot:T
  min max $Tmin $Tmax
  contour lines 0
  exit
 #  
frame series:rodDomain
contour
  min max  $Tmin $Tmax
  contour lines 0
  exit
grid
  plot shaded surfaces (3D) 0
  plot block boundaries 0
exit this menu
#
frame series:hexagonalContainer
contour
  min max  $Tmin $Tmax
  contour lines 0
  exit
grid
  plot shaded surfaces (3D) 0
  coarsening factor 4
  exit this menu
pause
  DISPLAY COLOUR BAR:0 0
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  DISPLAY LABELS:0 0
$pname = $name . "T.ps"; 
hardcopy file name:0 $pname
hardcopy save:0

x-r:0
x-r:0
x-r:0
x-r:0
x-r:0
x-r:0
x-r:0
x-r:0
x-r:0
x-r:0
x+r:0
y+r:0
y+r:0
y+r:0
x+r:0
x+r:0
bigger:0
set view:0 0.00483384 -0.00483384 0 1.30521 0.866025 0.17101 -0.469846 -0.5 0.296198 -0.813798 5.55112e-17 0.939693 0.34202


x-r 90
y-r 20
x+r 50
# 
previous
frame series:fluidDomain
contour



# bigger
  DISPLAY COLOUR BAR:0 0
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  DISPLAY LABELS:0 0
# 
if( $Tmax > $Tmin ){ $minMax ="min max $Tmin $Tmax"; }else{ $minMax="*"; }
#
previous


frame series:fluidDomain
contour
  contour lines 0
  $minMax
  plot:T
  exit
frame series:rodDomain
contour
  contour lines 0
  $minMax
  exit
frame series:hexagonalContainer
contour
  contour lines 0
  $minMax
  exit
grid
  plot shaded surfaces (3D) 0
  plot grid lines 0
  colour block boundaries black
  exit this menu
frame series:rodDomain
grid
  plot grid lines 0
  plot shaded surfaces (3D) 0
  colour block boundaries black
  exit this menu
# 
pause
$pname = $name . "T.ps"; 
hardcopy file name:0 $pname
hardcopy save:0
frame series:fluidDomain
plot:w
$pname = $name . "TandW.ps"; 
hardcopy file name:0 $pname
hardcopy save:0
# 
#* --- line plots ---
#
$num=201; $z0=.0001; $z1=.5; $z2=2.; $z3=3.999; 
$lines ="4 $num\n" . \
        "-2.5 0. $z0 2.5 0 $z0\n" .\
        "-2.5 0. $z1 2.5 0 $z1\n" .\
        "-2.5 0. $z2 2.5 0 $z2\n" .\
        "-2.5 0. $z3 2.5 0 $z3\n";
$lines .= "T_line_0\n";
$lines .= "add T_line_1\n";
$lines .= "add T_line_2\n";
$lines .= "add T_line_3\n";
$lines .= "add x0_line_0\n";
$lines .= "add x0_line_1\n";
$lines .= "add x0_line_2\n";
$lines .= "add x0_line_3\n";
frame series:fluidDomain
contour
  plot:T
  line plots
    set bogus value for points not interpolated
      -1234. 
    specify lines
     $lines
# 
      add w_line_0
      add w_line_1
      add w_line_2
      add w_line_3
      save results to a matlab file
      $mname = $name . "Fluid.m"; 
      $mname
      exit this menu
    exit this menu
  exit
frame series:hexagonalContainer
contour
  line plots
    set bogus value for points not interpolated
      -1234. 
    specify lines
     $lines
      save results to a matlab file
      $mname = $name . "Duct.m"; 
      $mname
      exit this menu
    exit this menu
  exit
frame series:rodDomain
contour
  line plots
    set bogus value for points not interpolated
      -1234. 
    specify lines
     $lines
      save results to a matlab file
      $mname = $name . "Pins.m"; 
      $mname
      exit this menu
    exit this menu
  exit
