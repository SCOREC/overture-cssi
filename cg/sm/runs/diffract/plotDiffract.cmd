#
#  plotStuff plotDiffract -show=diffract8l2r2.show -min=0. -max=2.4 -umin=0. -umax=4.5 
# 
# -- SOS: 110418
#  plotStuff plotDiffract -show=diffract8l2r4f4.show -name=diffract8l2r4noGrids -min=0. -max=2.4 -umin=0. -umax=4.5 -grids=0 
# -- SOS: 100201: 
#  plotStuff plotDiffract -show=diffract8l2r4f4.show 
#  plotStuff plotDiffract -show=diffract8l2r2f4.show 
#  plotStuff plotDiffract -show=diffract8l2r4.show -name=diffract8l2r4noGrids -min=0. -max=2.4 -umin=0. -umax=4.5 -grids=0 
#  plotStuff plotDiffract -show=diffract8l2r4g.show -name=diffract8l2r4noGridsg -min=0. -max=2.4 -umin=0. -umax=4.5 -grids=0 
# 
# -- rerun SOS-C with slip-wall BC and filter4 : 
# plotStuff plotDiffract -show=diffract8l2r4cf4a.show -name=diffract8l2r4cf4 -min=0. -max=2.4 ** for smog 
#
# -- filter4 : 
# plotStuff plotDiffract -show=diffract8l2r4cf4.show -name=diffract8l2r4cf4noGrids -min=0. -max=2.4 -grids=0 [no amr grids shown]
# plotStuff plotDiffract -show=diffract8l2r4cf4.show -name=diffract8l2r4cf4 -min=0. -max=2.4
# 
# plotStuff plotDiffract -show=diffract8l2r4c.show -name=diffract8l2r4c
# plotStuff plotDiffract -show=diffract8l2r4g.show -name=diffract8l2r4noGridsg -min=0. -max=2.4 -grids=0 [no amr grids shown]
# plotStuff plotDiffract -show=diffract8l2r4g.show -name=diffract8l2r4g -min=0. -max=2.4          ** for smog 
# plotStuff plotDiffract -show=diffract8l2r4g.show -name=diffract8l2r4gEvolve           [variable min max]
# 
# plotStuff plotDiffract -show=diffract64g.show -name=diffract64g -min=0. -max=2.4             ** for smog 
# plotStuff plotDiffract -show=diffract32g.show -name=diffract32g
# plotStuff plotDiffract -show=diffract8g.show -name=diffract8g -min=0. -max=2.4               ** for smog 
# 
# plotStuff plotDiffract -show="diffract16l2r4.show"
# 
# plotStuff plotDiffract -show=diffract64c.show 
# plotStuff plotDiffract -show=diffract16c.show 
# plotStuff plotDiffract -show=diffract16g.show 
# plotStuff plotDiffract -show=diffract32c.show -name=diffract32c
# plotStuff plotDiffract -show=diffract32g.show -name=diffract32g
# 
# plotStuff plotDiffract -show=diffract8g.show -name=diffract8g
# plotStuff plotDiffract -show=diffract64.show
#
$show="diffract16.show"; $min=0.; $max=-1.; $grids=1;  $umin=0.; $umax=-1.;
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "min=f"=>\$min,"max=f"=>\$max, "umin=f"=>\$umin,"umax=f"=>\$umax,"grids=i"=>\$grids );
#
$show
# 
derived types
specify velocity components
if( $show =~ /g.show/ ){ $vComp = "0 1 2"; $uComp="6 7 8"; }else{ $vComp = "4 5 6";  $uComp="0 1 2";}  # displacement and velocity components 
  $vComp
specify displacement components
  $uComp
speed
stressNorm
displacementNorm
exit
DISPLAY SQUARES:0 0
# grid
displacement
  plot grid lines 0
  plot non-physical boundaries 1
  colour boundaries by refinement level number
  raise the grid by this amount (2D) 0.1
  displacement scale factor 0.075
exit this menu
if( $grids eq 0 ){ $cmd="erase"; }else{ $cmd="#"; }
$cmd
contour
  plot:speed
  # optionally set the min and max
  if( $max > $min ){ $cmd ="min max $min $max"; }else{ $cmd="#"; }
  $cmd 
  if( $umax > $umin ){ $cmd ="plot:displacementNorm\n min max $umin $umax"; }else{ $cmd="#"; }
  $cmd 
  vertical scale factor 0.
  compute coarsening factor 0
  coarsening factor 1 (<0 : adaptive)
  adjust grid for displacement 1
  displacement scale factor 0.075
  exit
# 
#
DISPLAY AXES:0 0
x-:0
x-:0
y-:0
bigger 1.1
line width scale factor:0 4
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
#
@sol = (1,5,6,7,9);
@time = ("0p0","0p8","1p0","1p2","1p6" );
$num=@sol;
$cmd="";
for( $i=0; $i<$num; $i++ ){\
$plotName = $name . "Speed$time[$i].ps"; \
$cmd .= "solution: $sol[$i]\n plot:speed\n hardcopy file name:0 $plotName\n hardcopy save:0\n"; \
$plotName = $name . "uNorm$time[$i].ps"; \
$cmd .= "plot:displacementNorm\n hardcopy file name:0 $plotName\n hardcopy save:0\n"; \
}
$cmd .= "#";
$cmd



hardcopy file name:0 $plotName
hardcopy save:0
plot:displacementNorm



solution 1
$plotName = $name . "Speed0p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0

# 
solution: 5
$plotName = $name . "Speed0p8.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
solution: 6
$plotName = $name . "Speed1p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
solution: 7
$plotName = $name . "Speed1p2.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
solution: 9
$plotName = $name . "Speed1p6.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0


contour
  plot:speed
  plot:u
  plot:v
  plot:div
  plot:vorz
  plot contour lines (toggle)
  plot:div
  plot:vorz
  plot:speed
  erase and exit
displacement
  displacement scale factor .1
  exit this menu
displacement
  coarsening factor 4
  coarsening factor 2
