#
#  plotStuff plotDisk -show=diskTractionG4 -name=diskTractionG4 -solution=7
#  plotStuff plotDisk -show=diskTractionG8 -name=diskTractionG8 -solution=7
#  plotStuff plotDisk -show=diskTractionG16 -name=diskTractionG16 -solution=7
#
# -- annulus 
#   plotStuff plotDisk -show=annulusTTG4 -name=annulusTTG4 -solution=7
#   plotStuff plotDisk -show=annulusTTG8 -name=annulusTTG8 -solution=7
#
$show="diffract16.show"; $min=0.; $max=-1.; $grids=1;  $umin=0.; $umax=-1.; $solution=1; 
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "min=f"=>\$min,"max=f"=>\$max, "umin=f"=>\$umin,"umax=f"=>\$umax,\
            "solution=i"=>\$solution );
#
$show
# 
solution: $solution
derived types
specify displacement components
  0 1 2 
speed
# stressNorm
displacementNorm
exit
# ---- displacement ----
displacement
  # plot grid lines 0
  # plot non-physical boundaries 1
  # colour boundaries by refinement level number
  # raise the grid by this amount (2D) 0.1
  displacement scale factor 0.075
  coarsening factor 1
  colour grid lines from chosen name  
  grid colour 0 NAVYBLUE
  grid colour 1 NAVYBLUE
exit this menu
DISPLAY SQUARES:0 0
DISPLAY AXES:0 0
line width scale factor:0 3
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
# draw again to get dsf correct
displacement
exit 
#
$plotName = $name . "DisplacementMesh.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0

#
#   ------ CONTOURS -------------
contour
  plot:displacementNorm
  # optionally set the min and max
  if( $max > $min ){ $cmd ="min max $min $max"; }else{ $cmd="#"; }
  $cmd 
  if( $umax > $umin ){ $cmd ="plot:displacementNorm\n min max $umin $umax"; }else{ $cmd="#"; }
  $cmd 
  vertical scale factor 0.
  coarsening factor 1 
exit
x-
x-
DISPLAY AXES:0 0
$plotName = $name . "DisplacementNorm.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
pause
#
plot:uError
$plotName = $name . "u1Error.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#  
plot:pError
$plotName = $name . "pError.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0  

# ---- displacement ----
displacement
  # plot grid lines 0
  # plot non-physical boundaries 1
  # colour boundaries by refinement level number
  # raise the grid by this amount (2D) 0.1
  displacement scale factor 0.075
  coarsening factor 1
  colour grid lines from chosen name  
  grid colour 0 NAVYBLUE
  grid colour 1 NAVYBLUE
exit this menu
DISPLAY SQUARES:0 0
DISPLAY AXES:0 0
line width scale factor:0 3
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
# draw again to get dsf correct
displacement
exit 
#
$plotName = $name . "DisplacementMesh.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0




# 
#



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
