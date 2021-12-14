#
# -- annulus 
#   plotStuff plotAnnulus -show=annulusTTG2 -name=annulusTTG2 -solution=7
#   plotStuff plotAnnulus -show=annulusTTG4 -name=annulusTTG4 -solution=7
#   plotStuff plotAnnulus -show=annulusTTG8 -name=annulusTTG8 -solution=7
#   plotStuff plotAnnulus -show=annulusTTG16 -name=annulusTTG16 -solution=7
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
  # speed
  # stressNorm
  displacementNorm
exit
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
  # grid colour 1 NAVYBLUE
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


