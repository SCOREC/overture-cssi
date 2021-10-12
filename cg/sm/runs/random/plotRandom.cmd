# 
#   plotStuff plotRandom -show=randomSquare32tttt -name=randomSquaretttt
#   plotStuff plotRandom -show=randomSquare32tdtt -name=randomSquaretdtt
#   plotStuff plotRandom -show=randomSquare32dtdt -name=randomSquaredtdt
#
#   plotStuff plotRandom -show=randomDiskG4t -name=randomDiskG4t
#   plotStuff plotRandom -show=randomDiskG4d -name=randomDiskG4d
#
$show="diffract16.show"; $min=0.; $max=-1.; $grids=1;  $umin=0.; $umax=-1.; $solution=1; 
GetOptions( "show=s"=>\$show, "name=s"=>\$name, "min=f"=>\$min,"max=f"=>\$max, "umin=f"=>\$umin,"umax=f"=>\$umax,\
            "solution=i"=>\$solution );
#
$show
# 
# --- save sequence info ----
plot sequence:solutionNorms
  uNorm
  add vNorm
  add pNorm
  add divU
  save results to a matlab file
    $matlab = $name . ".m"; 
    $matlab
  exit this menu

solution: $solution
derived types
specify displacement components
  0 1 2 
speed
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