#
#  plotStuff plotObjects.cmd -show=fiveBodiesG4.show -name=fiveBodiesG4 -nDomains=5 -solution=41 -tp=2p0
#  plotStuff plotObjects.cmd -show=fiveBodiesG4.show -name=fiveBodiesG4 -nDomains=5 -solution=101 -tp=5p0
#
$show="fiveBodiesG4.show"; $vMax=""; $solution=1; $name="fiveBodiesG4";  $tp=""; $nDomains=5; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"vMax=f"=>\$vMax,"solution=i"=>\$solution,"tp=s"=>\$tp,"nDomains=i"=>\$nDomains );
#
$show
#
# solution: $solution
previous
#
frame series:fluidDomain
# 
stream lines
  streamline density 80
exit
#
# contour
#   plot:p
#   vertical scale factor 0.
#   plot contour lines (toggle)
#   if( $vMax ne "" ){ $cmd="min max 0 $vMax"; }else{ $cmd="#"; }
#   $cmd
#   exit
#
# --- plot solid domains ---
#
if( $vMax ne "" ){ $vcmd="min max 0 $vMax"; }else{ $vcmd="#"; } \
$cmd="#";
for( $d=1; $d <= $nDomains; $d++ ){\
$cmd .= "\n frame series:solidDomain$d\n derived types\n displacementNorm\n exit\n contour\n adjust grid for displacement 1\n plot:displacementNorm \n vertical scale factor 0.\n $vcmd\n  plot contour lines (toggle)\n exit"; }
$cmd
pause
#
DISPLAY COLOUR BAR:0 0
DISPLAY AXES:0 0
#
solution: $solution
$plotName = $name ."StreamLinesAndDisplacement$tp.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0

#
solution: 401
hardcopy file name:0 balloonG8t2p0stress.ps
hardcopy save:0
#
solution: 1801
hardcopy file name:0 balloonG8t9p0stress.ps
hardcopy save:0

pause
#
movie file name: $name
save movie files 1
show movie




# Line plots: 
contour
  line plots
    specify a boundary
    101 0 0 1
      add v1
      add v2
      add s11
      add s12
      add s22
      add u
      add v
      add x0
      add x1
      save results to a matlab file
      balloonInterface.m
      exit this menu
    exit this menu  
    clear all:0
    plot
exit



# -- plot speed ---
solution: $solution
frame series:fluidDomain
* 
derived types
  speed
exit
contour
  plot:speed
  vertical scale factor 0.
  plot contour lines (toggle)
  if( $vMax ne "" ){ $cmd="min max 0 $vMax"; }else{ $cmd="#"; }
  $cmd
  exit
frame series:solidDomain
* 
derived types
  speed
 specify velocity components
  0 1 2
exit
contour
  adjust grid for displacement 1
  * plot:vorz
  plot:speed
  vertical scale factor 0.
  if( $vMax ne "" ){ $cmd="min max 0 $vMax"; }else{ $cmd="#"; }
  $cmd
  plot contour lines (toggle)
exit
#
DISPLAY COLOUR BAR:0 0
DISPLAY AXES:0 0
#
pause
#
movie file name: $name
save movie files 1
show movie