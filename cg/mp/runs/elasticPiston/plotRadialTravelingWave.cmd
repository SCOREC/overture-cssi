#
#  plotStuff plotRadialTravelingWave.cmd -show=radialTravelingWaveG8.show -name=radialTravelingWaveG8
# 
#
$show="radialTravelingWave.show"; $vMax=""; $solution=1; $name="radialTravelingWave";
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"vMax=f"=>\$vMax,"solution=i"=>\$solution  );
#
$show
#
solution: $solution
frame series:fluidDomain
* 
contour
  plot:p
  vertical scale factor 0.
  plot contour lines (toggle)
  if( $vMax ne "" ){ $cmd="min max 0 $vMax"; }else{ $cmd="#"; }
  $cmd
  exit
frame series:solidDomain
* 
derived types
  stressNorm
exit
contour
  adjust grid for displacement 1
  plot:stressNorm
  vertical scale factor 0.
  if( $vMax ne "" ){ $cmd="min max 0 $vMax"; }else{ $cmd="#"; }
  $cmd
  plot contour lines (toggle)
exit
#
#
DISPLAY COLOUR BAR:0 0
DISPLAY AXES:0 0
#
solution: 11
hardcopy file name:0 radialTravelingWaveG8t0p1stress.ps
hardcopy save:0
#
solution: 21
hardcopy file name:0 radialTravelingWaveG8t0p2stress.ps
hardcopy save:0
#
solution: 31
hardcopy file name:0 radialTravelingWaveG8t0p3stress.ps
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