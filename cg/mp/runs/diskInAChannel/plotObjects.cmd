#
#  plotStuff plotObjects.cmd -show=fiveBodiesG4.show -name=fiveBodiesG4 -nDomains=5 -solution=41 -tp=2p0
#  plotStuff plotObjects.cmd -show=fiveBodiesG4.show -name=fiveBodiesG4 -nDomains=5 -solution=101 -tp=5p0
#
# plotStuff plotObjects.cmd -show=fiveBodiesG4d.show -name=fiveBodiesG4d -nDomains=5 -solution=41 -tp=2p0
# plotStuff plotObjects.cmd -show=fiveBodiesG4e.show -name=fiveBodiesG4e -nDomains=5 -solution=41 -tp=2p0
# plotStuff plotObjects.cmd -show=fiveBodiesG4f.show -name=fiveBodiesG4f -nDomains=5 -solution=41 -tp=2p0
#
# plotStuff plotObjects.cmd -show=fiveBodiesG8f.show -name=fiveBodiesG8f -nDomains=5 -solution=41 -tp=2p0
# plotStuff plotObjects.cmd -show=fiveBodiesG16f.show -name=fiveBodiesG16f -nDomains=5 -solution=41 -tp=2p0
#
# plotStuff plotObjects.cmd -show=fiveBodiesG8g.show -name=fiveBodiesG8g -nDomains=5 -solution=41 -tp=2p0
# plotStuff plotObjects.cmd -show=fiveBodiesG16g.show -name=fiveBodiesG16g -nDomains=5 -solution=41 -tp=2p0
#
# For paper:
#  plotStuff plotObjects.cmd -show=fiveBodiesG4f.show -name=fiveBodiesG4g -nDomains=5 -solution=61 -tp=6p0
#  plotStuff plotObjects.cmd -show=fiveBodiesG4f.show -name=fiveBodiesG4g -nDomains=5 -solution=81 -tp=8p0
#  plotStuff plotObjects.cmd -show=fiveBodiesG4f.show -name=fiveBodiesG4g -nDomains=5 -solution=101 -tp=10p0
#
#  plotStuff plotObjects.cmd -show=fiveBodiesG8g.show -name=fiveBodiesG8g -nDomains=5 -solution=81 -tp=8p0
#  plotStuff plotObjects.cmd -show=fiveBodiesG8g.show -name=fiveBodiesG8g -nDomains=5 -solution=101 -tp=10p0 -vMax=.13
# 
# -- larger nu
#     plotStuff plotObjects.cmd -show=fiveBodiesG4n0p1.show -name=fiveBodiesG4n0p1 -nDomains=5 -solution=41 -tp=2p0
#
$show="fiveBodiesG4.show"; $vMax=""; $solution=1; $name="fiveBodiesG4";  $tp=""; $nDomains=5; 
@uMax=(.15,.1,.2,.008,.09); 
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
  streamline density 100
  arrow size 0.02
  if( $vMax ne "" ){ $cmd="min max 0 $vMax"; }else{ $cmd="#"; }
  $cmd
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
$cmd="#";
for( $d=1; $d <= $nDomains; $d++ ){\
if( $uMax[$d-1] ne "" ){ $ucmd="min max 0 $uMax[$d-1]"; }else{ $ucmd="#"; } \
$cmd .= "\n frame series:solidDomain$d\n derived types\n displacementNorm\n exit\n contour\n adjust grid for displacement 1\n plot:displacementNorm \n vertical scale factor 0.\n $ucmd\n  plot contour lines (toggle)\n exit"; }
$cmd
# $cmd .= "\n frame series:solidDomain$d\n derived types\n stressNorm\n exit\n contour\n adjust grid for displacement 1\n plot:stressNorm \n vertical scale factor 0.\n $vcmd\n  plot contour lines (toggle)\n exit"; }
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