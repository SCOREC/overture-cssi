#
#  plotStuff plotDisk.cmd -show=diskInChannelScf10G2.show -name=diskInChannelScf10G2
#  plotStuff plotDisk.cmd -show=diskInChannelScf1G2.show -name=diskInChannelScf1G2 -solution=21
#
#  plotStuff plotDisk.cmd -show=diskInADeformingChannelScf50G2.show -name=diskInADeformingChannelScf50G2 -solution=21
#  plotStuff plotDisk.cmd -show=diskInADeformingChannelScf10G2 -name=diskInADeformingChannelScf10G2 -solution=41
#
# plotStuff plotDisk.cmd -show=ellipseG4Scf5 -name=ellipseG4Scf5 -solution=41
# plotStuff plotDisk.cmd -show=ellipseG4Scf2 -name=ellipseG4Scf2 -solution=41
#
# plotStuff plotDisk.cmd -show=ellipseG4LongScf5 -name=ellipseG4LongScf5 -solution=101
# plotStuff plotDisk.cmd -show=ellipseG4LongerScf10 -name=ellipseG4LongerScf10 -solution=301
#
#  plotStuff plotDisk.cmd -show=diskInAChannelScf10G8 -name=diskInAChannelScf10G8 -solution=21 -tp=1p0
#
#  plotStuff plotDisk.cmd -show=diskInAChannelScf10G2.show -name=diskInChannelScf10G2
#  plotStuff plotDisk.cmd -show=diskInAChannelScf10G4.show -name=diskInChannelScf10G4 -solution=16
#
$show="balloon.show"; $vMax=""; $solution=1; $name="balloon";  $tp=""; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"vMax=f"=>\$vMax,"solution=i"=>\$solution,"tp=s"=>\$tp );
#
$show
#
# solution: $solution
previous
#
frame series:fluidDomain
* 
stream lines
  exit
# contour
#   plot:p
#   vertical scale factor 0.
#   plot contour lines (toggle)
#   if( $vMax ne "" ){ $cmd="min max 0 $vMax"; }else{ $cmd="#"; }
#   $cmd
#   exit
#
frame series:solidDomain
* 
derived types
  stressNorm
  displacementNorm
exit
contour
  adjust grid for displacement 1
  # plot:stressNorm
  plot:displacementNorm
  vertical scale factor 0.
  if( $vMax ne "" ){ $cmd="min max 0 $vMax"; }else{ $cmd="#"; }
  $cmd
  plot contour lines (toggle)
exit
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