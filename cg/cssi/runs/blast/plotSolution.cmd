#
# Plot results from blast.cmd
#
#  plotStuff plotSolution.cmd -show=blastShapes.show -solution=9 -schlieren=rainbow -name=blastShapest0p4
#  plotStuff plotSolution.cmd -show=blastShapes.show -solution=21 -schlieren=rainbow -name=blastShapest1p0
#  plotStuff plotSolution.cmd -show=blastShapes.show -solution=9 -schlieren=gray -sMin=.2 -sMax=1 -name=blastShapes
#
#  plotStuff plotSolution.cmd -show=blastShapesG8.show -solution=9 -schlieren=rainbow -name=blastShapest0p4
#  plotStuff plotSolution.cmd -show=blastShapesG8.show -solution=21 -schlieren=rainbow -name=blastShapest1p0
#
#  plotStuff plotSolution.cmd -show=downTownl2r4.show -solution=9 -schlieren=rainbow -name=downTownSchlierenTimep08
#
# 
$show="blastShapes.show"; $sMin=.2; $sMax=-1.; $name=""; $solution=21; $time="10p0"; $schlieren="gray"
GetOptions( "show=s"=>\$show,"schlieren=s"=>\$schlieren,"name=s"=>\$name,"sMin=s"=>\$sMin,"sMax=s"=>\$sMax,"solution=i"=>\$solution,"time=s"=>\$time  );
#
$show
#
solution: $solution
# previous
# 
derived types
schlieren
   # flip sclieren scale for blue background and 
  if( $schlieren eq "rainbow" ){ $cmd="schlieren parameters\n -1. 15."; }else{ $cmd="#"; }
  $cmd
exit
DISPLAY SQUARES:0 0
DISPLAY COLOUR BAR:0 0
DISPLAY AXES:0 0
if( $show =~ /downTown/ ){ $cmd="DISPLAY LABELS:0 0\n set view:0 0.00309366 0 0 2.01397 1 0 0 0 1 0 0 0 1"; }else{ $cmd="#"; }
$cmd
contour
  plot:schlieren
  plot contour lines (toggle)
  # colour table:
  $schlieren
  # gray
  # rainbow
  vertical scale factor 0.
  coarsening factor 1 (<0 : adaptive)
  if( $sMin < $sMax ){ $cmd="min max $sMin $sMax"; }else{ $cmd="#"; }
  $cmd
exit
pause
$plotName = $name . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0


grid
  plot grid lines 0
  plot non-physical boundaries 1
  colour boundaries by refinement level number
exit this menu
