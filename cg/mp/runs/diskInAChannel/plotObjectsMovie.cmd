#
#  plotStuff plotObjectsMovie.cmd -show=fiveBodiesG4m.show -name=fiveBodiesG4 -nDomains=5 -vorMin=-3 -vorMax=3 -stride=2 
#
#  plotStuff plotObjectsMovie.cmd -show=fiveBodiesG8m.show -name=fiveBodiesG8 -nDomains=5 -vorMin=-10 -vorMax=10 -stride=2 
#  plotStuff plotObjectsMovie.cmd -show=fiveBodiesG8m.show -name=fiveBodiesG8 -nDomains=5 -plot=sl -stride=2 -slMax=1.5
# 
#  plotStuff plotObjectsMovie.cmd -show=fiveBodiesG8a.show -name=fiveBodiesG8 -nDomains=5 -vorMin=-10 -vorMax=10 -stride=2 
#
#  plotStuff plotObjectsMovie.cmd -show=fiveBodiesG4m.show -name=fiveBodiesG4SL -nDomains=5 -plot=sl -stride=2 
#
$show="fiveBodiesG4.show"; $vMax=""; $vorMin=-.5; $vorMax=.50001; $solution=1; $name="fiveBodiesG4";  $tp=""; $nDomains=5; 
$plot="contour"; $stride=1; $slMax=.75; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"vorMin=f"=>\$vorMin,"vorMax=f"=>\$vorMax,"vMax=f"=>\$vMax,"slMax=f"=>\$slMax,\
            "solution=i"=>\$solution,"tp=s"=>\$tp,"nDomains=i"=>\$nDomains,"stride=i"=>\$stride,"plot=s"=>\$plot );
#
$show
#
# solution: $solution
previous
#
frame series:fluidDomain
# 
derived types
  vorticity
exit
contour
  plot:vorticity
  if( $vorMax ne "" ){ $cmd="min max $vorMin $vorMax"; }else{ $cmd="#"; } 
  $cmd
  plot contour lines (toggle)
  vertical scale factor 0.
exit
# 
if( $plot eq "sl" ){ $cmd = "erase\n stream lines\n min max 0. $slMax\n streamline density 80\n exit"; }else{ $cmd="#"; }
$cmd
#
# --- plot solid domains ---
#
if( $vMax ne "" ){ $vcmd="min max 0 $vMax"; }else{ $vcmd="#"; } \
$cmd="#";
for( $d=1; $d <= $nDomains; $d++ ){\
$cmd .= "\n frame series:solidDomain$d\n derived types\n displacementNorm\n exit\n contour\n adjust grid for displacement 1\n plot:displacementNorm \n vertical scale factor 0.\n $vcmd\n  plot contour lines (toggle)\n exit"; }
$cmd
#
DISPLAY COLOUR BAR:0 0
DISPLAY AXES:0 0
pause
DISPLAY LABELS:0 0
set view:0 0.00906344 -0.00906344 0 1.26336 1 0 0 0 1 0 0 0 1
#
movie file name: $name
save movie files 1
solution: 1
stride: $stride
if( $name =~ /.*G8.*/ ){ $cmd = "set view:0 0.141315 -0.0134633 0 1.98305 1 0 0 0 1 0 0 0 1"; }else{ $cmd="#"; }
$cmd

show movie


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