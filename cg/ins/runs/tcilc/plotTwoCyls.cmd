#
#  plotStuff plotTwoCyls.cmd -show=tcilc4IMEX4.show -name=tcilc4 -vorMin=-10 -vorMax=10
#  plotStuff plotTwoCyls.cmd -show=tcilc8IMEX4.show -name=tcilc8 -vorMin=-20 -vorMax=20
#
#   plotStuff plotTwoCyls.cmd -show=twoCylG4O4.show -name=twoCylG4O4 -vorMin=-50 -vorMax=50
#   plotStuff plotTwoCyls.cmd -show=twoCylG4O4a.show -name=twoCylG4O4 -vorMin=-50 -vorMax=50
#   plotStuff plotTwoCyls.cmd -show=twoCylG4O4b.show -name=twoCylG4O4 -vorMin=-50 -vorMax=50
# 
#   plotStuff plotTwoCyls.cmd -show=twoCylG8O4.show -name=twoCylG8O4 -vorMin=-50 -vorMax=50
#
#  Under-resolved:
#   plotStuff plotTwoCyls.cmd -show=tcilc2IMEX4nu001.show -name=tcilc2IMEX -vorMin=-20 -vorMax=20
#   plotStuff plotTwoCyls.cmd -show=tcilc4IMEX4nu001.show -name=tcilc4IMEX -vorMin=-20 -vorMax=20
#   plotStuff plotTwoCyls.cmd -show=tcilc8IMEX4nu001.show -name=tcilc8IMEX -vorMin=-20 -vorMax=20
# 
# -- slow start
#   plotStuff plotTwoCyls.cmd -show=tcilc2IMEX4nu001ss.show -name=tcilc2IMEXss -vorMin=-20 -vorMax=20
#   plotStuff plotTwoCyls.cmd -show=tcilc4IMEX4nu001ss.show -name=tcilc4IMEXss -vorMin=-20 -vorMax=20
#   plotStuff plotTwoCyls.cmd -show=tcilc8IMEX4nu001.show -name=tcilc8IMEXss -vorMin=-20 -vorMax=20
# 
#  -- bewno
#   plotStuff plotTwoCyls.cmd -show=tcilc2IMEX4nu001b.show -name=tcilc2IMEXb -vorMin=-20 -vorMax=20
#   plotStuff plotTwoCyls.cmd -show=tcilc4IMEX4nu001b.show -name=tcilc4IMEXb -vorMin=-20 -vorMax=20
#   plotStuff plotTwoCyls.cmd -show=tcilc8IMEX4nu001b.show -name=tcilc8IMEXb -vorMin=-20 -vorMax=20
#   plotStuff plotTwoCyls.cmd -show=tcilc16IMEX4nu001b.show -name=tcilc16IMEXb -vorMin=-20 -vorMax=20
#
#  -- order=2
#   plotStuff plotTwoCyls.cmd -show=tcilc2IMEX2nu001.show -name=tcilc2IMEX2 -vorMin=-20 -vorMax=20
#   plotStuff plotTwoCyls.cmd -show=tcilc4IMEX2nu001.show -name=tcilc4IMEX2 -vorMin=-20 -vorMax=20
#   plotStuff plotTwoCyls.cmd -show=tcilc8IMEX2nu001.show -name=tcilc8IMEX2 -vorMin=-20 -vorMax=20
# 
# -- Smaller nu: nu=.0001
#   plotStuff plotTwoCyls.cmd -show=tcilc2IMEX4nu0001b.show -name=tcilc2IMEXc -vorMin=-30 -vorMax=30
#   plotStuff plotTwoCyls.cmd -show=tcilc4IMEX4nu0001b.show -name=tcilc4IMEXc -vorMin=-30 -vorMax=30
#   plotStuff plotTwoCyls.cmd -show=tcilc8IMEX4nu0001b.show -name=tcilc8IMEXc -vorMin=-30 -vorMax=30
#   plotStuff plotTwoCyls.cmd -show=tcilc16IMEX4nu0001b.show -name=tcilc16IMEXc -vorMin=-30 -vorMax=30
#
# For paper: 
#   plotStuff plotTwoCyls.cmd -show=twoCyl2IMEX44nu0001.show -name=twoCyl2IMEX44 -vorMin=-70 -vorMax=50
#   plotStuff plotTwoCyls.cmd -show=twoCyl4IMEX44nu0001.show -name=twoCyl4IMEX44 -vorMin=-70 -vorMax=50
#   plotStuff plotTwoCyls.cmd -show=twoCyl8IMEX44nu0001.show -name=twoCyl8IMEX44 -vorMin=-70 -vorMax=50
#   plotStuff plotTwoCyls.cmd -show=twoCyl16IMEX44nu0001.show -name=twoCyl16IMEX44 -vorMin=-70 -vorMax=50
#
#   plotStuff plotTwoCyls.cmd -show=twoCyl32IMEX44nu0001.show -name=twoCyl32IMEX44 -vorMin=-70 -vorMax=50
#
#   -- order=2 : 
#   plotStuff plotTwoCyls.cmd -show=twoCyl2IMEX22nu0001.show -name=twoCyl2IMEX22 -vorMin=-70 -vorMax=50
#   plotStuff plotTwoCyls.cmd -show=twoCyl4IMEX22nu0001.show -name=twoCyl4IMEX22 -vorMin=-70 -vorMax=50
#   plotStuff plotTwoCyls.cmd -show=twoCyl8IMEX22nu0001.show -name=twoCyl8IMEX22 -vorMin=-70 -vorMax=50
#   plotStuff plotTwoCyls.cmd -show=twoCyl16IMEX22nu0001.show -name=twoCyl16IMEX22 -vorMin=-70 -vorMax=50
#   plotStuff plotTwoCyls.cmd -show=twoCyl32IMEX22nu0001.show -name=twoCyl32IMEX22 -vorMin=-70 -vorMax=50
# 
$show="fiveBodiesG4.show"; $vMax=""; $vorMin=-.5; $vorMax=.50001; $solution=-1; $name="fiveBodiesG4";  $tp=""; $nDomains=5; 
# $pMin="-10"; $pMax="5";
$pMin="-7"; $pMax="4";
$sMin="0"; $sMax="3"; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"vorMin=f"=>\$vorMin,"vorMax=f"=>\$vorMax,"vMax=f"=>\$vMax,\
            "solution=i"=>\$solution,"tp=s"=>\$tp,"nDomains=i"=>\$nDomains,"pMin=f"=>\$pMin,"pMax=f"=>\$pMax,\
            "sMin=f"=>\$sMin,"sMax=f"=>\$sMax  );
#
$show
#
previous
# 
derived types
  vorticity
  speed
exit
contour
  plot:speed
  if( $sMax ne "" ){ $cmd="min max $sMin $sMax"; }else{ $cmd="#"; } 
  $cmd
  plot:p
  if( $pMax ne "" ){ $cmd="min max $pMin $pMax"; }else{ $cmd="#"; } 
  $cmd
  # 
  plot:vorticity
  if( $vorMax ne "" ){ $cmd="min max $vorMin $vorMax"; }else{ $cmd="#"; } 
  $cmd
  plot contour lines (toggle)
  vertical scale factor 0.
  coarsening factor 1 (<0 : adaptive)
exit
# pause
DISPLAY AXES:0 0
# set view:0 -0.0583647 -0.0115004 0 1.43487 1 0 0 0 1 0 0 0 1
set view:0 -0.253808 -0.00269676 0 1.43172 1 0 0 0 1 0 0 0 1
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
# ------ Just plot vorticity
#
solution: 7
$plotName= $name . "t3p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
solution: 11
plot:vorticity
$plotName= $name . "t5p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
solution: 16
plot:vorticity
$plotName= $name . "t7p5.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
exit


# ------ plot vorticity, pressure and speed
#
solution: 7
$plotName= $name . "t3p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
plot:p
$plotName= $name . "Pressuret3p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
plot:speed
$plotName= $name . "Speedt3p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
solution: 11
plot:vorticity
$plotName= $name . "t5p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
plot:p
$plotName= $name . "Pressuret5p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
plot:speed
$plotName= $name . "Speedt5p0.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
solution: 16
plot:vorticity
$plotName= $name . "t7p5.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
# 
plot:speed
$plotName= $name . "Speedt7p5.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0
#
plot:p
$plotName= $name . "Pressuret7p5.ps"; 
hardcopy file name:0 $plotName
hardcopy save:0






# ZOOM NEAR INTERPOLATION BOUNDARY 
solution: 16 
contour
  plot:u
  wire frame (toggle)
  exit
hardcopy vertical resolution:0 2048
hardcopy horizontal resolution:0 2048
line width scale factor:0 3
set view:0 0.0686027 -0.0103694 0 2.86197 1 0 0 0 1 0 0 0 1
bigger:0
smaller:0
hardcopy save:0




















pause
#
DISPLAY COLOUR BAR:0 0
DISPLAY AXES:0 0
DISPLAY LABELS:0 0
set view:0 0 0.00564972 0 1.30467 1 0 0 0 1 0 0 0 1
pause
#
movie file name: $name
save movie files 1
solution: 1

show movie
