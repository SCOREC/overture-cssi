#
# Plot results from the shock hitting a beam
#   plotStuff plotBeamInAChannel -show=<s> -var=[r|p|s|rcs]
# 
# rcs= reverse colour schlieren]
#
#   plotStuff plotBeamInAChannel -show=beam8.show
#
#   plotStuff plotBeamInAChannel -show=beam32.show -var=p -solution=8 -name=shockBeamPressureT0p7
#   plotStuff plotBeamInAChannel -show=beam32.show -var=p -solution=14 -name=shockBeamPressureT1p3
#
#   plotStuff plotBeamInAChannel -show=beam32.show -var=s -vMin=.4 -vMax=1 -solution=8 -res=2048 -name=shockBeamSchlierenT0p7
#   plotStuff plotBeamInAChannel -show=beam32.show -var=s -vMin=.4 -vMax=1 -solution=14 -res=2048 -name=shockBeamSchlierenT1p3
#
#   plotStuff plotBeamInAChannel -show=beam32.show -var=rcs -vMin=-1 -vMax=-.1 -solution=8 -res=2048 -name=shockBeamColourSchlierenT0p7
#   plotStuff plotBeamInAChannel -show=beam32.show -var=rcs -vMin=-1 -vMax=-.1 -solution=14 -res=2048 -name=shockBeamColourSchlierenT1p3
#
$show="beam4.show"; $name="beam4"; $root=""; $vMin=0.; $vMax=-1.;  $res=1024; $solution=2; 
# $root = "dE_M1p5"; $root = "dELin_M1p5";
$var="p"; # component names
* ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"var=s"=>\$var,"solution=i"=>\$solution,"vMin=f"=>\$vMin,"vMax=f"=>\$vMax,\
            "res=i"=>\$res,"solution=i"=>\$solution );
if( $var eq "r" ){ $plotVar="rho"; }
if( $var eq "p" ){ $plotVar="p"; }
if( $var eq "s" ){ $plotVar="schlieren"; }
if( $var eq "rcs" ){ $plotVar="schlieren"; } # reverse colour schlieren
#
$show
#
derived types
  schlieren
  # flip sclieren scale for blue background and 
  if( $var eq "rcs" ){ $cmd="schlieren parameters\n -1. 15."; }else{ $cmd="#"; }
  $cmd
exit
solution: $solution
plot:$plotVar
# plot:p
contour
  vertical scale factor 0.
  coarsening factor 1 (<0 : adaptive)
  # schlieren: 
  plot contour lines (toggle)
  if( $var eq "s" ){ $cmd="gray\n min max .4 .1"; }else{ $cmd="#"; }
  $cmd 
  if( $vMin < $vMax ){ $cmd="min max $vMin $vMax"; }else{ $cmd="#"; }
  $cmd 
  # p: 
  # min max .5 1.0001 
exit
#
DISPLAY AXES:0 0
hardcopy vertical resolution:0 $res
hardcopy horizontal resolution:0 $res
if( $res eq 2048 ){ $cmd="line width scale factor:0 2"; }else{ $cmd="#"; }
$cmd
 # plot beam centerline:
forcing regions
exit
x-
#
pause
$plotName = $name . ".ps"; 
hardcopy file name:0 $plotName
hardcopy save:0

# 
bigger
DISPLAY COLOUR BAR:0 0
DISPLAY AXES:0 0
#
# -- MOVIE --
#
DISPLAY LABELS:0 0
solution: 1
save movie files 1
movie file name:beam32
pause
show movie