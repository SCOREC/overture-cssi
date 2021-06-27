#
#   Solid objects in the shape of "R" "P" "I" in a fluid channel
#
#
# usage: ogen [noplot] rpiGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i]
# 
# examples:
#     ogen -noplot rpiGrid -interp=e -order=2 -factor=4 
#     ogen -noplot rpiGrid -interp=e -order=2 -factor=6 
#     ogen -noplot rpiGrid -interp=e -order=2 -factor=8
#
$prefix="rpiGrid"; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=0; $xb=8.5; $ya=-2.; $yb=2.; $ml=0; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name);
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.1/$factor;
#
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
#
create mappings
#
#  Background 
#
  rectangle
    mappingName
      backGround
   set corners
    $xa $xb $ya $yb
   lines
    $nx = int( ($xb-$xa)/$ds +1.5 ); 
    $ny = int( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
    boundary conditions
      1 2 3 4
  exit
#
$fluidGridNames="#";
$solidGridNames="#"; 
$solidDomains="#";
$shareInterface=99;
$bcInterface=99;
#
# I-beam, beam (upper), beam (lower), split-ring, cross
#
# --- build R-beam grid 
#
$count=0; 
$degree=3; # degree of nurbs
#
#   ---- R-Beam ----
#
$edgeWidth=0.5; $middleWidth=0.25; 
$slantLength=0.5; $bottomHeight=0.5; 
$topHeight=0.25;
$cx=2; $cy=$edgeWidth;
$coreRadius=.25; $xCore=.5*$middleWidth; $yCore=0.; # offset core from center
# --- define curve ---
include $ENV{CG}/mp/runs/diskInAChannel/solidRBeamCurve.h
#  --- build grids ---
include $ENV{CG}/mp/runs/diskInAChannel/solidBodyGrids.h
#
#   ---- P-Beam ----
#
$edgeWidth=0.5; $middleWidth=0.25; 
$slantLength=0.5; $bottomHeight=0.5; 
$topHeight=0.25;
$cx=4.25; $cy=$edgeWidth;
$coreRadius=.25; $xCore=.5*$middleWidth; $yCore=0.; # offset core from center
# --- define curve ---
include $ENV{CG}/mp/runs/diskInAChannel/solidPBeamCurve.h
#  --- build grids ---
include $ENV{CG}/mp/runs/diskInAChannel/solidBodyGrids.h
#
#   ---- I-Beam ----
#
$centerWidth=$edgeWidth; $edgeHeight=$edgeWidth;
$centerHeight=$height-2*$edgeHeight; 
$edgeWidth=$width; 
$angle=0; 
$cx=6.5; $cy=0;
$coreRadius=.15; $xCore=0; $yCore=.5*$centerHeight+.25*$edgeHeight; # offset core from center
# --- define curve ---
include $ENV{CG}/mp/runs/diskInAChannel/solidIBeamCurve.h
#  --- build grids ---
include $ENV{CG}/mp/runs/diskInAChannel/solidBodyGrids.h
# 
exit 
#
generate an overlapping grid
  backGround
  $fluidGridNames
  $solidGridNames
  done
  change parameters
    $solidDomains
    #
    specify a domain
     fluidDomain
       backGround
       $fluidGridNames
    done
    # choose implicit or explicit interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
    exit
  # pause
  # open graphics
  compute overlap
exit
#
save an overlapping grid
$name
rpiGrid
exit
