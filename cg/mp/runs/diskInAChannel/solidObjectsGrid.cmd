#
#   Solid objects in a fluid channel
#
#
# usage: ogen [noplot] solidObjectsGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i]
# 
# examples:
#     ogen -noplot solidObjectsGrid -interp=e -order=2 -factor=4 
#     ogen -noplot solidObjectsGrid -interp=e -order=2 -factor=6 
#     ogen -noplot solidObjectsGrid -interp=e -order=2 -factor=8
#
# -- longer domain
#     ogen -noplot solidObjectsGrid -interp=e -order=2 -factor=8 -xb=8 -prefix=fiveSolidsGrid
#
$prefix="solidObjectsGrid"; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-3.5; $xb=5; $ya=-2.; $yb=2.; $ml=0; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"prefix=s"=> \$prefix);
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
# --- build I-beam grid 
#
$count=0; 
$degree=3; # degree of nurbs
#
#   ---- I-Beam ----
#
$centerHeight=.5; $centerWidth=.35; $edgeHeight=.25;   $edgeWidth=1.; $angle=45; 
$cx=1.5; $cy=.75;
$coreRadius=.05; $xCore=0; $yCore=-.3; # offset core from center
# --- define curve ---
include solidIBeamCurve.h
#  --- build grids ---
include solidBodyGrids.h
# 
#  ---- SKIP : Square ---- 
#
$edgeHeight=1.; $edgeWidth=1.; $angle=0;
$coreRadius=.1; $xCore=.25; $yCore=0; # offset core from center
$cx=1.5; $cy=0;
#- include solidBeamCurve.h
#  --- build grids ---
#- include solidBodyGrids.h
# 
# ---- Beam -----
# 
$edgeHeight=.25; $edgeWidth=1.5; $angle=-10;
$coreRadius=.05; $xCore=.25; $yCore=0; # offset core from center
$cx=-1.5; $cy=.85;
include solidBeamCurve.h
#  --- build grids ---
include solidBodyGrids.h
# 
# ---- Beam -----
# 
$edgeHeight=.35; $edgeWidth=1.5; $angle=10;
$coreRadius=.05; $xCore=-.25; $yCore=0; # offset core from center
$cx=-1.5; $cy=-1.;
include solidBeamCurve.h
#  --- build grids ---
include solidBodyGrids.h
#
# ---- split-ring ---
#
$H=.6; $W=.5; $w=.25; $d=.2; $cx=0; $cy=0;  $angle=200;
$coreRadius=.05; $xCore=-.35; $yCore=-.35; # offset core from center
$cx=1.5; $cy=-.85;
include splitRingCurve.h
#  --- build grids ---
include solidBodyGrids.h
#
# ---- cross ---
#
$w=.15; $h=.75; $cx=0; $cy=0; $angle=200;
$coreRadius=.1; $xCore=0; $yCore=0; # core
include solidCrossCurve.h
#  --- build grids ---
include solidBodyGrids.h
#
## $gridNames = "iBeamExteriorGrid\n iBeamInnerBackGround\n solidCore\n iBeamInteriorGrid"; 
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
solidObjectsGrid
exit
