#
# Grid for a radial traveling waves FSI example (cgmp)
#
#   Notes: If ra=0 then the inner region has a background grid to remove the polar singularity 
#
# usage: ogen [noplot] radialTravelingWave -factor=<num> -order=[2/4/6/8] -interp=[e/i] 
# Options:
#    -ra,rb
# 
# examples:
#   ogen -noplot radialTravelingWaveGrid -interp=e -factor=2
#   ogen -noplot radialTravelingWaveGrid -interp=e -factor=4
#
#
$prefix="radialTravelingWaveGrid";
$order=2; $factor=1; $interp="e"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids"; $ml=0; 
$name=""; 
$R=1.; $Rbar=1.2;
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"ra=f"=> \$ra,"rb=f"=> \$rb,"rc=f"=> \$rc,\
            "interp=s"=> \$interp,"name=s"=> \$name,"per=i"=>\$per,"dsx=f"=>\$dsx,"prefix=s"=>\$prefix );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $per eq 1 ){ $suffix .= "p"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.1/$factor;
$pi = 4.*atan2(1.,1.);
# 
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
$bcInterface=100;  # bc for interfaces
$shareInterface=100;        # share value for interfaces
#
create mappings
#
# --- make a nurbs for a circle ----
#   **** Later we could make different curves for the interface ****
$radX=$R; $radY=$R;
include ellipseCurve.h
# 
# -- Make a hyperbolic grid --
#
  $nr = intmg( 5 + $order/2 );
  hyperbolic
    backward
    # Fixed radial distance
    $nDistInterface=.2;
    $nrm = intmg($nDistInterface/$ds+4.5);
    # $nDistInterface=($nr-3)*$ds;
    # $nrm=$nr-1; 
    distance to march $nDistInterface
    lines to march $nrm
    $nThetaInterface = intmg($arcLength/$ds+1.5);
    points on initial curve $nThetaInterface
    uniform dissipation 0.05
    volume smooths $numberOfVolumeSmooths
    equidistribution 0 (in [0,1])
    #
    # spacing: geometric
    # geometric stretch factor 1.05 
    #
    generate
    boundary conditions
      -1 -1 100 0 0 0
    share 
       0 0 100 0 0 0
    name fluidInterface
  exit
#
# --- inner domain as a hyperbolic mapping
  hyperbolic
    forward
    $nDist = $Rbar-$R;
    $nrm = intmg($nDist/$ds+4.5);
    distance to march $nDist
    lines to march $nrm
    points on initial curve $nThetaInterface
    uniform dissipation 0.05
    volume smooths $numberOfVolumeSmooths
    equidistribution 0. (in [0,1])
    #
    # spacing: geometric
    # geometric stretch factor 1.05 
    #
    generate
    boundary conditions
      -1 -1 100 4 0 0
    share 
       0 0 100 0 0 0
    name solidRegion
  exit
#
#  --- Inner fluid background grid ---
#
 rectangle
  $dsr=$ds*.8; # make spacing a litle smaller 
  $xb=$R-$nDistInterface+2*$ds; $xa=-$xb; $ya=$xa; $yb=$xb; 
  set corners
    $xa $xb $ya $yb
  lines
    $nx = int( ($xb-$xa)/$dsr +1.5 ); 
    $ny = int( ($yb-$ya)/$dsr +1.5 ); 
    $nx $ny
  boundary conditions
    0 0 0 0 
  share 
    0 0 0 0 
  mappingName
    innerBackGround
  exit
#
exit this menu
#
generate an overlapping grid
  solidRegion
  done choosing mappings
  change parameters
    # choose implicit or explicit interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      2 2 2 2 2 2
    exit
  # open graphics
  compute overlap
exit
#
# save an overlapping grid
save a grid (compressed)
$name
radialTravelingWaveGrid
exit

