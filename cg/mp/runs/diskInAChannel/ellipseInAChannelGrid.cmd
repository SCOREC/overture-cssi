#
# Grid for an elastic ellipse in a channel
#
#
# usage: ogen [noplot] ellipseInAChannelGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] 
# Options:
#    -ra,rb
# 
# examples:
#   ogen -noplot ellipseInAChannelGrid -interp=e -factor=2
#   ogen -noplot ellipseInAChannelGrid -interp=e -factor=4
#
# -- rotated 
# ogen -noplot ellipseInAChannelGrid -radX=.6 -radY=1.2 -corex=0 -corey=.5 -ra=.2 -theta=10. -prefix=rotatedEllipseGrid -interp=e -factor=4
#
#
$prefix="ellipseInAChannelGrid";
$order=2; $factor=1; $interp="e"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids"; $ml=0; 
$name=""; 
# Fixed radial distance
$nDistInterface=.1;
#
$xa=-3; $xb=3; $ya=-2; $yb=2; 
$radX=1.5; $radY=.6; 
$ra=.2; $rb=$ra+$nDistInterface; $corex=.5; $corey=0.;   # ellipse has a hollow core
$numberOfVolumeSmooths=0; 
$theta=0.;
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,\
            "radX=f"=> \$radX,"radY=f"=> \$radY,\
            "ra=f"=> \$ra,"rb=f"=> \$rb,"corex=f"=> \$corex,"corey=f"=> \$corey,\
            "interp=s"=> \$interp,"name=s"=> \$name,"per=i"=>\$per,\
            "xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "prefix=s"=>\$prefix,"theta=f"=>\$theta );
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
$theta=$theta*$pi/180.;
# 
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }
$minRad=min($radX,$radY);
# $nDistInterface=.1*$minRad;
$rb=$ra+$nDistInterface;
#
$bcInterface=100;  # bc for interfaces
$shareInterface=100;        # share value for interfaces
#
create mappings
#
# --- make a nurbs for a circle ----
#   **** Later we could make different curves for the interface ****
$cx=0.; $cy=0.; 
# include ellipseCurve.h
include rotatedEllipseCurve.h
# 
# -- Make a hyperbolic grid --
#
  $nr = intmg( 5 + $order/2 );
  hyperbolic
    forward
    # *wdh* May 14, 2018 $nrm = intmg($nDistInterface/$ds+4.5);
    # Make sure the number of grid cells doubles when ds is halved.
    $nDistExtraFactor=1.2; # make grid a bit finer in the normal direction
    $nrm = intmg(($nDistInterface*$nDistExtraFactor)/$ds+.5);
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
    # use fourth order interpolant to define the mapping: *wdh* May 8, 2018
    fourth order
    # evaluate as nurbs 1
    boundary conditions
      -1 -1 100 0 0 0
    share 
       0 0 100 0 0 0
    name fluidInterface
  exit
#
# -- Solid ellipse near interface ----
#
  hyperbolic
    backward
    # Fixed radial distance
    # *wdh* May 14, 2018 $nrm = intmg($nDistInterface/$ds+4.5);
    # Make sure the number of grid cells doubles when ds is halved.
    $nDistExtraFactor=1.2; # make grid a bit finer in the normal direction
    $nrm = intmg(($nDistInterface*$nDistExtraFactor)/$ds+.5);
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
    # use fourth order interpolant to define the mapping: *wdh* May 8, 2018
    fourth order
    # evaluate as nurbs 1
    boundary conditions
      -1 -1 100 0 0 0
    share 
       0 0 100 0 0 0
    name solidInterface
  exit
#
#   Solid annular core 
  Annulus
  center: $corex $corey
  inner and outer radii
    $innerRad=$ra; $outerRad=$rb; 
    $innerRad $outerRad
  lines
    $nTheta = intmg( 2.*$pi*($innerRad+$outerRad)*.5/$ds + 3.5 );
    $nr = intmg( $nDistExtraFactor*($outerRad-$innerRad)/$ds + 1.5 );
    $nTheta $nr
  boundary conditions
    -1 -1 5 0 
  share
     0  0 5 0 
  mappingName
    solidCore
exit
#
#  --- solid backGround Grid ---
#
 rectangle
  set corners
    $rmax=max($radX,$radY);
    $xas=-$rmax; $xbs=-$xas; $yas=-$rmax; $ybs=-$yas; 
    $xas $xbs $yas $ybs
    $refineFactor=1.;
  lines
    $nx = int( ($xbs-$xas)/$ds*$refineFactor + .5 ); 
    $ny = int( ($ybs-$yas)/$ds*$refineFactor + .5 ); 
    $nx $ny
  boundary conditions
    0 0 0 0 
  share 
    0 0 0 0 
  mappingName
    solidBackGround
  exit
#
#  --- fluid Channel Grid ---
#
 rectangle
  set corners
    $xa $xb $ya $yb
  lines
    $nx = int( ($xb-$xa)/$ds + .5 ); 
    $ny = int( ($yb-$ya)/$ds + .5 ); 
    $nx $ny
  boundary conditions
    1 2 3 4 
  share 
    0 0 0 0 
  mappingName
    fluidChannel
  exit
#
exit this menu
#
generate an overlapping grid
  solidBackGround
  solidCore
  solidInterface
  fluidChannel
  fluidInterface
  done choosing mappings
  change parameters
    specify a domain
     solidDomain
       solidBackGround
       solidCore
       solidInterface
    done
    specify a domain
     fluidDomain
       fluidChannel
       fluidInterface
    done
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
ellipseInAChannelGrid
exit
