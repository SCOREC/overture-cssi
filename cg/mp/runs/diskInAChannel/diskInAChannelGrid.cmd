#
# Grid for an elastic disk in a channel
#
#
# usage: ogen [noplot] diskInAChannelGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] 
# Options:
#    -ra,rb
# 
# examples:
#   ogen -noplot diskInAChannelGrid -interp=e -factor=2
#   ogen -noplot diskInAChannelGrid -interp=e -factor=4
#
#
$prefix="diskInAChannelGrid";
$order=2; $factor=1; $interp="e"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids"; $ml=0; 
$name=""; 
$xa=-3; $xb=3; $ya=-1.75; $yb=1.75; 
$ra=.4; $rb=.8;   # disk is hollow 
$numberOfVolumeSmooths=0; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"ra=f"=> \$ra,"rb=f"=> \$rb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"per=i"=>\$per,\
            "xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "prefix=s"=>\$prefix );
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
$radX=$rb; $radY=$rb;
include ellipseCurve.h
# 
# -- Make a hyperbolic grid --
#
  $nr = intmg( 5 + $order/2 );
  hyperbolic
    forward
    # Fixed radial distance
    $nDistInterface=.2;
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
# -- Solid annulus ----
#
Annulus
  inner and outer radii
    $innerRad=$ra; $outerRad=$rb; 
    # $innerRad $outerRad
    $outerRad $innerRad
  lines
    # $nTheta = intmg( 2.*$pi*($innerRad+$outerRad)*.5/$ds + 3.5 );
    $nTheta = $nThetaInterface;
    $nr = intmg( $nDistExtraFactor*($outerRad-$innerRad)/$ds + 1.5 );
    $nTheta $nr
  boundary conditions
#    -1 -1 5 100 
    -1 -1 100 5
  share
#     0  0 5  100 
     0  0 100 5  
  mappingName
    solidDisk
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
  solidDisk
  fluidChannel
  fluidInterface
  done choosing mappings
  change parameters
    specify a domain
     solidDomain
       solidDisk
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
diskInAChannelGrid
exit

