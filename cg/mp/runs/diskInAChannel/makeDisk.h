#
#   Create grids for a deforming disk
# Input:
#      $count, $cx,$cy, $ra,$b
#
$count = $count + 1;
$curveName="curve$count"; $fluidName="fluidInterface$count"; $share=99+$count; $solidName="solidDisk$count";
#
# --- make a nurbs for a circle ----
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
      -1 -1 $share 0 0 0
    share 
       0 0 $share 0 0 0
    name $fluidName
  exit
#
# -- Solid annulus ----
#
Annulus
  inner and outer radii
    $innerRad=$ra; $outerRad=$rb; 
    # $innerRad $outerRad
    $outerRad $innerRad
  center: $cx $cy
  lines
    # $nTheta = intmg( 2.*$pi*($innerRad+$outerRad)*.5/$ds + 3.5 );
    $nTheta = $nThetaInterface;
    $nr = intmg( $nDistExtraFactor*($outerRad-$innerRad)/$ds + 1.5 );
    $nTheta $nr
  boundary conditions
#    -1 -1 5 $share 
    -1 -1 $share 5
  share
#     0  0 5  $share
     0  0 $share 5  
  mappingName
    $solidName
exit
