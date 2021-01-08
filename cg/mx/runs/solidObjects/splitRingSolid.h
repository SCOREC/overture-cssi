#
#  ----- Build a split ring ------
#
#       13                            12
#         +---------------------------+ (H,W)
#         |                           |
#         |            w          9   |
#         |   +------------------+    |---
#         |   |8                 |    |  d 
#      14 |   |               10 +----+---
#        1| w |                        11
#         |   |                        
#         |   |               5  +----+4
#         |   |                  |    |
#         |   +------------------+    |
#         |    7       w        6     |
#         |                           |
#       2 +---------------------------+ (H,-W)
#      (-W,H)                         3
#
#
#
# Optional input:
#   $numPointsPerSegment = number of extra points per segment (default=3) increase to make sharper corners
#   $arcLengthScaleFactor : scale arcLenght by this value for computing number of grid point default =1
#   $numRadial :  number of point in radial direction is $numRadial + $order/2  (default=6)
#
# number of control points per segment -- increase to make corners sharper: 
if( $numPointsPerSegment eq "" ){ $numPerSegment=3; }else{ $numPerSegment=$numPointsPerSegment;}   # 
if( $arcLengthScaleFactor eq "" ){ $arcLengthScaleFactor=1; };
if( $numRadial eq "" ){ $nr0=6; }else{ $nr0=$numRadial;}   # 
#
$degree=3;  # degree of Nurbs 
$H=.5; $W=.5; $w=.2; $d=.2;
# $cx=0.; $cy=0.;  # center for the I-beam
$numberOfVolumeSmooths=20;
#
$nc=14;  # number of control points 
# $numPerSegment=3; # number of control points per segment -- increase to make corners sharper
#
    $x1 =-$W;    $y1 =0;
    $x2 = $x1;   $y2 =-$H;
    $x3 = $W;    $y3 =$y2;
    $x4 = $x3;   $y4 =$y2+$w+$d;
    $x5 = $W-$w; $y5 =$y4;
    $x6 = $x5;   $y6 =$y5-$d;
    $x7 =-$W+$w; $y7 =$y6;
    $x8 = $x7;   $y8 =$H-$w;
    $x9 =$W-$w;  $y9 =$y8;
    $x10=$x9;    $y10=$H-$w-$d;
    $x11= $W;    $y11=$y10;
    $x12=$x11;   $y12=$H;
    $x13=-$W;    $y13=$y12;
    $x14=$x1;    $y14=$y1;
# 
@xc=($x1,$x2,$x3,$x4,$x5,$x6,$x7,$x8,$x9,$x10,$x11,$x12,$x13,$x14);
@yc=($y1,$y2,$y3,$y4,$y5,$y6,$y7,$y8,$y9,$y10,$y11,$y12,$y13,$y14);
$radX=0; $radY=0.; # for inner background grid  **FIX ME**
$x0=$xc[0]; $y0=$yc[0];
$cmd="#";  $ns=0;  $arcLength=0.; 
for( $ic=0; $ic<$nc-1; $ic++ ){ $numpt=$numPerSegment; if( $ic == $nc-2 ){ $numpt=$numPerSegment+1;} \
for( $i=0; $i<$numpt; $i++ ){ $s=($i)/($numPerSegment); $ns=$ns+1;  \
   $x=(1.-$s)*$xc[$ic]+$s*$xc[$ic+1]; \
   $y=(1.-$s)*$yc[$ic]+$s*$yc[$ic+1];   \
   $radX = max($radX,abs($x)); $radY = max($radY,abs($y));  \
   $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 ); $x0=$x; $y0=$y; \
   $cmd .= "\n $x $y 1."; } }\
  $knots="#"; for( $i=$degree-1; $i<$ns-($degree-1); $i++ ){ $s=$i/($ns-2); $knots .= "\n $s"; } \
#
nurbs (curve)
  periodicity
   2
  enter control points
    $degree
    $ns
    $knots
    $cmd 
 parameterize by chord length
 #
 lines
  $lines=intmg($arcLength/$ds + 1.5 );
  $lines
 mappingName
   splitRingBoundaryInitial
exit
# -- interpolate the initial NURBS so that we have an arclength parameterization --
nurbs (curve)
  interpolate from a mapping
    splitRingBoundaryInitial
 lines
  $lines=intmg($arcLength/$ds + 1.5 );
  $lines
  mappingName
   splitRingBoundary
exit 
# 
# -- Make a hyperbolic grid (OUTSIDE) --
#
  $nr = intmg( $nr0 + $order/2 );
  hyperbolic
    forward
    $nDist=($nr-2.5)*$ds;
    distance to march $nDist
    $nrm=$nr-1; 
    lines to march $nrm
    $nTheta = int($arcLengthScaleFactor*$arcLength/$ds+1.5);
    points on initial curve $nTheta
    uniform dissipation 0.075
    volume smooths $numberOfVolumeSmooths
    equidistribution 0 (in [0,1])
    #
    spacing: geometric
    geometric stretch factor 1.05 
    #
    generate
    # open graphics
    boundary conditions
      -1 -1 7 0 0 0
    share 
       0 0 100 0 0 0
    name splitRingGridBase
  exit
# -- Make a hyperbolic grid (INSIDE) --
#
  $nr = intmg( $nr0-1 + $order/2 );
  hyperbolic
    backward
    $nDist=($nr-3)*$ds;
    distance to march $nDist
    $nrm=$nr-1; 
    lines to march $nrm
    # $nTheta = int($arcLengthScaleFactor*$arcLength/$ds+1.5);
    points on initial curve $nTheta
    uniform dissipation 0.05
    volume smooths $numberOfVolumeSmooths
    equidistribution 0 (in [0,1])
    #
    # spacing: geometric
    # geometric stretch factor 1.05 
    #
    generate
    # open graphics
    boundary conditions
      -1 -1 7 0 0 0
    share 
       0 0 100 0 0 0
    name splitRingGridBaseSolid
  exit
