# ----- I-Beam ----
#    Start on a flat edge since the BC is clamped. 
#    NOTE: go in a counter clockwise direction
#
#         P5                 P4    
#          +------------------+    
#          |                  |    
#          |           P2     |    
#        P6+-----+     +------+    
#             ^  |     |      P3      CH = centerHeight
#             |  |     |              CW = centerWidth 
#             CH |  C  |P1, P14    
#             |  |     |           
#             \/ |  CW |       P12 
#          +-----+     +------+    
#       EH |                  |       EH = edgeWith
#          |                  |       EW = edgeHeight
#          +------------------+    
#          <------EW -------->     
#                                  
#
$degree=3;  # degree of Nurbs 
$centerHeight=.75; $centerWidth=.25;
$edgeHeight=.25;   $edgeWidth=1.;
$cx=0.; $cy=0.;  # center for the I-beam
#
## $w=.15; $h=.75    # w=halfwidth and h=length of the arm
$nc=14;  # number of control points 
$numPerSegment=2; # number of control points per segment -- increase to make corners sharper
#
    $armWidth =.5*( $edgeWidth-$centerWidth);     # arm length f
    $x1 = $cx +.5*$centerWidth; $y1=$cy;                   # center-right point on I-beam
    $x2 = $x1;                  $y2=$y1+.5*$centerHeight;
    $x3 = $x2 +$armWidth;       $y3=$y2;
    $x4 = $x3;                  $y4=$y3+$edgeHeight;
    $x5 = $x4 -$edgeWidth;      $y5=$y4;
    $x6 = $x5;                  $y6=$y5-$edgeHeight;
    $x7 = $x6 +$armWidth;       $y7=$y6;
    $x8 = $x7;                  $y8=$y7-$centerHeight;
    $x9 = $x8 -$armWidth;       $y9=$y8;
    $x10= $x9;                  $y10=$y9-$edgeHeight;
    $x11= $x10+$edgeWidth;      $y11=$y10;
    $x12= $x11;                 $y12=$y11+$edgeHeight;
    $x13= $x12-$armWidth;       $y13=$y12;
    $x14= $x1;                  $y14=$y1;
# 
@xc=($x1,$x2,$x3,$x4,$x5,$x6,$x7,$x8,$x9,$x10,$x11,$x12,$x13,$x14);
@yc=($y1,$y2,$y3,$y4,$y5,$y6,$y7,$y8,$y9,$y10,$y11,$y12,$y13,$y14);
$radX=0; $radY=0.; # for inner background grid  **FIX ME**
$cmd="#";  $ns=0;  $arcLength=0.; $x0=$xc[0]; $y0=$yc[0];
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
   iBeamBoundaryInitial
exit
# -- interpolate the initial NURBS so that we have an arclength parameterization --
nurbs (curve)
  interpolate from a mapping
    iBeamBoundaryInitial
  mappingName
   iBeamBoundary
exit 
# 
# -- Make a hyperbolic grid --
#
  $nr = intmg( 6 + $order/2 );
  hyperbolic
    forward
    $nDist=($nr-2)*$ds;
    distance to march $nDist
    $nrm=$nr-1; 
    lines to march $nrm
    $nTheta = int($arcLength/$ds+1.5);
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
    name iBeamGridBase
  exit

  
