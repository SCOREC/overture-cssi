# ----- Solid I-Beam ----
#
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
#             CH |  C  |P1, P14       C = center = (cx,cy)
#             |  |     |           
#             \/ |  CW |       P12 
#          +-----+     +------+    
#       EH |        X         |       EH = edgeWith
#          |                  |       EW = edgeHeight
#          +------------------+       X = core-center (xCore,yCore), offset from center
#          <------EW -------->     
#                                  
#
if( $centerHeight eq "" ){ $centerHeight=.75; }
if( $centerWidth eq "" ){ $centerWidth=.25; }
if( $edgeHeight eq "" ){ $edgeHeight=.25; }   
if( $edgeWidth eq "" ){ $edgeWidth=1.; } 
# center of Ibeam
if( $cx eq "" ){ $cx=0.; } 
if( $cy eq "" ){ $cy=0.; } 
#
$nc=14;  # number of control points 
$numPerSegment=3; # number of control points per segment -- increase to make corners sharper
#
    $armWidth =.5*( $edgeWidth-$centerWidth);     # arm length 
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
