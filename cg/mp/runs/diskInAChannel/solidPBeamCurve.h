# ----- Solid P-Beam ----
#
#    Start on a flat edge since the BC is clamped. 
#    NOTE: go in a counter clockwise direction
#
#         P9                 P8  
#         +-------------------+                  ^
#         |                    \                 |
#         |                     \                | SL
#         |                      + P7          ^ v
#         |                      |             |
#      P1 |         C            |             | TH
#         |                      |             | 
#         |                      + P6        ^ v
#         |                     /            |
#         |      P4            /             | SL
#         |      +------------+ P5         ^ v
#         |      |                         |
#         |      |                         | SL
#         |      |                       ^ v
#         |      |                       |
#         |      |                       | BH
#         +------+                       v
#         P2     P3               
#
#         <------>                       EW = edgeWidth
#            EW                          MW = middleWidth
#                <---->                  SH = slantLength
#                  MW                    BH = bottomHeight
#                     <-->               TH = topHeight
#                      SL                 C = center = (cx,cy)
#                        <------>
#                           EW
#
if( $edgeWidth eq "" ){ $edgeWidth=1.0; }
if( $middleWidth eq "" ){ $middleWidth=1.0; }
if( $slantLength eq "" ){ $slantLength=1.0; }   
if( $bottomHeight eq "" ){ $bottomHeight=1.0; } 
if( $topHeight eq "" ){ $topHeight=1.0; } 
# center of Ibeam
if( $cx eq "" ){ $cx=0.; } 
if( $cy eq "" ){ $cy=0.; } 
#
$nc=10;  # number of control points 
$numPerSegment=3; # number of control points per segment -- increase to make corners sharper
#
$topX = $cx-.5*$middleWidth-$edgeWidth; $topY = $cy+.5*$topHeight-$slantLength;
$width = 2*$edgeWidth+$middleWidth+$slantLength;
$height = $bottomHeight+3*$slantLength+$topHeight;
$x1 = $topX;             $y1= $cy;
$x2 = $x1;               $y2=$y1-.5*$topHeight-2*$slantLength-$bottomHeight;
$x3 = $x2+$edgeWidth;    $y3=$y2;
$x4 = $x3;               $y4=$y3+$bottomHeight+$slantLength;
$x5 = $x4+$middleWidth+$edgeWidth;  $y5 =$y4;
$x6 = $x5+$slantLength; $y6 =$y5+$slantLength;
$x7 = $x6;              $y7 =$y6+$topHeight;
$x8 = $x7-$slantLength; $y8 =$y7+$slantLength;
$x9 = $topX;             $y9 =$y8;
$x10= $x1;               $y10=$y1;
# 
@xc=($x1,$x2,$x3,$x4,$x5,$x6,$x7,$x8,$x9,$x10);
@yc=($y1,$y2,$y3,$y4,$y5,$y6,$y7,$y8,$y9,$y10);
