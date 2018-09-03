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
#         |   |        C                         C=center =(cx,cy)
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
#
if( $H ne "" ){ $H=.6;}
if( $W ne "" ){ $W=.5; }
if( $w ne "" ){ $w=.25;}
if( $d ne "" ){ $d=.2; }
# center 
if( $cx eq "" ){ $cx=0.; } 
if( $cy eq "" ){ $cy=0.; } 
#
$nc=14;  # number of control points 
$numPerSegment=3; # number of control points per segment -- increase to make corners sharper
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
# Shift to specified center (cx,cy)
for( $ic=0; $ic<$nc; $ic++ ){ $xc[$ic] += $cx; $yc[$ic] += $cy; } 
