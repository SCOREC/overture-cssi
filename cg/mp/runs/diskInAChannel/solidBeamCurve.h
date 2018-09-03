# ----- Solid Beam ----
#
#    Start on a flat edge since the BC is clamped. 
#    NOTE: go in a counter clockwise direction
#
#         P4                     
#          +------------------------+  P3  
#          |                        |    
#          |                        |    
#      EH  |           C            |     C = center = (cx,cy)
#          |                 X      |     X = core-center (xCore,yCore), offset from center
#          |           P6           |     EH = edge height
#          +-----------+------------+     Ew = edge width 
#         P5           P1           P2
#           <----------EW ---------->     
#                                  
#
if( $edgeHeight eq "" ){ $edgeHeight=.25; }   
if( $edgeWidth eq "" ){ $edgeWidth=1.; } 
# center 
if( $cx eq "" ){ $cx=0.; } 
if( $cy eq "" ){ $cy=0.; } 
#
$nc=6;  # number of control points 
$numPerSegment=3; # number of control points per segment -- increase to make corners sharper
#
    $x1 = $cx;                  $y1=$cy-.5*$edgeHeight;
    $x2 = $cx+.5*$edgeWidth;    $y2=$y1;
    $x3 = $x2;                  $y3=$cy+.5*$edgeHeight;
    $x4 = $cx-.5*$edgeWidth;    $y4=$y3;
    $x5 = $x4;                  $y5=$y1;
    $x6 = $x1;                  $y6=$y1;
# 
@xc=($x1,$x2,$x3,$x4,$x5,$x6);
@yc=($y1,$y2,$y3,$y4,$y5,$y6);
