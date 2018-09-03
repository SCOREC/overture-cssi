# ----- Cross ----
#    Start on a flat edge since the BC is clamped. 
#    NOTE: go in a counter clockwise direction
#
#   h        11----10
#             |    | 
#             |    | 
#             |    | 
#  +w  13-----12   9------8
#       |                 |
#      1,14     0         |  w=halfwidth,  h=length of the arm
#       |                 |
#  -w   2-----3    6------7
#             |    | 
#             |    | 
#             |    | 
#             4----5 
#            -w    w
if( $w eq "" ){ $w=.15; }
if( $h eq "" ){ $h=.75; }   # w=halfwidth and h=length of the arm
$nc=14;  # number of control points 
if( $numPerSegment eq "" ){ $numPerSegment=2; } # number of control points per segment -- increase to make corners sharper
@xc=(-$h,-$h,-$w,-$w, $w, $w, $h, $h, $w, $w,-$w,-$w,-$h,-$h);
@yc=( 0.,-$w,-$w,-$h,-$h,-$w,-$w, $w, $w, $h, $h, $w, $w, 0.);
#
# Shift to specified center (cx,cy)
for( $ic=0; $ic<$nc; $ic++ ){ $xc[$ic] += $cx; $yc[$ic] += $cy; } 
