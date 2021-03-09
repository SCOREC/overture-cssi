#
#  ----- Define a center-curve for the curved channel ------
#
$cmd="#";  $ns=100;  $arcLength=0.; $degree=3;
$xStart=0.; $xEnd=20; # Start and end positions in x
$amp = .2; $pi=4.*atan2(1.,1.);
$freq=2.;
#  
for( $i=0; $i<$ns; $i++ )\
{ $s=($i)/($ns-1);                 \
  $x=$xStart + $s*($xEnd-$xStart); \
  $y=$amp*cos(2*$pi*$freq*$s);           \
  if( $i==0 ){ $x0=$x; $y0=$y };   \
  $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 ); $x0=$x; $y0=$y; \
  $cmd .= "\n $x $y 1."; \
} 
$knots="#"; for( $i=$degree-1; $i<$ns-($degree-1); $i++ ){ $s=$i/($ns-2); $knots .= "\n $s"; }
#
nurbs (curve)
  enter control points
    $degree
    $ns
    $knots
    $cmd 
 parameterize by chord length
 #
 lines
  $lines=int($arcLength/$ds + 1.5 );
  $lines
 mappingName
   curvedChannelInitial
 # open graphics
exit
# -- interpolate the initial NURBS so that we have an arclength parameterization --
nurbs (curve)
  interpolate from a mapping
    curvedChannelInitial
 lines
  $lines=int($arcLength/$ds + 1.5 );
  $lines
  mappingName
   curvedChannelCurve
exit 
