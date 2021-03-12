#
#  Curve for circle with noise
#    NOTE: go in a counter clockwise direction
#
$cmd="#";
$ra=.75;  # inner radius
$ns=10*$factor+1;  # number of spline points
$arcLength=0.;
$freq1=15; $amp1=.02;
# $freq2=71; $amp2=.007;
$freq2=61; $amp2=.007;
# $freq2=91; $amp2=.007;
$radX=$ra + $amp1 + $amp2; $radY=$radX;
for( $i=0; $i<$ns; $i++ ){ $s=2.*$pi*($i-1.)/($ns-1.);  \
  $r=$ra + $amp1*cos($freq1*($s+.5)) + $amp2*cos($freq2*($s+.3));								\
  $xx=$r*cos($s); \
  $yy=$r*sin($s);				\
 if( $i > 0 ){ $arcLength=$arcLength + sqrt( ($xx-$x0)**2 + ($yy-$y0)**2 );} $x0=$xx; $y0=$yy; \
$cmd .= "\n $xx $yy"; }
#
spline
  #
  enter spline points
    $ns
    $cmd
  lines
    $ns
    periodicity
      2
  mappingName
    curveBoundary
 exit
# for( $i=0; $i<$ns; $i++ ){ $s=2.*$pi*($i-1.)/($ns-1.);  $r=$ra; \
#  $xx=$r*cos($s) + $amp1*cos($freq1*($s+.5)) + $amp2*cos($freq2*($s+.3)); \
#  $yy=$r*sin($s) + $amp1*sin($freq1*($s+.5)) + $amp2*sin($freq2*($s+.3));				\
#  if( $i > 0 ){ $arcLength=$arcLength + sqrt( ($xx-$x0)**2 + ($yy-$y0)**2 );} $x0=$xx; $y0=$yy; \
