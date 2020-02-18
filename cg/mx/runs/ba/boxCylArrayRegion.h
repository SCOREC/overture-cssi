#
# Circular region with co-centrix arrays of squares and circles
# 
forcing options...
define material region...
material file: $matFile2
  $xc=.0; $yc=0; $zc=0; $radius=.7;
  cylinder: $xc $yc $zc $radius 0 1 (x0,y0,z0, radius, za,zb)
# 
material file: $matFile3
  $cmd ="#";
  $pi = 4*atan2(1.,1.);
  # outer ring of squares 
  $xc0=0; $yc0=0; $w=.1; $h=.1;  $rad=.55; $numBox=20;
  for( $i=0; $i<$numBox; $i++ ){ $theta=2*$pi*($i)/$numBox; $xc=$xc0+$rad*cos($theta); $yc=$yc0+$rad*sin($theta); \
    $xa=$xc-$w*.5; $xb=$xc+$w*.5; $ya=$yc-$h*.5; $yb=$yc+$h*.5; \
    $cmd .= "\n box: $xa $xb $ya $yb 0 1 (xa,xb,ya,yb,za,zb)";  } 
 $cmd
 # Inner ring of cylinders
  $xc0=0; $yc0=0; $r=.06; $rad=.3; $numCyl=10;
  for( $i=0; $i<$numCyl; $i++ ){ $theta=2*$pi*($i)/$numCyl; $xc=$xc0+$rad*cos($theta); $yc=$yc0+$rad*sin($theta); \
    $cmd .= "\n cylinder: $xc $yc $zc $r 0 1 (x0,y0,z0, radius, za,zb)";  } 
 $cmd
#   $xc0=0.2; $yc0=0; $w=.1; $h=.1;  $rad=.3; $numBox=10;
#   for( $i=0; $i<$numBox; $i++ ){ $theta=2*$pi*($i)/$numBox; $xc=$xc0+$rad*cos($theta); $yc=$yc0+$rad*sin($theta); \
#     $xa=$xc-$w*.5; $xb=$xc+$w*.5; $ya=$yc-$h*.5; $yb=$yc+$h*.5; \
#     $cmd .= "\n box: $xa $xb $ya $yb 0 1 (xa,xb,ya,yb,za,zb)";  } 
#  $cmd
#  $xa=$xc-$w*.5; $xb=$xc+$w*.5; $ya=$yc-$h*.5; $yb=$yc+$h*.5;
#  box: $xa $xb $ya $yb 0 1 (xa,xb,ya,yb,za,zb)
continue
