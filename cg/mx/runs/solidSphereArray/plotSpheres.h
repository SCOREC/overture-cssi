#
# Set up 3D contour plot 
# 
if( $yContour eq "" ){ $yContour=-.65; }
if( $zContour eq "" ){ $zContour=-.65; }
#    $ypos=-.65; $zpos=-.65;  # y and z-positions of contour plane 
plot:Ey
contour
  plot contour lines (toggle)
  # if( $grid =~ /3d/ || $grid =~ /Sphere/ || $grid =~ /Ellipsoid/ )
  $cmd="delete contour plane 2\n delete contour plane 1\n delete contour plane 0"; 
  $cmd
  add contour plane  0  0  1  0     0   $zContour
  add contour plane  0  1  0  0 $yContour  0 
exit
grid
 plot block boundaries 0
  plot shaded surfaces (3D) 0
  toggle grid 0 0
  coarsening factor 2
exit this menu
