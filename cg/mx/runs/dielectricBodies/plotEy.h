if( $az==0 ){ $cmd="plot:Ey"; }else{ $cmd="plot:Ez"; }
$cmd 
contour
  plot contour lines (toggle)
  # vertical scale factor 0.2
  # min max -1.1 1.1
  # plot a contour plane in 3d 
  if( $grid =~ /3d/ || $grid =~ /Sphere/ || $grid =~ /Ellipsoid/ ){ $cmd="delete contour plane 2\n delete contour plane 1\n delete contour plane 0\n add contour plane  0.00000e+00  0.00000e+00  1.00000e+00 0 0 0"; }else{ $cmd="#"; }
  $cmd
exit
  
