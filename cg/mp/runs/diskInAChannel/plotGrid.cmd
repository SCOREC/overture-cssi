#
#  plotStuff plotGrid -show=diskInAChannelGride2.order2.hdf -name=diskInAChannelGride2
#  plotStuff plotGrid -show=diskInADeformingChannelGride2.order2.hdf -name=diskInADeformingChannelGride2
#
#  plotStuff plotGrid -show=rotatedEllipseGride4.order2 -name=rotatedEllipseGride4
#  plotStuff plotGrid -show=rotatedEllipseGride2.order2 -name=rotatedEllipseGride2
#
#  plotStuff plotGrid -show=solidObjectsGride4.order2 -name=solidObjectsGride4
#  plotStuff plotGrid -show=solidObjectsGride8.order2 -name=solidObjectsGride8
#
$show="diskInAChannelGride2.order2.hdf";
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name  );
#
$show
 colour boundaries by chosen name
 colour grid lines from chosen name
#
  # bigger
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
#
  $cmd="colour grid lines from chosen name\n grid colour 0 RED\n grid colour 1 BLUE\n grid colour 2 GREEN\n set view:0 -0.00279057 -5.76444e-05 0 2.04923 1 0 0 0 1 0 0 0 1";
  if( $name =~ "rotated.*" ){ $cmd ="grid colour 3 BLUE\n grid colour 4 GREEN\n grid colour 2 RED\n grid colour 0 VIOLETRED\n grid colour 1 MEDIUMVIOLETRED\n  set view:0 -0.00279057 0.00112178 0 2.62498 1 0 0 0 1 0 0 0 1"; }
  $cmd
  if( $name =~ "solid.*" ){ \
      $cmd =  "grid colour 0 BLUE\n" . \
        "grid colour 1 GREEN\n" . \
        "grid colour 2 GREEN\n" . \
        "grid colour 3 GREEN\n" . \
        "grid colour 4 GREEN\n" . \
        "grid colour 5 GREEN\n" . \
        "grid colour 6 RED\n" . \
        "grid colour 9 RED\n" . \
        "grid colour 12 RED\n" . \
        "grid colour 15 RED\n" . \
        "grid colour 18 RED\n" . \
        "grid colour 8 VIOLETRED\n" . \
        "grid colour 11 VIOLETRED\n" . \
        "grid colour 14 VIOLETRED\n" . \
        "grid colour 17 VIOLETRED\n" . \
        "grid colour 20 VIOLETRED\n" . \
        "grid colour 7 ORANGERED\n" . \
        "grid colour 10 ORANGERED\n" . \
        "grid colour 13 ORANGERED\n" . \
        "grid colour 16 ORANGERED\n" . \
        "grid colour 19 ORANGERED\n" . \
        "set view:0 0.109844 0.0282565 0 2.00447 1 0 0 0 1 0 0 0 1"; }else{ $cmd="#"; }
$cmd
#
  line width scale factor:0 3
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  plot interpolation points 1
  if($name =~ "solid.*" ) { $cmd="line width scale factor:0 2\n  plot interpolation points 0"; }else{ $cmd="#"; }
  $cmd
  # colour interpolation points 1
#
  hardcopy file name:0 $name.ps
  hardcopy save:0
pause
#
  set view:0 -0.0016676 -0.0115571 0 5.74392 1 0 0 0 1 0 0 0 1
  $plotName = $name . "Zoom.ps";
  hardcopy file name:0 $plotName
  hardcopy save:0
