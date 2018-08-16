#
#  plotStuff plotGrid -show=diskInAChannelGride2.order2.hdf -name=diskInAChannelGride2
#  plotStuff plotGrid -show=diskInADeformingChannelGride2.order2.hdf -name=diskInADeformingChannelGride2
#
#  plotStuff plotGrid -show=rotatedEllipseGride4.order2 -name=rotatedEllipseGride4
#  plotStuff plotGrid -show=rotatedEllipseGride2.order2 -name=rotatedEllipseGride2
#
$show="diskInAChannelGride2.order2.hdf";
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name  );
#
$show
#
  # bigger
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
#
 colour boundaries by chosen name
 colour grid lines from chosen name
  $cmd="colour grid lines from chosen name\n grid colour 0 RED\n grid colour 1 BLUE\n grid colour 2 GREEN\n set view:0 -0.00279057 -5.76444e-05 0 2.04923 1 0 0 0 1 0 0 0 1";
  if( $name =~ "rotated.*" ){ $cmd ="grid colour 3 BLUE\n grid colour 4 GREEN\n grid colour 2 RED\n grid colour 0 VIOLETRED\n grid colour 1 MEDIUMVIOLETRED\n  set view:0 -0.00279057 0.00112178 0 2.62498 1 0 0 0 1 0 0 0 1"; }
  $cmd
#
  line width scale factor:0 3
  plot interpolation points 1
  # colour interpolation points 1
#
  hardcopy file name:0 $name.ps
  hardcopy save:0
  hardcopy save:0
