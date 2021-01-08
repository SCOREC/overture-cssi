#
#   --- plot grid for an array of solid ellipsoid ---
#
#     plotStuff plotEllipsoidArray.cmd -show=fourSolidEllipsoidsGride3.order2.hdf -name=fourSolidEllipsoidsGrid.ps
#     plotStuff plotEllipsoidArray.cmd -show=eightSolidEllipsoidsGride3.order2.hdf -name=eightSolidEllipsoidsGride3.ps
# 
#     plotStuff plotEllipsoidArray.cmd -show=36SolidEllipsoidsGride4.order2.hdf -name=36SolidEllipsoidsGrid.ps -ns=36 
#     plotStuff plotEllipsoidArray.cmd -show=125SolidEllipsoidsGride4.order2.hdf -name=125SolidEllipsoidsGrid.ps -ns=125 
#
$show="fourSolidEllipsoidsGride3.order2.hdf"; $name="fourSolidEllipsoidsGride3.ps"; $ns=8; 
#
GetOptions( "show=s"=>\$show,"name=s"=>\$name,"ns=i"=>\$ns );
#
$show
  DISPLAY AXES:0 0
  DISPLAY SQUARES:0 0
  plot block boundaries 0
# 
  toggle grid 0 0
  # Note : we only colour outer domain grids -- 3 per ellipsoid 
  $j=1; $i=1; 
  for( $k=0; $k<$ns; $k++ ){ \
     $cmd .= "grid colour $j COPPER\n";  $j++; \
     $cmd .= "grid colour $j COPPER\n";  $j++; \
     $cmd .= "grid colour $j COPPER\n";  $j++; $i++; if( $i > $ns ){ last; }\
     $cmd .= "grid colour $j CADETBLUE\n";  $j++; \
     $cmd .= "grid colour $j CADETBLUE\n";  $j++; \
     $cmd .= "grid colour $j CADETBLUE\n";  $j++;  $i++; if( $i > $ns ){ last; }\
     $cmd .= "grid colour $j JADE\n";  $j++; \
     $cmd .= "grid colour $j JADE\n";  $j++; \
     $cmd .= "grid colour $j JADE\n";  $j++;  $i++; if( $i > $ns ){ last; }\
     $cmd .= "grid colour $j LIGHTSTEELBLUE\n";  $j++; \
     $cmd .= "grid colour $j LIGHTSTEELBLUE\n";  $j++; \
     $cmd .= "grid colour $j LIGHTSTEELBLUE\n";  $j++;  $i++; if( $i > $ns ){ last; }\
     $cmd .= "grid colour $j BRONZE\n"; $j++; \
     $cmd .= "grid colour $j BRONZE\n"; $j++; \
     $cmd .= "grid colour $j BRONZE\n"; $j++;  $i++; if( $i > $ns ){ last; }\
     $cmd .= "grid colour $j MEDIUMSEAGREEN\n";  $j++; \
     $cmd .= "grid colour $j MEDIUMSEAGREEN\n";  $j++; \
     $cmd .= "grid colour $j MEDIUMSEAGREEN\n";  $j++;  $i++; if( $i > $ns ){ last; }\
     $cmd .= "grid colour $j ORANGE\n"; $j++; \
     $cmd .= "grid colour $j ORANGE\n"; $j++; \
     $cmd .= "grid colour $j ORANGE\n"; $j++;  $i++; if( $i > $ns ){ last; }\
     $cmd .= "grid colour $j BRASS\n";  $j++; \
     $cmd .= "grid colour $j BRASS\n";  $j++; \
     $cmd .= "grid colour $j BRASS\n";  $j++; $i++; if( $i > $ns ){ last; }\
     $cmd .= "grid colour $j SKYBLUE\n";  $j++; \
     $cmd .= "grid colour $j SKYBLUE\n";  $j++; \
     $cmd .= "grid colour $j SKYBLUE\n";  $j++; $i++; if( $i > $ns ){ last; }\
     $cmd .= "grid colour $j MEDIUMTURQUOISE\n";  $j++; \
     $cmd .= "grid colour $j MEDIUMTURQUOISE\n";  $j++; \
     $cmd .= "grid colour $j MEDIUMTURQUOISE\n";  $j++; $i++; if( $i > $ns ){ last; }\
     $cmd .= "grid colour $j MEDIUMAQUAMARINE\n";  $j++; \
     $cmd .= "grid colour $j MEDIUMAQUAMARINE\n";  $j++; \
     $cmd .= "grid colour $j MEDIUMAQUAMARINE\n";  $j++; $i++; if( $i > $ns ){ last; }\
  }
  $cmd.="#";
echo to terminal 0
  $cmd
echo to terminal 1
  coarsening factor 2
  set view:0 0.014107 -0.00202235 0 1.27542 0.766044 0.219846 -0.604023 0 0.939693 0.34202 0.642788 -0.262003 0.719846
  # open graphics 
#
  line width scale factor:0 2
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  hardcopy file name:0 $name
  plot
  hardcopy save:0

  grid colour 1 GREEN
  grid colour 2 GREEN
  grid colour 3 GREEN
  grid colour 4 CADETBLUE
  grid colour 5 CADETBLUE
  grid colour 6 CADETBLUE
  grid colour 7 JADE
  grid colour 8 JADE
  grid colour 9 JADE
  grid colour 10 LIGHTSTEELBLUE
  grid colour 11 LIGHTSTEELBLUE
  grid colour 12 LIGHTSTEELBLUE
  grid colour 13 BRONZE
  grid colour 14 BRONZE
  grid colour 15 BRONZE
  grid colour 16 MEDIUMSEAGREEN
  grid colour 17 MEDIUMSEAGREEN
  grid colour 18 MEDIUMSEAGREEN
  grid colour 19 ORANGE
  grid colour 20 ORANGE
  grid colour 21 ORANGE
  grid colour 22 BRASS
  grid colour 23 BRASS
  grid colour 24 BRASS
#
  coarsening factor 2
  set view:0 0.014107 -0.00202235 0 1.27542 0.766044 0.219846 -0.604023 0 0.939693 0.34202 0.642788 -0.262003 0.719846
pause 
#
  line width scale factor:0 2
  hardcopy vertical resolution:0 2048
  hardcopy horizontal resolution:0 2048
  hardcopy file name:0 $name
  plot
  hardcopy save:0
