#
#    plotStuff plotSphereArray -show=ssa27G2O2.show -name=ssa27
# 
#    plotStuff plotSphereArray -show=ssa64G2.show -name=ssa64
#
$show="spiralWire.show"; $name=""; 
# ----------------------------- get command line arguments ---------------------------------------
GetOptions( "show=s"=>\$show,"name=s"=>\$name );
# -------------------------------------------------------------------------------------------------
$show
previous
# solution: 6
derived types
  E field norm
exit
plot:eFieldNorm
# 
grid
  toggle grid 0 0
  plot block boundaries 0
  plot shaded surfaces (3D) 0
  coarsening factor 2
#  colour grid lines by grid number
  colour grid lines black
  set view:0 0 0 0 1 0.896724 0.344593 -0.277744 -0.130018 0.804953 0.578917 0.423062 -0.483017 0.766624
  exit this menu
contour
  contour lines 0
  delete contour plane 0
  delete contour plane 0
  add contour plane  0.00000e+00  1.00000e+00  0.00000e+00 -2.55110e+00 -2.46824e+00  3.49787e+00 
  min max 0 1 
exit
pause
DISPLAY AXES:0 0
DISPLAY SQUARES:0 0
DISPLAY COLOUR BAR:0 0
DISPLAY LABELS:0 0
set view:0 0.0847458 -0.00219197 0 1.0153 0.704945 0.282896 -0.650402 -0.00934022 0.920635 0.390313 0.7092 -0.269074 0.651639
movie file name: twentySevenSphereScat
save movie files 1


show movie





















  toggle grid 0 0
  plot grid lines 0
  plot block boundaries 0
  grid colour 1 BRASS
  grid colour 2 BRASS
  grid colour 3 BRASS
exit this menu
# 
derived types
E field norm
exit
plot:eFieldNorm
contour
  contour lines 0
exit
