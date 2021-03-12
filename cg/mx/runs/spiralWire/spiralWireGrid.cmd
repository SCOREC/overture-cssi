#
#  Ogen: Grid for a  SPIRAL WIRE 
#
# usage: ogen [-noplot] spiralWireGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] -nrExtra=<> -ml=<> -ns=<>
#
#  nrExtra: extra lines to add in the radial direction on the sphere grids 
#  ml = number of (extra) multigrid levels to support
#  ns : number of spheres, 1 or 2 
#  xa, xb, ya, yb, za, zb : bounds on the channel
# 
# examples:
#     ogen -noplot spiralWireGrid -factor=1 -order=2
#     ogen -noplot spiralWireGrid -factor=1 -interp=e -order=2
#     ogen -noplot spiralWireGrid -factor=1 -order=4
# 
#     ogen -noplot spiralWireGrid -factor=2 -interp=e -order=4
#
$prefix="spiralWireGrid"; 
# Background: 
$xa=-2; $xb=2; $ya=-1.75; $yb=1.75; $za=-1.75; $zb=1.75; 
#
$nrExtra=3; $loadBalance=0; $ml=0;
$numGhost=-1; # set to over-ride default
$numGhostNurbs=2; # For conversion to Nurbs
#
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids"; $dse=0.; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=i"=> \$factor,"nrExtra=i"=> \$nrExtra,"interp=s"=> \$interp,\
            "loadBalance=i"=>\$loadBalance,"xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "za=f"=>\$za,"zb=f"=>\$zb,"ml=i"=>\$ml,"numGhost=i"=>\$numGhost );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; $numGhostNurbs=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; $numGhostNurbs=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; $dse=1.; }
# 
$suffix = ".order$order";
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $ns eq 1 ){ $prefix = "oneSphereInAChannel"; }
$name = $prefix . "$interp$factor" . $suffix . ".hdf";
# 
# $ds=.1/$factor;
# ********************* NOTE : FIX ME NOT STEUP FOR FACTOR YET ***
$ds=.025;
$pi=4.*atan2(1.,1.);
# 
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
#
# ------ create mappings -------
#
create mappings
  read iges file
    SpiralWire.igs
    continue
    choose all
# 
    CSUP:refine plot
    refine surface 5
    refine surface 0
    CSUP:determine topology
    deltaS 0.01
    maximum area .001
    merge tolerance 0.01
    compute topology
    exit
    CSUP:mappingName SpiralWire.igs
    exit
#
#
#
  builder...
    target grid spacing .025 .025 (tang,norm)((<0 : use default)
#
#   Inner Cap
#  
    create surface grid...
      surface grid options...
      initial curve:points on surface
      choose point on surface 7 7.633909e-01 1.583710e-01 -1.506164e-03 6.177111e-01 3.411315e-01
      choose point on surface 7 7.694557e-01 1.231155e-01 -8.488371e-04 1.453908e-01 1.605104e-02
      choose point on surface 7 7.798769e-01 9.525114e-02 -1.257311e-03 4.885886e-01 3.778090e-02
      choose point on surface 7 7.955344e-01 7.103412e-02 -1.001720e-03 1.039155e-01 8.721846e-01
      choose point on surface 7 8.220276e-01 4.664295e-02 -1.044733e-03 4.441570e-01 4.283682e-02
      choose point on surface 7 8.572396e-01 2.821526e-02 -3.508096e-04 1.704907e-02 3.711218e-01
      choose point on surface 4 8.842202e-01 1.982088e-02 4.402358e-04 6.494436e-01 2.541835e-02
      choose point on surface 4 9.162447e-01 1.336756e-02 1.899146e-04 5.340318e-01 4.432188e-01
      choose point on surface 4 1.040082e+00 -6.361089e-03 1.404762e-03 2.136144e-01 7.252616e-01
      choose point on surface 4 1.087355e+00 -1.078098e-02 1.765331e-03 9.278660e-02 1.429520e-01
      choose point on surface 4 1.148179e+00 -8.591129e-04 6.322735e-04 3.018321e-02 8.246916e-01
      choose point on surface 4 1.192766e+00 2.347606e-02 1.542410e-03 5.188504e-02 1.371897e-03
      choose point on surface 4 1.229451e+00 6.216334e-02 1.970997e-03 6.532666e-01 2.723166e-01
      choose point on surface 4 1.248927e+00 9.930017e-02 1.769312e-03 9.390988e-01 6.090116e-02
      choose point on surface 4 1.258054e+00 1.317437e-01 1.680958e-03 3.752150e-02 9.182196e-01
      done
      forward and backward
      BC: left (forward) outward splay
      BC: right (forward) outward splay
      BC: left (backward) outward splay
      BC: right (backward) outward splay
      points on initial curve 26
      $lines = 12; #  + $order-2; 
      lines to march $lines, $lines (forward,backward)  
      # lines to march 12, 12 (forward,backward)  
      outward splay .5 .5 (left, right for outward splay BC)
      normal blending 6 6 (lines: left, right)
      generate
      smoothing...
        GSM:BC: top smoothed
        GSM:BC: bottom smoothed
        GSM:BC: right smoothed
        GSM:BC: left smoothed
        GSM:smooth grid
      name wireCapInnerSurface
      exit
#
#   Wire surface
#
    create surface grid...
      choose edge curve 12 1.015398e+00 2.006594e-01 2.500601e-01 
      choose edge curve 3 1.015045e+00 2.007773e-01 -2.499750e-01 
      done
      marching options...
      uniform dissipation 100
      equidistribution 0.5(in [0,1])
      forward and backward
      # lines to march 667, 3 (forward,backward)  
      lines to march 668, 4 (forward,backward)  
      points on initial curve 63
      name spiralWireSurface
      generate
      smoothing...
      GSM:number of iterations 10
      GSM:smooth grid
      exit
#
#   Wire cap on outer spiral tip
#
    create surface grid...
      initial curve:points on surface
      choose point on surface 6 2.736656e+00 -1.510493e-01 1.106056e-03 4.781924e-01 4.320019e-02
      choose point on surface 2 2.743331e+00 -1.159358e-01 -3.089706e-04 9.407925e-03 3.753841e-01
      choose point on surface 6 2.762835e+00 -6.832717e-02 1.421712e-04 7.253652e-01 2.675979e-01
      choose point on surface 2 2.799684e+00 -2.745385e-02 -2.274045e-04 1.023971e-02 6.110983e-01
      choose point on surface 6 2.844728e+00 -4.576368e-03 5.542039e-04 5.140550e-02 4.388778e-01
      choose point on surface 2 2.888725e+00 3.510975e-03 -1.194887e-03 1.055568e-01 2.284118e-01
      choose point on surface 2 2.931267e+00 3.618484e-03 -1.177694e-03 1.525775e-01 1.399015e-01
      choose point on surface 2 3.054499e+00 -2.917317e-03 -1.033690e-03 1.706604e-01 5.944569e-01
      choose point on surface 2 3.102185e+00 -8.671826e-03 -4.604870e-04 6.391208e-02 5.088326e-02
      choose point on surface 2 3.154793e+00 -2.831301e-02 -1.495935e-03 9.958754e-02 5.087137e-01
      choose point on surface 2 3.198982e+00 -6.617821e-02 -6.060278e-04 6.017430e-01 3.721354e-01
      choose point on surface 2 3.220721e+00 -1.025024e-01 -1.951789e-03 5.575616e-01 3.278702e-01
      choose point on surface 2 3.232818e+00 -1.458436e-01 -1.255885e-03 4.221257e-02 2.775204e-01
      choose point on surface 2 3.233592e+00 -2.036583e-01 -1.554692e-03 4.892464e-01 4.375964e-01
      done
      marching options...
      forward and backward
      BC: left (forward) outward splay
      BC: right (forward) outward splay
      BC: left (backward) outward splay
      BC: right (backward) outward splay
      outward splay .5 .5 (left, right for outward splay BC)
      normal blending 6 6 (lines: left, right)
      lines to march $lines $lines (forward,backward)  
      points on initial curve 29
      generate
      smoothing...
        GSM:BC: top smoothed
        GSM:BC: bottom smoothed
        GSM:BC: right smoothed
        GSM:BC: left smoothed
        GSM:smooth grid
      name wireCapOuterSurface
      exit
#
#   Wire volume grid 
#
    active grid:spiralWireSurface
    create volume grid...
      backward
      $linesToMarch = 7; $dsn = .025;
      if( $order eq 4 ){ $dsn = .025*.7; $linesToMarch=9; }
      target grid spacing -1, $dsn (tang,normal, <0 : use default)      
      lines to march $linesToMarch
      geometric stretch factor 1.05
      points on initial curve 63, 669
      generate
      # open graphics
      name spiralWire
    exit
#
#    Wire cap outer volume grid
#
    active grid:wireCapOuterSurface
    create volume grid...
      target grid spacing -1, $dsn (tang,normal, <0 : use default)      
      lines to march $linesToMarch
      geometric stretch factor 1.05
#
      # lines to march 7
      # geometric stretch factor 1.05
      backward
      generate
      name wireCapOuter
      exit
#
#    Wire cap inner volume grid
#
    active grid:wireCapInnerSurface
    create volume grid...
      target grid spacing -1, $dsn (tang,normal, <0 : use default)      
      lines to march $linesToMarch
      geometric stretch factor 1.05
      # lines to march 7
      # geometric stretch factor 1.05
      backward
      generate
      name wireCapInner
      exit
#
#
#
    assign BC and share values
      boundary condition: 5
      shared boundary flag: 5
      plot lines on non-physical boundaries 0
      set BC and share 2 0 2 5 5
      set BC and share 1 0 2 5 5
      set BC and share 0 0 2 5 5
    exit
#     build a box grid
#       x bounds: -3.5 3.5
#       y bounds: -3.5 3.
#       z bounds: -.5 .5
#       lines: 281, 221, 41
#       name box1
#       bc 1 1 1 1 1 1 (l r b t b f)
#       share 0 0 0 0 0 0 (l r b t b f)
#     exit
  exit
#
# Here is the back ground grid
#
Box
  set corners
  #     $xa=-4.; $xb=4.; $ya=-4.5; $yb=4.; $za=-1.; $zb=1.; 
    $xa=-2.; $xb=2.; $ya=-3.5; $yb=4.; $za=-3.5; $zb=3.5; 
    $xa $xb $ya $yb $za $zb
  lines
    $nx = intmg( ($xb-$xa)/$ds +1.5);
    $ny = intmg( ($yb-$ya)/$ds +1.5);
    $nz = intmg( ($zb-$za)/$ds +1.5);
    $nx $ny $nz
  boundary conditions
    1 2 -1 -1 -1 -1 
  mappingName
    backGround
  exit
#
# Subroutine to convert to Nurbs Mapping and rotate/shift
#
#   NOTE: Always first rotate about y and x axes to make channel axis the x-axis
#      1. rotate 90 degrees about y axis
#      2. rotate $angleXaxis about the x-axis
#
#  ROTATE about the point (x0,y0,z0)
#
sub convertToNurbs\
{ local($old,$new,$angle,$rotationAxis,$x0,$y0,$z0,$angleXaxis,$share)=@_; \
  $cmds = "nurbs \n" . \
   "interpolate from mapping with options\n" . \
   " $old \n" . \
   " parameterize by index (uniform)\n" . \
   " number of ghost points to include\n $numGhostNurbs\n" . \
   " choose degree\n" . \
   "  3 \n" . \
   " # number of points to interpolate\n" . \
   " #  11 21 5 \n" . \
   "done\n" . \
   "# First rotate about y and x axes to make channel axis the x-axis \n" . \
   "rotate \n" . \
   " 90 1 \n" . \
   " 0. 0. 0.\n" . \
   "rotate \n" . \
   " $angleXaxis 0 \n" . \
   " 0. 0. 0.\n" . \
   "# user defined additional rotation: \n" . \
   "rotate \n" . \
   " $angle $rotationAxis \n" . \
   " $x0 $y0 $z0\n" . \
   "share \n" . \
   "  0 0 0 0 $share 0 \n" . \
   "mappingName\n" . \
   " $new\n" . \
   "exit"; \
}
#
# -------- CONVERT TO NURBS ------
$theta1=0.; # rotate spiral by this many degrees
$angle=$theta1; $rotationAxis=2; $x0=0.; $y0=0.; $z0=0.; $angleXaxis=90; $share=1; 
convertToNurbs(spiralWire,spiralWire1,$angle,$rotationAxis,$x0,$y0,$z0,$angleXaxis,$share);
$cmds
convertToNurbs(wireCapOuter,wireCapOuter1,$angle,$rotationAxis,$x0,$y0,$z0,$angleXaxis,$share);
$cmds
convertToNurbs(wireCapInner,wireCapInner1,$angle,$rotationAxis,$x0,$y0,$z0,$angleXaxis,$share);
$cmds
exit
#
#  Generate the overlapping grid 
#
generate an overlapping grid
  backGround
  spiralWire1
  wireCapOuter1
  wireCapInner1
  done choosing mappings
# 
  change parameters
    * improve quality of interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
  exit
 # change the plot
 # open graphics
 compute overlap
exit
* save an overlapping grid
save a grid (compressed)
$name
spiralWireGrid
exit
