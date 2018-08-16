#
# Create the initial grid for HALF a deforming disk in a channel
# Use this grid with cgmp for a fluid-structure example.
#
# Usage:
#         ogen [noplot] halfDiskDeform [options]
# where options are
#     -factor=<num>     : grid spacing is .1 divided by this factor
#     -interp=[e/i]     : implicit or explicit interpolation
#     -name=<string>    : over-ride the default name  
#     -case=[inner|outer] : only build a grid for the inner or outer domain
#     -nExtra          : add extra lines in the normal direction on the boundary fitted grids
#     -ae=<>, -be=<>   : scale factors for major and minor axes (default 1 for a circle)
#
# Examples:
#
#      ogen noplot halfDiskDeform -factor=.5
#      ogen noplot halfDiskDeform -factor=1
#      ogen noplot halfDiskDeform -factor=2
#      ogen noplot halfDiskDeform -factor=1 -interp=e
#      ogen noplot halfDiskDeform -factor=2 -interp=e
#      ogen noplot halfDiskDeform -factor=4 -interp=e
# -- fixed radial distance for grids next to the annulus
#     ogen noplot halfDiskDeform -fixedRadius=.25 -interp=e -factor=2 
#     ogen noplot halfDiskDeform -fixedRadius=.25 -interp=e -factor=4
# Bigger outer domain: 
#   ogen noplot halfDiskDeform -factor=1 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="halfDiskDeformBig1e.hdf"
#   ogen noplot halfDiskDeform -factor=2 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="halfDiskDeformBig2e.hdf"
#   ogen noplot halfDiskDeform -factor=4 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="halfDiskDeformBig4e.hdf"
#   ogen noplot halfDiskDeform -factor=8 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -nExtra=5 -name="halfDiskDeformBig8e.hdf"
#   ogen noplot halfDiskDeform -factor=16 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -nExtra=5 -name="halfDiskDeformBig16e.hdf"
#   ogen noplot halfDiskDeform -factor=32 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -nExtra=5 -name="halfDiskDeformBig32e.hdf"
#   ogen noplot halfDiskDeform -factor=64 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -nExtra=5 -name="halfDiskDeformBig64e.hdf"
#
# -- Bigger outer domain + fixed radial distance:
#    ogen -noplot halfDiskDeform -fixedRadius=.25 -factor=1 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="halfDiskDeformFixedBig1e.hdf"
#    ogen -noplot halfDiskDeform -fixedRadius=.25 -factor=2 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="halfDiskDeformFixedBig2e.hdf"
#    ogen -noplot halfDiskDeform -fixedRadius=.25 -factor=4 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="halfDiskDeformFixedBig4e.hdf"
#    ogen -noplot halfDiskDeform -fixedRadius=.25 -factor=8 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="halfDiskDeformFixedBig8e.hdf"
#    ogen -noplot halfDiskDeform -fixedRadius=.25 -factor=16 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="halfDiskDeformFixedBig16e.hdf"
# 
#  -- Non-matching grid spacing at the interface: (for testing non-matching interfaces)
#   ogen noplot halfDiskDeform -factor=2 -factor2=1 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="halfDiskDeformBig2x1e.hdf"
#   ogen noplot halfDiskDeform -factor=2 -factor2=3 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -name="halfDiskDeformBig2x3e.hdf"
#  -- add a refinement grid
#   ogen noplot halfDiskDeform -factor=2 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5 -refineInner=1 -name="halfDiskDeformBigRefineInner2e.hdf"
# outer-domain only:
#     ogen noplot halfDiskDeform -case=outer -prefix=halfDiskDeformOuter -interp=e -factor=1   [creates halfDiskDeformOutere1.hdf
#     ogen noplot halfDiskDeform -case=outer -factor=1 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5  -name=halfDiskDeformOuterBige1.hdf
#     ogen noplot halfDiskDeform -case=outer -factor=2 -interp=e -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.5  -name=halfDiskDeformOuterBige2.hdf
# 
#     ogen noplot halfDiskDeform -case=outer -factor=2 -interp=e -xa=-4. -xb=8.0 -ya=-3.5 -yb=3.5  -name=halfDiskDeformOuterBige2a.hdf
#
# -- Ellipse:
# ogen noplot halfDiskDeform -interp=e -xa=-3. -xb=3. -ya=-3. -yb=3. -ae=1.2 -be=.8 -factor=1 -name="ellipseDeformBig1e.hdf" 
# ogen noplot halfDiskDeform -interp=e -xa=-3. -xb=3. -ya=-3. -yb=3. -ae=1.2 -be=.8 -factor=2 -name="ellipseDeformBig2e.hdf" 
# ogen noplot halfDiskDeform -interp=e -xa=-3. -xb=3. -ya=-3. -yb=3. -ae=1.2 -be=.8 -factor=4 -name="ellipseDeformBig4e.hdf" 
# ogen noplot halfDiskDeform -interp=e -xa=-3. -xb=3. -ya=-3. -yb=3. -ae=1.2 -be=.8 -factor=8 -name="ellipseDeformBig8e.hdf" 
#
#
$ae=1.; $be=1.; 
$factor=1; $name=""; $case=""; 
$factor2=-1;   # by default factor2=factor
$interp="i"; $interpType = "implicit for all grids"; 
$order=2; $orderOfAccuracy = "second order"; $ng=2; 
$xa=-2.; $xb=2.; $ya=0.; $yb=2.; $nExtra=0; 
$refineInner=0; $refineOuter=0; $fixedRadius=-1; 
$prefix = "halfDiskDeform"; 
# 
# get command line arguments
GetOptions("name=s"=> \$name,"order=i"=>\$order,"factor=f"=> \$factor,"interp=s"=> \$interp,"case=s"=> \$case,\
           "xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,"nExtra=i"=>\$nExtra,"factor2=f"=> \$factor2,\
           "refineInner=i"=>\$refineInner,"refineOuter=i"=>\$refineOuter,"fixedRadius=f"=>\$fixedRadius,\
           "ae=f"=>\$ae,"be=f"=>\$be,"prefix=s"=> \$prefix );
#
if( $factor2 < 0 ){ $factor2=$factor; }
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
$suffix = ".order$order"; 
if( $fixedRadius ne -1 ){ $prefix .= "Fixed"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
#
$bcInterface0=100;  # bc for interfaces
$bcInterface1=101;  
$shareInterface=100;        # share value for interfaces
#
$Pi=4.*atan2(1.,1.);
#
$ds0 = .1; 
# target grid spacing:
$ds = $ds0/$factor;
$ds2 = $ds0/$factor2;
#
    $width=.25; if( $factor<1. ){ $width = .25*$factor; }
    $width2=.25; if( $factor2<1. ){ $width2 = .25*$factor2; }
    $rad=1.; 
# 
create mappings
#
  rectangle
    $nx=int( ($xb-$xa)/$ds+1.5 ); 
    $ny=int( ($yb-$ya)/$ds+1.5 ); 
    set corners
      $xa $xb $ya $yb 
    lines
      $nx $ny 
    boundary conditions
      1 2 3 4
    share
      0 0 3 0
    mappingName
      outerSquare
    exit
#
  rectangle
 # make the inner square bigger for deforming grid problems
 # $xb = $rad-$width/$factor + 2.*$ds ;
    $xb = $rad + .5;
    $xa=-$xb; $yb=$xb; 
    $nx=int( ($xb-$xa)/$ds2+1.5 ); 
    $ny=int( ($yb-$ya)/$ds2+1.5 ); 
    set corners
      $xa $xb $ya $yb 
    lines
      $nx $ny 
    boundary conditions
      0 0 3 0 
    share
      0 0 3 0 
    mappingName
      innerSquare
    exit
#
# Create a start curve for the interface
#
#    $innerRadius=$rad; $outerRadius=$innerRadius + $width/$factor; 
#    $averageRadius=($innerRadius+$outerRadius)/2.;
# 
  spline
    $n=51;
    enter spline points
      $n 
    $x0=0.; $y0=1.*$rad;
    $commands="";
    for( $i=0; $i<$n; $i++ ){ $theta=$Pi*$i/($n-1.); $x0=$ae*$rad*cos($theta); $y0=$be*$rad*sin($theta); \
                              $commands = $commands . "$x0 $y0\n"; }
    $commands
    lines
      $nr = int( 2.*$Pi*$rad/$ds+1.5 );
      $nr
 # pause
    exit
#  
   $width = $width + $ds0*$nExtra;
# 
  hyperbolic
 # add a few extra points as the boundary deforms it gets longer
    $stretchFactor=1.25; 
    $dist = $width/$factor;     
    $ns = int( $width/$ds0 +1.5 );  if( $ns<3 ){ $ns=3; }
    if( $fixedRadius ne -1 ){ $dist=$fixedRadius; $ns = int( $dist/$ds + 2.5 ); }
    distance to march $dist 
    lines to march $ns
    points on initial curve $nr
    BC: left fix y, float x and z
    BC: right fix y, float x and z
    generate
 # -- set the order of data point interpolation: 
    fourth order
 # second order
    boundary conditions
       3 3 $bcInterface0 0
    share
       3 3 $shareInterface 0 
    name outerInterface
 # pause
    exit
#  
  hyperbolic
 # add a few extra points as the boundary deforms it gets longer
    $stretchFactor=1.25; 
    $dist = $width2/$factor2;     
    $ns = int( $width2/$ds0 +1.5 ); if( $ns<3 ){ $ns=3; }
    if( $fixedRadius ne -1 ){ $dist=$fixedRadius; $ns = int( $dist/$ds + 2.5 ); }
    backward
    distance to march $dist 
    lines to march $ns
    $nr2 = int( 2.*$Pi*$rad/$ds2+1.5 ); 
    points on initial curve $nr2
    BC: left fix y, float x and z
    BC: right fix y, float x and z
    generate
 # -- set the order of data point interpolation: 
    fourth order
 # second order
    boundary conditions
       3 3 $bcInterface0 0
    share
       3 3 $shareInterface 0 
    name innerInterface
    exit
#
  $innerGrids="innerSquare\n innerInterface"; 
# 
#
  exit this menu
#
generate an overlapping grid
if( $case eq "inner" ){ $gridList=$innerGrids; }\
elsif( $case eq "outer" ){ $gridList="outerSquare\n outerInterface"; }\
else{ $gridList="outerSquare\n outerInterface\n $innerGrids";  }
  $gridList
  done choosing mappings
# 
  change parameters 
 # define the domains -- these will behave like independent overlapping grids
   if( $case eq "" ){ $cmd="specify a domain\n outerDomain\n outerSquare\n outerInterface\n done";}else{ $cmd="*"; }
      $cmd
   if( $case eq "" ){ $cmd="specify a domain\n innerDomain\n $innerGrids \n done";}else{ $cmd="*"; }
      $cmd
    order of accuracy
     $orderOfAccuracy
    interpolation type
      $interpType
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
    exit 
#
  # open graphics
  compute overlap
#
exit
#
save an overlapping grid
  $name
  halfDiskDeform
exit
