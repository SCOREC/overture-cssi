#
# Create the initial grid for a flexible half-disk in a channel
# Use this grid with cgmp for a fluid-structure example.
#
# Usage:
#         ogen [noplot] halfDiskInAChannelGrid [options]
# where options are
#     -factor=<num>          : grid spacing factor
#     -interp=[e/i]          : implicit or explicit interpolation
#     -prefix=<string>       : over-ride the default prefix in the grid name
#     -radius=<num>          : radius of the half disk
#     -xa, xb, ya, yb        : channel dimensions
#     -beamX=<float>         : location of the beam, from the start of the channel
#
# Examples:
#
#    ogen -noplot halfDiskInAChannelGrid -interp=e -factor=2
# 
#
$prefix = "halfDiskInAChannelGrid";
$factor=1; $name="";
$interp="i"; $interpType = "implicit for all grids";
$order=2; $orderOfAccuracy = "second order"; $ng=2;
$nExtra=0;
$radius=.5;
$xa=-1.5; $xb=2.0;
$ya=0.; $yb=1.; 
#
# get command line arguments
GetOptions("name=s"=> \$name,"order=i"=>\$order,"factor=f"=> \$factor,"interp=s"=> \$interp,\
           "nExtra=i"=>\$nExtra,"width=f"=> \$width, "beamThickness=f"=> \$beamThickness, \
           "radius=f"=> \$radius, "xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb, \
           "tStretch=f"=> \$tStretch,"beamX"=>\$beamX,"prefix=s"=> \$prefix );
#
sub min{ local($n,$m)=@_; if( $n<$m ){ return $n; }else{ return $m; } }                  
#
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }else{ $interpType = "implicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
#
#
$ds0=.1;
$ds=$ds0/$factor;
$pi = 4.*atan2(1.,1.);
# 
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
#
#     --- create the individual mappings -----
#
create mappings 
#
  rectangle 
    $nx=int( ($xb-$xa)/$ds+1.5 );
    $ny=int( ($yb-$ya)/$ds+1.5 );
    set corners 
      $xa $xb $ya $yb
    lines 
      $nx,$ny 
    boundary conditions 
      1 2 3 4 
    share 
      0,0,3,0 
    mappingName
      backGroundFluid
    exit 
# 
#
# --- make the ellipse using a spline
#    NOTE: go in a counter clockwise direction
#
$radX=$radius; $radY=$radius; $cx=0.; $cy=0.; 
$scmd="#";
$ns=10*$factor; $arcLength=0.;
for( $i=1; $i<=$ns; $i++ ){ $s=$pi*($i-1.)/($ns-1.); $x=$radX*cos($s)+$cx; $y=$radY*sin($s)+$cy; \
   if( $i > 1 ){ $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 );} $x0=$x; $y0=$y; \
   $scmd .= "\n $x $y"; }
# 
spline
  #
  enter spline points
    $ns
    $scmd
  lines
    $ns
  mappingName
    curveBoundary
 exit
#
#  -- fluid interface grid ----
# 
  hyperbolic
    forward
    $nDistInterface=.2; 
    $nDistExtraFactor=1.2; # make grid a bit finer in the normal direction
    $nrm = intmg(($nDistInterface*$nDistExtraFactor)/$ds+.5);
    #
    distance to march $nDistInterface
    lines to march $nrm
    $nThetaInterface = intmg($arcLength/$ds+1.5);
    points on initial curve $nThetaInterface
    # volume smooths 50
    # This next line is important to keep the ghost lines
    # near the corner of good quality:
    apply boundary conditions to start curve 1
    generate 
    boundary conditions
      3 3 100 0 
    share
      3 3 100 0 
    boundary condition options...
    BC: left fix y, float x and z
    BC: right fix y, float x and z
    # normal blending 5 5 (lines: left, right)
    #
    generate 
    name diskInterfaceFluid
    exit
#
#
# -------------------- SOLID GRIDS ---------------------
#
# -- Solid half annulus ----
#
Annulus
  inner and outer radii
    $deltar =.2; 
    $rb = $radius - $deltar; # flip outside to inside to interface is at r2=0
    $innerRad=$radius; $outerRad=$rb; 
    $innerRad $outerRad
  angles: 0. .5
  lines
    $nTheta = $nThetaInterface;
    $nr = intmg( $nDistExtraFactor*($deltar)/$ds + 1.5 );
    $nTheta $nr
  boundary conditions
    3 3 100 5 
  share
     0  0 100 0 
  mappingName
    diskInterfaceSolid
  exit
# 
#  Solid background grid 
# 
#-   rectangle 
#-     $xar=-$radius; $xbr=-$xar; $yar=$ya; $ybr=$ya+$radius; 
#-     $nx=int( ($xbr-$xar)/$ds+1.5 );
#-     $ny=int( ($ybr-$yar)/$ds+1.5 );
#-     set corners 
#-       $xar $xbr $yar $ybr
#-     lines 
#-       $nx $ny 
#-     share 
#-       0,0,3,0 
#-     boundary conditions 
#-       0 0 3 0 
#-     mappingName
#-       backGroundSolid
#-     exit 
  exit this menu 
# 
generate an overlapping grid
  backGroundFluid
  diskInterfaceFluid
  #  backGroundSolid
  diskInterfaceSolid
  done
  #
  change parameters
    specify a domain
      solidDomain
      # backGroundSolid
      diskInterfaceSolid
      done
    specify a domain
      fluidDomain
      backGroundFluid
      diskInterfaceFluid
      done
    order of accuracy
     $orderOfAccuracy
    interpolation type
      $interpType
    ghost points
      all 
      $ng $ng $ng $ng $ng $ng
  exit
  # open graphics
  compute overlap
  exit
save an overlapping grid
$name
halfDiskInAChannel
exit
