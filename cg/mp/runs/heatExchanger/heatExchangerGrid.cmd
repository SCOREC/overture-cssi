#**************************************************************************
#
#  Create a grid for a simple heat-exchanger:
#     A curved tube (containing a fluid) that passes through a solid block. 
#     This grid is used by cgmp to solve a thermal hydraulics fluid-solid heat transfer problem.
#
#  Usage:
#    ogen noplot heatExchangerGrid -factor=<> -interp=[e,i] -order=[2,4,6,8] -match=[yes,no]
#
# ogen -noplot heatExchangerGrid -factor=1 
# ogen -noplot heatExchangerGrid -factor=2 -interp=e
# ogen -noplot heatExchangerGrid -factor=4 -interp=e
#
#**************************************************************************
$prefix= "heatExchangerGrid"; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; 
$xa =-2.; $xb=2.; $ya=0.; $yb=3.; $za=-.5; $zb=.5; 
$radius=.20;          # radius of the tube
$nr=6;                # number of lines in the radial direction 
$pi=4.*atan2(1.,1.); 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"prefix=s"=> \$prefix );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids";  $nr=$nr+1; }  # add extra grid lines to nr for explicit interp. 
# 
$suffix = ".order$order"; 
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
#
$ds=.05/$factor;    # target grid spacing for fluid
$dss=$ds*1.5;       # target grid spacing for the solid (make larger)
$rDist=($nr-3)*$ds; # normal distance for fluid grids  (allow extra points for stretching)
$rDists=($nr-3)*$ds; # normal distance for solid grids  (allow extra points for stretching)
# 
$bcInflow=10;       # bc for inflow boundaries on the tube
$bcOutflow=11;      # bc for outflow boundaries on the tube
$bcInterface1=100;  # bc for interfaces are numbered starting from 100 
# 
create mappings
#
# 
# - define the cross-section for the solid tube --
  annulus
    inner and outer radii
      $deltaRad=$rDists; $outerRad=$radius+$deltaRad+$dss; 
      $radius $outerRad
    lines
      $nTheta= int( 2.*$pi*$radius/$dss +1.5 );
      $nTheta $nr
    boundary conditions
      -1 -1 1 0
    mappingName
      tube-outer-cross-section-unstretched
    exit
# -- cluster grid points near the tube wall 
  stretch coordinates
    Stretch r2:itanh
    $dsWall=$dss*.25; 
    STP:stretch r2 itanh: position and min dx 0 $dsWall
    stretch grid
    STRT:name tube-outer-cross-section
    exit
# 
# - define the cross-section for the fluid tube --
  annulus
    inner and outer radii
      $deltaRad=$rDist; $innerRad=$radius-$deltaRad; 
      $innerRad $radius
    lines
 # $nTheta= int( 2.*$pi*($radius+$innerRad)*.5/$ds +1.5 );
      $nTheta $nr
    boundary conditions
      -1 -1 0 1
    mappingName
      tube-inner-cross-section-unstretched
    exit
# -- cluster grid points near the tube wall 
  stretch coordinates
    Stretch r2:itanh
    $dsWall=$dss*.25; 
    STP:stretch r2 itanh: position and min dx 1. $dsWall
    stretch grid
    STRT:name tube-inner-cross-section
    exit
#
# --- Here is the curve that defines center line of the pipe ---
#       The end of the pipe should straighten out to meet the end walls
# 
  nurbs (curve)
    set range dimension
    3
    enter points 
       $txa=-1.; $txb=1.; $tya=0.; $tyb=2.; 
       11
      $txa .0  0. 
      $txa .1  0. 
      $txa .4  0. 
      $txa .8  0. 
      $txa 1.2  0. 
       0  $tyb  0. 
       $txb 1.2 0. 
       $txb .8  0. 
       $txb .4  0. 
       $txb .1  0. 
       $txb .0  0. 
# 
    lines
      101
    mappingName
      sweep-curve
   exit
#
#  -- outer grid on the tube --
# 
  sweep
    choose reference mapping
     tube-outer-cross-section
    choose sweep curve
      sweep-curve
    lines
      $arcLength = $txb-$txa+ 2.*($tyb-$tya); # guess arclength of the sweep-curve 
      $nSweep = int( $arcLength/$ds +1.5 );
      $nTheta $nr $nSweep
    boundary condition
      -1 -1 $bcInterface1 0 2 2
    share 
      0 0 $bcInterface1 0 2 2
    mappingName
      tube-outer
   exit
# 
#  -- inner grid on the tube --
# 
  sweep
    choose reference mapping
     tube-inner-cross-section
    choose sweep curve
      sweep-curve
    lines
      $nSweep = int( $arcLength/$ds +1.5 );
      $nTheta $nr $nSweep
    boundary condition
      -1 -1 0 $bcInterface1 $bcInflow $bcOutflow
    share 
      0 0 0 $bcInterface1 2 2
    mappingName
      tube-inner
   exit
#
# -- make a grid to go down the core of the tube ---
# 
#  -- core cross section --
rectangle
  set corners
 # Note: we need to make the core grid larger for explicit interpolation
    if( $interp eq "e" ){ $dsPlus = 2.*$ds; }else{ $dsPlus=$ds; }
    $xac=-($innerRad+$dsPlus); $xbc=-$xac; $yac=$xac; $ybc=$xbc; 
    $xac $xbc $yac $ybc
  lines
    $nxc = int( ($xbc-$xac)/$ds +2.5 ); 
    $nyc = int( ($ybc-$yac)/$ds +2.5 ); 
    $nxc $nyc
  boundary conditions
    0 0 0 0 
  mappingName
    tube-core-cross-section
exit
# 
#  -- tube core grid --
# 
  sweep
    choose reference mapping
     tube-core-cross-section
    choose sweep curve
      sweep-curve
    lines
      $nxc $nyc $nSweep
    boundary condition
      0 0 0 0 $bcInflow $bcOutflow
    share 
      0 0 0 0 2 2 
    mappingName
      tube-core
   exit
#
Box
  set corners
    $xa $xb $ya $yb $za $zb
  lines
    $nx = int( ($xb-$xa)/$dss +1.5);
    $ny = int( ($yb-$ya)/$dss +1.5);
    $nz = int( ($zb-$za)/$dss +1.5);
    $nx $ny $nz
  boundary conditions
    1 1 2 2 3 3 
  share
    0 0 2 0 0  0  
  mappingName
    solid
  exit
exit
#
# ----- generate the overlapping grid ---------
#
generate an overlapping grid
  solid
  tube-outer
# 
  tube-core
  tube-inner
  done choosing mappings
# 
  change parameters 
 # define the domains -- these will behave like independent overlapping grids
    specify a domain
 # domain name:
      solidDomain
 # grids in the domain:
      solid
      tube-outer
      done
    specify a domain
 # domain name:
      tubeDomain 
 # grids in the domain:
      tube-core
      tube-inner
      done
    order of accuracy
     $orderOfAccuracy
    interpolation type
      $interpType
    exit 
  compute overlap
 # pause
exit
# save an overlapping grid
save a grid (compressed)
$name
heatExchanger
exit
