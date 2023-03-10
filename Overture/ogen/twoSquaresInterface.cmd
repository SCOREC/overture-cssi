#**************************************************************************
#
#    Grid for an interface calculation between two squares
#
# usage: ogen [-noplot] twoSquaresInterface -factor=<num> -order=[2/4/6/8] -interp=[e/i] -name= -bc=[d|p]
#              -angle=<f> -prefix=<s> -orientation=[horizontal|vertical]
# 
#    factor : grid resolution factor 
#    yFactor : if set, use this grid resolution factro in y 
# 
# Examples:
#        ogen -noplot twoSquaresInterface -factor=.25 -yFactor=.25 -interp=e -name="twoSquaresInterface0.hdf"
#        ogen -noplot twoSquaresInterface -factor=1 -yFactor=.5 -interp=e
#        ogen -noplot twoSquaresInterface -factor=2 -yFactor=.5 -interp=e
#        ogen -noplot twoSquaresInterface -factor=2 -yFactor=2 -interp=e
# 
#        ogen -noplot twoSquaresInterface -factor=4 -yFactor=4 -interp=e
#        ogen -noplot twoSquaresInterface -factor=16 -yFactor=16 -interp=e
#        ogen -noplot twoSquaresInterface -factor=32 -yFactor=32 -interp=e
#        ogen -noplot twoSquaresInterface -factor=64 -yFactor=64 -interp=e
# 
#        ogen -noplot twoSquaresInterface -interp=e -order=2 -factor=1 
#        ogen -noplot twoSquaresInterface -interp=e -order=2 -factor=2 
#        ogen -noplot twoSquaresInterface -interp=e -order=2 -factor=4 
#        ogen -noplot twoSquaresInterface -interp=e -order=2 -factor=8 
#        ogen -noplot twoSquaresInterface -interp=e -order=2 -factor=16
# 
#        ogen -noplot twoSquaresInterface -interp=e -order=4 -factor=1 
#        ogen -noplot twoSquaresInterface -interp=e -order=4 -factor=2 
#        ogen -noplot twoSquaresInterface -interp=e -order=4 -factor=4 
#        ogen -noplot twoSquaresInterface -interp=e -order=4 -factor=8 
#        ogen -noplot twoSquaresInterface -interp=e -order=4 -factor=16
# 
#        ogen -noplot twoSquaresInterface -factor=.5 -yFactor=.5 -interp=e -name="twoSquaresInterfacee0.order2"
# 
#  periodic in y:
#    ogen -noplot twoSquaresInterface -factor=2 -yFactor=.5 -bc=p -interp=e -name="twoSquaresInterfacenp2.hdf"
#    ogen -noplot twoSquaresInterface -factor=32 -yFactor=4 -bc=p -interp=e -name="twoSquaresInterfacenp32.hdf"
#
# rotated:
#    ogen -noplot twoSquaresInterface -factor=1 -order=2 -angle=45 -name="twoSquaresInterfaceRotated1.order2.hdf"
#    ogen -noplot twoSquaresInterface -factor=1 -order=4 -angle=45 -name="twoSquaresInterfaceRotated1.order4.hdf"
# 
# nonSquares:
#    ogen -noplot twoSquaresInterface -interp=e -angle=0.0 -prefix=twoNonSquaresInterface -order=4 -factor=2 
# non-matching:
#    ogen -noplot twoSquaresInterface -factor=1 -yFactor=.5 -yFactorRight=2. -name="twoSquaresInterface1to2.order2.hdf"
#
# with a refinement at the interface
#    ogen -noplot twoSquaresInterface -factor=1 -yFactor=1 -refineLeft=1 -name="twoSquaresInterface1RefineLeft.order2.hdf"
#    ogen -noplot twoSquaresInterface -factor=1 -yFactor=1 -refineLeft=1 -refineRight=1 -name="twoSquaresInterface1Refine.order2.hdf"
# 
#**************************************************************************
$order=2;  $orderOfAccuracy = "second order"; 
$prefix="twoSquaresInterface"; 
$interp="i"; $interpType = "implicit for all grids";
$order=2; $interp="i"; $name=""; $bc="d";
$xa=-1; $xb=1; # far left and far right
$ya=-1; $yb=1; # top and bottom for orientation=vertical
$ml=0;  # to-do: add MG levels
$numGhost=-1; # if >0 use this many ghost 
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$angle=0;
$factor=1; $yFactor=-1; 
$xFactorRight=1;  # additional grid resolution factor for the right domain (for non-matching grid lines)
$yFactorRight=1;  # additional grid resolution factor for the right domain (for non-matching grid lines)
$refineLeft=0; $refineRight=0; 
$orientation="horizontal"; # or vertical
# 
# get command line arguments
GetOptions("order=i"=>\$order,"factor=f"=> \$factor,"yFactor=f"=> \$yFactor,"interp=s"=> \$interp,\
           "outerRad=f"=> \$outerRad,"xa=f"=> \$xa,"xb=f"=> \$xb,"bc=s"=>\$bc,"angle=f"=> \$angle,\
           "name=s"=>\$name,"xFactorRight=f"=> \$xFactorRight,"yFactorRight=f"=> \$yFactorRight,\
           "refineLeft=i"=>\$refineLeft,"refineRight=i"=>\$refineRight,"numGhost=i"=>\$numGhost,\
           "prefix=s"=>\$prefix,"ml=i"=>\$ml,"orientation=s"=>\$orientation );
# 
if( $yFactor eq -1 ){ $yFactor=$factor; }
#
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
if( $bc eq p ){ $suffix = "p"; }else{ $suffix=""; } 
$suffix .= ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; } 
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $name eq "" ){ $name = $prefix . "$interp$factor" . $suffix . ".hdf"; }
#
# 
#
# $bcLeft = "1 100  1  1"; 
# $bcRight= "100 1  1  1"; 
# *new* Sept. 17, 2016: 
if( $orientation eq "horizontal" ){\
  $bcLeft = "1 100  3  4"; $bcRight= "100 6 7 8";  $shareLeft="0 100 0 0"; $shareRight="100 0 0 0"; }\
else{\
  $bcLeft = "1 2 3 100"; $bcRight= "5 6 100 8";  $shareLeft="0 0 0 100"; $shareRight="0 0 100 0"; }\
# 
if( $orientation eq "horizontal" ){\
  if( $bc eq "p" ){ $bcLeft ="1 100 -1 -1"; $bcRight="100 6 -1 -1"; $shareLeft="0 100 0 0"; $shareRight="100 0 0 0"; }\
}else{\
  if( $bc eq "p" ){ $bcLeft ="-1 -1 1 100"; $bcRight="-1 -1 100 6"; $shareLeft="0 0 0 100"; $shareRight="0 0 100 0";}\
    }
#
#**************************************************************************
#
#
# domain parameters:  
$dsx = .1/$factor; # target grid spacing in the x direction
$dsy= .1/$yFactor;          # target grid spacing in the y direction
#
create mappings
#
if( $angle eq "0" ){ $leftSquareBase="leftSquare"; $leftSquareRotated="leftSquareRotated"; }else{ $leftSquareBase="leftSquare0";  $leftSquareRotated="leftSquare";}
if( $angle eq "0" ){ $rightSquareBase="rightSquare"; $rightSquareRotated="rightSquareRotated"; }else{ $rightSquareBase="rightSquare0";  $rightSquareRotated="rightSquare";}
#
#  here is the left grid
#
  rectangle
    if( $orientation eq "horizontal" ){\
         $xar=$xa; $xbr=0.0; $yar= 0.; $ybr=1.; }else{\
         $xar=0.;  $xbr=1.0; $yar=$ya; $ybr=0.; }
    set corners
     $xar $xbr $yar $ybr 
    lines
      $nx=int( ($xbr-$xar)/$dsx+1.5 );
      $ny=int( ($ybr-$yar)/$dsy+1.5 );
      $nx $ny
    boundary conditions
      $bcLeft
    share
      # for now interfaces are marked with share>=100 
      $shareLeft
    mappingName
     $leftSquareBase 
     #  leftSquare0
     #*      leftSquare
  exit
#
# -- right grid ----
#
  rectangle
    if( $orientation eq "horizontal" ){\
         $xar=0.; $xbr=$xb; $yar= 0.; $ybr=1.; }else{\
         $xar=0.; $xbr=1.0; $yar=0.0; $ybr=$yb; }  
    # $xar= 0.0;  $xbr=$xb;
    # $yar= 0.;   $ybr=1.; 
    set corners
     $xar $xbr $yar $ybr 
    lines
      $dsx = $dsx/$xFactorRight;
      $dsy = $dsy/$yFactorRight;
      $nx=int( ($xbr-$xar)/$dsx+1.5 );
      $ny=int( ($ybr-$yar)/$dsy+1.5 );
      $nx $ny
    boundary conditions
      $bcRight
    share
      # for now interfaces are marked with share>=100 
      $shareRight
    mappingName
      $rightSquareBase
      # rightSquare0
      #* rightSquare
  exit
#
  rotate/scale/shift
    transform which mapping?
      $leftSquareBase
    rotate
     $angle
    0. .5 0.
    mappingName
      $leftSquareRotated
    exit
#
  rotate/scale/shift
    transform which mapping?
      $rightSquareBase
    rotate
     $angle
    0. .5 0.
    mappingName
     $rightSquareRotated
     # rightSquare
    exit
#
  reparameterize
    transform which mapping?
      $leftSquareBase 
      # leftSquare0
    set corners
      .5 1. .25 .75
    lines
      $nxr = $nx;  $nyr=$ny; 
      $nxr $nyr
    boundary conditions
     $shareLeft
    mappingName
      refinedLeftSquare
    exit
# 
    $leftGrids="leftSquare"; 
    if( $refineLeft eq 1 ){ $leftGrids .="\n refinedLeftSquare"; }
# 
#
  reparameterize
    transform which mapping?
      $rightSquareBase
      # rightSquare0
    set corners
      .0 .5 .1 .6
    lines
      $nxr = $nx;  $nyr=$ny; 
      $nxr $nyr
    boundary conditions
     $shareRight
    mappingName
      refinedRightSquare
    exit
# 
    $rightGrids="rightSquare"; 
    if( $refineRight eq 1 ){ $rightGrids .="\n refinedRightSquare"; }
# 
  exit this menu
#
generate an overlapping grid
  $leftGrids
  $rightGrids
  done 
  change parameters 
 # define the domains -- these will behave like independent overlapping grids
    specify a domain
 # domain name:
      leftDomain 
 # grids in the domain:
      $leftGrids
      done
    specify a domain
 # domain name:
      rightDomain 
 # grids in the domain:
      $rightGrids
      done
 # choose implicit or explicit interpolation
    interpolation type
      $interpType
    ghost points 
      all 
      $ng $ng $ng $ng
    order of accuracy
     $orderOfAccuracy
    exit 
# 
#    display intermediate results
# pause
  compute overlap
# pause
  exit
#
save an overlapping grid
  $name
  twoSquaresInterface
exit





# ------------ old way --------------

generate an overlapping grid
  leftSquare
  rightSquare
  done 
  change parameters
    prevent interpolation
      all
      all
    done
    order of accuracy
     $orderOfAccuracy
    ghost points
      all
      3 3 3 3 3 3 
  exit
# pause
  compute overlap
# pause
  exit
#
save an overlapping grid
  $name
  twoSquaresInterface
exit
