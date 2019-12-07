#
#   Solid objects in a big domain
#
#
# usage: ogen [noplot] solidObjectsGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] -numGhost=<i>
# 
# examples:
#     ogen -noplot solidObjectsGrid -interp=e -order=2 -factor=4 
#     ogen -noplot solidObjectsGrid -interp=e -order=2 -factor=8
#     ogen -noplot solidObjectsGrid -interp=e -order=4 -factor=4 
#     ogen -noplot solidObjectsGrid -interp=e -order=4 -factor=8
#
$prefix="solidObjectsGrid"; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$name=""; $xa=-2.5; $xb=3.5; $ya=-2.5; $yb=2.5; $ml=0; 
$numGhost=-1;  # if this value is set, then use this number of ghost points
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=> \$xa,"xb=f"=> \$xb,"ya=f"=> \$ya,"yb=f"=> \$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"numGhost=i"=> \$numGhost);
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=3; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $numGhost ne -1 ){ $ng = $numGhost; } # overide number of ghost
if( $numGhost ne -1 ){ $suffix .= ".ng$numGhost"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.05/$factor;
#
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
create mappings
#
#  Coarse background 
#
#-  rectangle
#-    mappingName
#-      backGround
#-   set corners
#-    $dsc=2.*$ds; # coarse spacing 
#-    $width=4;
#-    $xac=$xa-$width; $xbc=$xb+$width; $yac=$ya-$width; $ybc=$yb+$width; 
#-    $xac $xbc $yac $ybc
#-   lines
#-    $nx = int( ($xbc-$xac)/$dsc +1.5 ); 
#-    $ny = int( ($ybc-$yac)/$dsc +1.5 ); 
#-    $nx $ny
#-  exit
  # Inner refinement 
  rectangle
    mappingName
     backGround
   set corners
    $xa $xb $ya $yb
   lines
    $nx = int( ($xb-$xa)/$ds +1.5 ); 
    $ny = int( ($yb-$ya)/$ds +1.5 ); 
    $nx $ny
    boundary conditions
      # 0 0 0 0
      1 2 -1 -1 
  exit
#
$gridNames="#"; 
#
# -----------------------------------------
# -------- build I-beam grid --------------
# -----------------------------------------
#
$iBeamNames="#"; 
$count=0;
$centerHeight=.75; $centerWidth=.25;
$edgeHeight=.25;   $edgeWidth=1.;
# -- we only build the hyperbolic grids once -- they are rotated/translated below 
include iBeamSolid.h
# 
# ---------------------- IBEAM 1 -------------------------------
$xr=.0; $yr=0; # center for rotation
$xShift=2; $yShift=0; $angle=0; $xMin=-$radX; $xMax=-$xMin; $yMin=-$radY; $yMax=-$yMin;
include buildIBeamGrids.h
# 
# ---------------------- IBEAM 2 ------------------------
$xShift=0; $yShift=1.5; $angle=90; $xMin=-$radY; $xMax=-$xMin; $yMin=-$radX; $yMax=-$yMin;
include buildIBeamGrids.h
# 
# ---------------------- IBEAM 3 ------------------------
$xShift=0; $yShift=-1.5; $angle=90; $xMin=-$radY; $xMax=-$xMin; $yMin=-$radX; $yMax=-$yMin;
include buildIBeamGrids.h
# 
# open graphics
#   
# ---------------------------------------
# ------- build base split-ring ---------
# ---------------------------------------
#
$splitRingNames="#"; 
$count=0;
include splitRingSolid.h
#
# ------------------ Split Ring 1 ---------------
$xr=.0; $yr=0; # center for rotation
$xShift=0;  $yShift=0; $angle=0; $xMin=-$radX; $xMax=-$xMin; $yMin=-$radY; $yMax=-$yMin;
include buildSplitRingGrids.h 
#
# ------------------ Split Ring 2 ---------------
$xr=.0; $yr=0; # center for rotation
$xShift=2; $yShift=1.5; $angle=90; $xMin=-$radY; $xMax=-$xMin; $yMin=-$radX; $yMax=-$yMin;
include buildSplitRingGrids.h 
#
# ------------------ Split Ring 3 ---------------
$xr=.0; $yr=0; # center for rotation
$xShift=2; $yShift=-1.5; $angle=-90; $xMin=-$radY; $xMax=-$xMin; $yMin=-$radX; $yMax=-$yMin;
include buildSplitRingGrids.h 
#
#
exit 
#
generate an overlapping grid
  backGround
#  innerBackGround
  $gridNames
  $iBeamNames
  $splitRingNames
  done
  change parameters
    # ---- Specify domains -----
    specify a domain
      # domain name:
      outerDomain 
      # grids in the background domain:
      backGround
      $gridNames
     done
#
    specify a domain
      # domain name:
      iBeamDomain
      # grids in the domain:
      $iBeamNames
    done
#
    specify a domain
      # domain name:
      splitRingDomain
      # grids in the domain:
      $splitRingNames
    done
#
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
    exit
   # pause
  # open graphics
  compute overlap
exit
#
save an overlapping grid
$name
triangle
exit
