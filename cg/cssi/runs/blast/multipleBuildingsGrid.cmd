#***************************************************
# ogen command file: multiple buildings
#   Usage:
#     ogen [-noplot] multipleBuildingsGrid.cmd -wakeGrid=[0|1]
#
#    -wakeGrid=[0|1] : 1=add a grid in the wake
#
#  ogen -noplot multipleBuildingsGrid.cmd -prefix=eightBuildingsWithWake 
#  ogen -noplot multipleBuildingsGrid.cmd -prefix=eightBuildings -wakeGrid=0 -factor=1
#  ogen -noplot multipleBuildingsGrid.cmd -prefix=eightBuildings -wakeGrid=0 -factor=2  [ 4.7 M pts]
#
#***************************************************
#* Boundary conditions for cgins:
#   1= noSlipWall
#   2= slipWall
#   3=inflow
#   4=outflow
#**************************************************************************
# scale number of grid points in each direction by the following factor
# $factor=1;
# Here we get twice as many points:
#
$prefix="multipleBuildingsGrid";     # prefix for grid name 
$order=2; $factor=1; $interp = "e";  $ml=0; # default values
$orderOfAccuracy = "second order"; $ng=2; $wakeGrid=1; 
$factor=2.**(1./3.); printf(" factor=$factor\n");
#
# get command line arguments
GetOptions( "prefix=s"=>\$prefix,"order=i"=>\$order,"factor=f"=> \$factor,"blf=f"=>\$blf,"refinementBox=i"=>\$refinementBox,\
            "interp=s"=> \$interp,"name=s"=> \$name,"ml=i"=>\$ml,"sharpnessLB=f"=> \$sharpnessLB,\
            "widthX=f"=> \$widthX,"widthY=f"=> \$widthY,"widthZ=f"=> \$widthZ,\
            "xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,"za=f"=>\$za,"zb=f"=>\$zb,\
            "xac=f"=>\$xac,"xbc=f"=>\$xbc,"yac=f"=>\$yac,"ybc=f"=>\$ybc,"zac=f"=>\$zac,"zbc=f"=>\$zbc,\
            "wakeGrid=s"=> \$wakeGrid );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=3; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=4; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }else{ $interpType = "implicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
#
# Define a subroutine to convert the number of grid points
sub getGridPoints\
{ local($n1,$n2,$n3)=@_; \
  $nx=int(($n1-1)*$factor+1.5); $ny=int(($n2-1)*$factor+1.5); $nz=int(($n3-1)*$factor+1.5);\
}
#
#*************************************************************************
#
create mappings 
# 
#
#************************************************************************
#   Make the roundedCylinder buildings
#************************************************************************
#
include $ENV{'CG'}/cns/runs/blast/buildRoundedCylinder.h
#
#************************************************************************
#   Make the poly-building - using a smoothedPolygon as the cross-section
#************************************************************************
#
include $ENV{'CG'}/cns/runs/blast/buildPolyBuilding.h
#
#**************************************************************************
#   Now take the basic building and scale/shift it to create new buildings
#**************************************************************************
#
# ============== rounded building 1 ==============================
  rotate/scale/shift
    transform which mapping?
    roundedCylinderGrid
    shift
      .25 0 0.
    scale
     .5 1. 1.
   mappingName
    roundedCylinderGrid1
  exit
#
  rotate/scale/shift
    transform which mapping?
    roundedCylinderTop
    shift
      .25 0 0.
    scale
     .5 1. 1.
   mappingName
    roundedCylinderTop1
  exit
#
# ============== rounded building 2 ==============================
  rotate/scale/shift
    transform which mapping?
    roundedCylinderGrid
    shift
      .85 0. .6 
    scale
     1. 1.5 .75 
   mappingName
    roundedCylinderGrid2
  exit
#
  rotate/scale/shift
    transform which mapping?
    roundedCylinderTop
    shift
      .85 0. .6 
    scale
     1. 1.5 .75 
   mappingName
    roundedCylinderTop2
  exit
# ============== rounded building 3 ==============================
  rotate/scale/shift
    transform which mapping?
    roundedCylinderGrid
    shift
      .85 0. -.6 
    scale
     1. 1.25 1.
   mappingName
    roundedCylinderGrid3
  exit
#
  rotate/scale/shift
    transform which mapping?
    roundedCylinderTop
    shift
      .85 0. -.6 
    scale
     1. 1.25 1.
   mappingName
    roundedCylinderTop3
  exit
# ===================================================================
# ============== poly building 1 ==============================
  rotate/scale/shift
    transform which mapping?
    polyBuilding
    shift
      -.45 0 .5 
    scale
     1. .75 .75 
   mappingName
    polyBuilding1
  exit
#
  rotate/scale/shift
    transform which mapping?
    polyTopBox
    shift
      -.45 0 .5 
    scale
     1. .75 .75 
   mappingName
    polyTopBox1
  exit
# ============== poly building 2 ==============================
  rotate/scale/shift
    transform which mapping?
    polyBuilding
    shift
      -.75 0 -1.0 
    scale
     .55 1.25 .55  
   mappingName
    polyBuilding2
  exit
#
  rotate/scale/shift
    transform which mapping?
    polyTopBox
    shift
      -.75 0 -1.0 
    scale
     .55 1.25 .55  
   mappingName
    polyTopBox2
  exit
# ============== poly building 3 ==============================
  rotate/scale/shift
    transform which mapping?
    polyBuilding
    shift
      .25 0 1.25
    scale
     .75 1. .75    
   mappingName
    polyBuilding3
  exit
#
  rotate/scale/shift
    transform which mapping?
    polyTopBox
    shift
      .25 0 1.25  
    scale
     .75 1. .75 
   mappingName
    polyTopBox3
  exit
# ============== poly building 4 ==============================
  rotate/scale/shift
    transform which mapping?
    polyBuilding
    shift
      3.00 0 -.5 
    scale
     .55 1.25 1. 
   mappingName
    polyBuilding4
  exit
#
  rotate/scale/shift
    transform which mapping?
    polyTopBox
    shift
      3.00 0 -.5 
    scale
     .55 1.25 1. 
   mappingName
    polyTopBox4
  exit
#
#
# ==================================================================
#   ** build the tower **
include $ENV{'CG'}/cns/runs/blast/buildTower.h
#
#
# Now shift and scale the tower 
#
  rotate/scale/shift
    transform which mapping?
    towerPod
    scale
     .75 .75 .75 
    shift
      -.95 0 -.125
   mappingName
    towerPod1
  exit
  rotate/scale/shift
    transform which mapping?
    towerPodTop
    scale
     .75 .75 .75 
    shift
      -.95 0 -.125
   mappingName
    towerPodTop1
  exit
  rotate/scale/shift
    transform which mapping?
    tower
    scale
     .75 .75 .75 
    shift
      -.95 0 -.125
   mappingName
    tower1
  exit
  rotate/scale/shift
    transform which mapping?
      towerSpike
    scale
     .75 .75 .75 
    shift
      -.95 0 -.125
   mappingName
    towerSpike1
  exit
  rotate/scale/shift
    transform which mapping?
      towerSpikeCap
    scale
     .75 .75 .75 
    shift
      -.95 0 -.125
   mappingName
    towerSpikeCap1
  exit
# ===================================================================
#
#
# Here is the fine box around the buildings
#
Box
  set corners
#     -.5 1.5 0. 2.0 -.5 1.5 
   $xa=-2; $xb=2; $ya=0; $yb=2.5; $za=-1.5; $zb=1.5; # fine grid around buildings 
   if( $wakeGrid eq "0" ){ $xa=-1.75; $xb=2.25; }
   $xa $xb $ya $yb $za $zb
   # -2. 2. 0. 2.5 -1.5 1.5 
  lines
    # 81 65 65
    getGridPoints(81,65,65);
    $nx $ny $nz
  boundary conditions
    if( $wakeGrid eq "1" ){ $cmd="3 0 1 2 2 2"; }else{ $cmd="3 4 1 2 2 2"; }
    $cmd 
    # 3 0 1 2 2 2
  share
    0 0 2 0 0 0
  mappingName
    backGround
  exit
#
# Here are coarser boxes to extend the domain
#
Box
  set corners
   2.00  5.00   0. 2.5  -1.5 1.5 
  lines
 # 33 33 33
    getGridPoints(33,33,33);
    $nx $ny $nz
  boundary conditions
    0 4 1 2 2 2
  share
    0 0 2 0 0 0
  mappingName
    backGround1
  exit
#*
exit
generate an overlapping grid
  if( $wakeGrid eq "1" ){ $cmd="backGround1"; }else{ $cmd="#"; }
  $cmd
  # backGround1
  backGround
  roundedCylinderTop1
  roundedCylinderGrid1
  roundedCylinderTop2
  roundedCylinderGrid2
  roundedCylinderTop3
  roundedCylinderGrid3
#
  polyBuilding1
  polyTopBox1
  polyBuilding2
  polyTopBox2
  polyBuilding3
  polyTopBox3
  polyBuilding4
  polyTopBox4
#
  tower1
  towerPod1
#  towerPodTop1
  towerSpike1
  towerSpikeCap1
  done
#
  change the plot
    toggle grid 0 0
    toggle grid 1 0
    plot block boundaries (toggle) 1
   exit this menu
# display intermediate results
  change parameters
    ghost points
      all
      2 2 2 2 2 2
  exit
  # pause
  compute overlap 
  # pause
  exit
save a grid (compressed)
$name
multipleBuildingsGrid
exit

