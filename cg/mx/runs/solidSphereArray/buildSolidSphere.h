#
#  Output:
#   $innerNames
#   $outerNames
#
$sphereName="outerSphere$count"; 
$northPoleName="outerNorthPole$count";
$southPoleName="outerSouthPole$count"; 
#
$outerNames = "$sphereName\n $northPoleName\n $southPoleName"; 
#
$sphereShare=100 + 3*($count-1); 
#
# ------------- BUILD OUTER SPHERES -------------------------
$radiusDir=1; 
include $ENV{Overture}/sampleGrids/solidSphereInABox.h
#
#   ******** inner-sphere ********
$sphereName="innerSphere$count"; 
$northPoleName="innerNorthPole$count";
$southPoleName="innerSouthPole$count"; 
$innerBoxName = "innerBox$count"; 
#
$innerNames = "$innerBoxName\n $sphereName\n $northPoleName\n $southPoleName"; 
#
$sphereShare=100 + 3*($count-1);  # reset this so the inner sphere has the same corresponding share values
# 
# ------------- BUILD INNER SPHERES -------------------------
$radiusDir=-1;
include $ENV{Overture}/sampleGrids/solidSphereInABox.h
#
# Here is the inner box
#
Box
  set corners
    $dSphere = ( $innerRad + ($order-2 + $dse*1.5 )*$ds ); # size of box to enclose inner sphere
    $xa = $xSphere - $dSphere; $xb=$xSphere + $dSphere;
    $ya = $ySphere - $dSphere; $yb=$ySphere + $dSphere;
    $za = $zSphere - $dSphere; $zb=$zSphere + $dSphere;
    # 
    $xa $xb $ya $yb $za $zb
  lines
    $nx = int( ($xb-$xa)/$ds +1.5);
    $ny = int( ($yb-$ya)/$ds +1.5);
    $nz = int( ($zb-$za)/$ds +1.5);
    $nx $ny $nz
  boundary conditions
    0 0 0 0 0 0
  mappingName
    $innerBoxName
  exit
#**********************************
