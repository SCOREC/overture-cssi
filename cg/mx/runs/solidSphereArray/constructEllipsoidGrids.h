#
#   construct the ellipsoid grids
#
# -- master hyperbolic grids --
$ellipsoidInner="ellipsoidInnerMaster$masterCount"; 
$southPoleInner="southPoleInnerMaster$masterCount";
$northPoleInner="northPoleInnerMaster$masterCount";
$ellipsoidOuter="ellipsoidOuterMaster$masterCount"; 
$southPoleOuter="southPoleOuterMaster$masterCount";
$northPoleOuter="northPoleOuterMaster$masterCount";
# 
$xShift=$xe[$count]; $yShift=$ye[$count]; $zShift=$ze[$count]; $angle1=$rot1[$count];  $angle2=$rot2[$count];
$rotationAxis1=0;  $rotationAxis2=1;
#
include $ENV{CG}/mx/runs/solidSphereArray/buildEllipsoidGridsFromMaster.h
$count=$count+1;

