#
#   Create Hyperbolic grids for the ellipsoid with semi-axes ($a,$b,$c) 
#
#
# ------ Create Ellipsoid Surface Mappings for Hyperbolic Volume Grids ----
#            Ellipsoid is defined by ($a,$b,$c) 
$masterCount = $masterCount+1;
$interfaceShare=100 + $masterCount;   # share flag for matching interfaces
# 
$a = $ae[$masterCount]; $b = $be[$masterCount]; $c = $ce[$masterCount];
printf("buildMasterEllipsoidGrids: masterCount=$masterCount (a,b,c)=($a,$b,$c)\n");
#
include $ENV{CG}/mx/runs/solidSphereArray/ellipsoidSurfacePatch.h
#
#
$rFactor=.75; 
$rDist = $rFactor*($nr-1)*$ds; 	# note $rFactor -- to make radial spacing a bit smaller 
#
#    ---- Construct INNER hyperbolic grids for the ellipsoid body and two patches on the poles ----
#
$nrSave=$nr;
$nr = $nrMin + ($order-2);
if( $interp eq "e" ){ $nr = $nrMin + $order +$nrExtra; }
$rDist = $rFactor*($nr-1)*$ds; 	# note $rFactor -- to make radial spacing a bit smaller 
$directionToMarch="backward"; 
#
$ellipsoidBodyName     ="ellipsoidInnerMaster$masterCount";
$ellipsoidSouthPoleName="southPoleInnerMaster$masterCount";
$ellipsoidNorthPoleName="northPoleInnerMaster$masterCount";
# 
include $ENV{CG}/mx/runs/solidSphereArray/ellipsoidVolumeGrids.h
#
#    ---- Construct OUTER hyperbolic grids for the ellipsoid body and two patches on the poles ----
#
$nr=$nrSave; 
$rDist = $rFactor*($nr-1)*$ds; 	# note $rFactor -- to make radial spacing a bit smaller 
$directionToMarch="forward"; 
#
$ellipsoidBodyName     ="ellipsoidOuterMaster$masterCount";
$ellipsoidSouthPoleName="southPoleOuterMaster$masterCount";
$ellipsoidNorthPoleName="northPoleOuterMaster$masterCount ";
# 
include $ENV{CG}/mx/runs/solidSphereArray/ellipsoidVolumeGrids.h 
#
