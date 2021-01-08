#
#
#  Construct ellipsoid grids for solidEllipsoidArrayGrid.cmd
#  Input:
#    $angle1, $rotationAxis, $angle2, $rotationAxis2, $xShift, $yShift, $zShift
#
$numGhost=$ng+1; # N.B. to avoid negative volumes in the ghost points interpolate ghost too in Nurbs.
#
# -- name the grids 
$ellipsoidInnerGrid="ellipsoidInner$numDomains"; 
$southPoleInnerGrid="southPoleInner$numDomains"; 
$northPoleInnerGrid="northPoleInner$numDomains"; 
$backGroundInner   ="backGroundInner$numDomains"; 
#
$ellipsoidOuterGrid="ellipsoidOuter$numDomains"; 
$southPoleOuterGrid="southPoleOuter$numDomains"; 
$northPoleOuterGrid="northPoleOuter$numDomains"; 
#
# Here is where we create the ellipsoid grids for this domain from the master ellipsoid grids
#
printf("convertToNurbs: source=$ellipsoidInner target=$ellipsoidInnerGrid\n");
#
convertToNurbs($ellipsoidInner,$ellipsoidInnerGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
convertToNurbs($southPoleInner,$southPoleInnerGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
convertToNurbs($northPoleInner,$northPoleInnerGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
convertToNurbs($ellipsoidOuter,$ellipsoidOuterGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
convertToNurbs($southPoleOuter,$southPoleOuterGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
convertToNurbs($northPoleOuter,$northPoleOuterGrid,$angle1,$rotationAxis1,$angle2,$rotationAxis2,$xShift,$yShift,$zShift);
$cmds
#
# Here is the INNER background box
#   **** WARNING : THIS MAY FAIL IF THE GRID IS ROTATED --- FIX ME ****
#
Box
  set corners
    # Choose bounds on inner background grid 
    $delta = 2*$ds; 
    $aEllipse=$a; $bEllipse=$b; $cEllipse=$c;
    # If the ellipsoid is rotated, make the inner background a big enough square
    if( $angle1 !=0 || $angle2 !=0 ){ $aEllipse=max($a,$b,$c); $bEllipse=$aEllipse; $cEllipse=$aEllipse; } # 
    $xai=-$aEllipse+$delta+$xShift; $xbi=$aEllipse-$delta+$xShift;
    $yai=-$bEllipse+$delta+$yShift; $ybi=$bEllipse-$delta+$yShift;
    $zai=-$cEllipse+$delta+$zShift; $zbi=$cEllipse-$delta+$zShift; 
    #
    $xai $xbi $yai $ybi $zai $zbi
  lines
    $nx = intmg( ($xbi-$xai)/$ds +1.5);
    $ny = intmg( ($ybi-$yai)/$ds +1.5);
    $nz = intmg( ($zbi-$zai)/$ds +1.5);
    $nx $ny $nz
  boundary conditions
    0 0 0 0 0 0
  mappingName
    $backGroundInner
  exit
#
$domainGrids[0] .= "\n $ellipsoidOuterGrid\n $southPoleOuterGrid\n $northPoleOuterGrid";
#
$domainName[$numDomains]="ellipsoid$numDomains";
$domainGrids[$numDomains] = "$backGroundInner\n $ellipsoidInnerGrid\n $southPoleInnerGrid\n $northPoleInnerGrid";
$numDomains = $numDomains + 1;
