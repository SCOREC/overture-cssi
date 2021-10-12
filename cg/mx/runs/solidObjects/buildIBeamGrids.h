# ---------------------- IBEAM Grids -------------------------------
#   Input: 
#     $xShift=-1; $yShift=0; $angle=0; $xMin=-$radX; $xMax=-$xMin; $yMin=-$radY; $yMax=-$yMin;
#
# $xr=.0; $yr=0; # center for rotation
# $xShift=-1; $yShift=0; $angle=0; $xMin=-$radX; $xMax=-$xMin; $yMin=-$radY; $yMax=-$yMin;
$count=$count+1; $gridName="iBeam$count"; $gridNames .= "\n" . $gridName;
$mapName="iBeamGridBase";
# include transform.h
include $ENV{CG}/mx/runs/solidObjects/transformToNurbs.h
# 
# Inner background grid 
$gridName="iBeamSolidCore$count"; $iBeamNames .= "\n" . $gridName;
rectangle
    mappingName
      $gridName
   set corners
    $dss=$ds*.8; # inner grid spacing a bit smaller 
    # (radX,radY) : defined in iBeamSolid.h
    $xas=$xMin+$xShift; $xbs=$xMax+$xShift; 
    $yas=$yMin+$yShift; $ybs=$yMax+$yShift; 
    $xas $xbs $yas $ybs
   lines
    $nx = int( ($xbs-$xas)/$dss +1.5 ); 
    $ny = int( ($ybs-$yas)/$dss +1.5 ); 
    $nx $ny
    boundary conditions
      0 0 0 0
  exit
$mapName="iBeamGridBaseSolid";
$gridName="iBeamSolid$count"; $iBeamNames .= "\n" . $gridName;
# include transform.h
include $ENV{CG}/mx/runs/solidObjects/transformToNurbs.h
