#
#   Create grids for a deforming disk
# Input:
#      $count, $thickness, $chord, $theta, $cx, $cy, $r1, $r2, $x1, $x2
# $chord $rotation $cx $cy $r1 $x1 $r2 $x2
# $xCamber $nPoints $camber $toc $cutOff
  # $thickness: thickness of beam
  # $chord: length of beam
  # $theta: angle of beam
  # $cx, $cy: beam center
  # $x1, $r1: normalized core distance from leading edge, radius of core
#
$count = $count + 1;
$curveName="curve$count"; $fluidName="fluidInterface$count"; $share=99+$count; $solidName="solidInterface$count"; $coreName="core$count"; $solidBackGroundName="solidBackGround$count";
$coreBName="coreB$count";
$fluidDomain=$fluidDomain . $fluidName . "\n";
$solidGrids=$solidName . "\n" . $coreName . "\n" . $solidBackGroundName . "\n";
if( $r2 ne "" ){$solidGrids=$solidGrids . $coreBName . "\n";} else {$r2=0.; $x2=0.;}
$specifyDomainCmd=$specifyDomainCmd \
  . "specify a domain\nsolidDomain$count\n" \
  . $solidGrids \
  . "done\n";
$solidDomain=$solidDomain . $solidGrids;
$vertices="";
$sharpnessCmd="";
$tStretchCmd="";
for($i=0; $i <= $nPoints; $i++) { \
  $xoc=$cutOff*$i/$nPoints; \
  if ($xoc < $xCamber) { \
    $yc = $camber*$chord/($xCamber*$xCamber) \
      *(2*$xCamber*($xoc)-($xoc*$xoc)); \
    $dy = 2*$camber*$chord/($xCamber*$xCamber) \
      *($xCamber-$xoc); \
  } else { \
    $yc = $camber*$chord/((1-$xCamber)*(1-$xCamber))* \
      ((1-2*$xCamber)+2*$xCamber*($xoc)-($xoc*$xoc)); \
    $dy = 2*$camber*$chord/((1-$xCamber)*(1-$xCamber))* \
      ($xCamber-$xoc); \
  } \
  $theta = atan2($dy,1);                    \
  $yt = $toc*$chord/.2*(+0.2969*sqrt($xoc) \
                        -0.1260*$xoc \
                        -0.3516*($xoc)**2 \
                        +0.2843*($xoc)**3 \
                        -0.1015*($xoc)**4); \
  $x = $xoc*$chord-$yt*sin($theta); \
  $y = $yc+$yt*cos($theta); \
  $vertices=$vertices . "$x, $y\n"; \
  $sharpnessCmd=$sharpnessCmd . "$sharpness\n"; \
  if (($i eq 0) || ($i eq $nPoints)) { \
    $tStretchCmd=$tStretchCmd . "$tStretch, $tStretchWeight\n"; \
  } else { \
    $tStretchCmd=$tStretchCmd . "0, $tStretchWeight\n"; \
  } \
  if ($i eq $nPoints) { \
    $y = 0; \
    $vertices=$vertices . "$x, $y\n"; \
    $sharpnessCmd=$sharpnessCmd . "$sharpness\n"; \
    $tStretchCmd=$tStretchCmd . "$tStretch, $tStretchWeight\n"; \
  } \
 }
for($i=$nPoints; $i >= 0; $i--) { \
  $xoc=$cutOff*$i/$nPoints; \
  if ($xoc < $xCamber) { \
    $yc = $camber*$chord/($xCamber*$xCamber) \
      *(2*$xCamber*($xoc)-($xoc*$xoc)); \
    $dy = 2*$camber*$chord/($xCamber*$xCamber) \
      *($xCamber-$xoc); \
  } else { \
    $yc = $camber*$chord/((1-$xCamber)*(1-$xCamber))* \
      ((1-2*$xCamber)+2*$xCamber*($xoc)-($xoc*$xoc)); \
    $dy = 2*$camber*$chord/((1-$xCamber)*(1-$xCamber))* \
      ($xCamber-$xoc); \
  } \
  $theta = atan2($dy,1); \
  $yt = $toc*$chord/.2*(+0.2969*sqrt($xoc) \
                        -0.1260*$xoc \
                        -0.3516*($xoc)**2 \
                        +0.2843*($xoc)**3 \
                        -0.1015*($xoc)**4); \
  $x = $xoc*$chord+$yt*sin($theta); \
  $y = $yc-$yt*cos($theta); \
  $vertices=$vertices . "$x, $y\n"; \
  $sharpnessCmd=$sharpnessCmd . "$sharpness\n"; \
  if (($i eq 0) || ($i eq $nPoints)) { \
    $tStretchCmd=$tStretchCmd . "$tStretch, $tStretchWeight\n"; \
  } else { \
    $tStretchCmd=$tStretchCmd . "0, $tStretchWeight\n"; \
  } \
}
$totalVertices = 2*($nPoints+1)+1;
#
# --- make a nurbs for a circle ----
#   **** Later we could make different curves for the interface ****
  smoothedPolygon  
    # 
    $spacingFactor=1.2;  # add extra points around beam
    $nl=int($spacingFactor*(($toc*$chord+0.01)+($chord*2+0.02))/$ds);
    $nx=$nl;
#    $ny=int(7);
#    $ndist=($ny-1)*$ds;
    # wdh:
    $ny=int(7);
    $ndist=($ny-3)*$ds;
    periodicity
      2
    vertices 
      $totalVertices
      $vertices
    curve or area (toggle)
    lines 
    $nx
    sharpness
      $sharpnessCmd
    t-stretch
      $tStretchCmd
*    n-stretch
*      1 5 0
    mappingName
      $curveName
    exit 
  transform Mappings...
  rotate/scale/shift
    shift
    $cx $cy
    rotate
    $rotation
    $cx $cy
    exit
#
#  -- fluid beam grid ----
# 
  hyperbolic
    backward
    $dsFluid=.25*$toc*$chord*$ds*5;
    distance to march $ndist
    $ny = int($ndist/$dsFluid);
    lines to march $ny
    points on initial curve $nx
    # wdh: 
    # geometric stretch factor 1.1
    volume smooths 50
    # This next line is important to keep the ghost lines
    # near the corner of good quality:
    apply boundary conditions to start curve 1
    march along normals 1
    BC: top outward splay
    generate
    boundary conditions
      -1 -1 $share 0 0 0
    share 
       0 0 $share 0 0 0
    generate
    name $fluidName
  exit
#
# -- Solid ellipse near interface ----
#
  hyperbolic
    $nx=$nl;
    $nr0 = int(6);
    $solidFactor = $spacingFactor*(1. + 2./$factor); 
    $nr=int($nr0*$solidFactor);
    # wdh $ndist=$ds*$ny*0.1;
    # $ndist=$ds*$ny*0.175;
    $ndist=min( $toc*$chord*.1, ($nr0-3)*$ds);
    forward
    distance to march $ndist
    lines to march $nr
    points on initial curve $nx
    volume smooths 200
    # This next line is important to keep the ghost lines
    # near the corner of good quality:
    apply boundary conditions to start curve 1
    march along normals 1
    boundary conditions
      -1 -1 $share 0 0 0
    share 
       0 0 $share 0 0 0
    name $solidName
    generate
    exit
#
#   Solid annular core 
#
$xoc = $x1;
if ($xoc < $xCamber) {                       \
$yc = $camber*$chord/($xCamber*$xCamber)     \
  *(2*$xCamber*($xoc)-($xoc*$xoc)); \
} else { \
$yc = $camber*$chord/((1-$xCamber)*(1-$xCamber))* \
  ((1-2*$xCamber)+2*$xCamber*($xoc)-($xoc*$xoc)); \
}
$xcore=$xOffset*$chord+$cx+$x1*$chord*cos($rotation*$pi/180)-$yc*sin($rotation*$pi/180);
$ycore=$yOffset*$chord+$cy+$x1*$chord*sin($rotation*$pi/180)+$yc*cos($rotation*$pi/180);
  Annulus
  center: $xcore $ycore
  inner and outer radii
    $innerRad=$r1;
    $outerRad=$r1+$ndist;
    $innerRad $outerRad
  lines
    $nTheta = intmg( 2.*$pi*($toc*$chord)/$ds + 3.5 );
    $nDistExtraFactor=1.0;
    $nr = intmg( $nDistExtraFactor*($toc*$chord)/$ds + 1.5 );
    $nTheta $nr
  boundary conditions
    -1 -1 5 0 
  share
     0  0 5 0 
  mappingName
    $coreName
exit
#
#   Solid annular core 
#
  Annulus
    $xcore=$cx+cos($rotation*$pi/180)*$chord*($x2);
    $ycore=$cy+sin($rotation*$pi/180)*$chord*($x2)+$camber*$chord;
  center: $xcore $ycore
  inner and outer radii
    $innerRad=$r2;
    $outerRad=$r2+$ndist;
    $innerRad $outerRad
  lines
    $nTheta = intmg( 2.*$pi*($toc*$chord)/$ds + 3.5 );
    $nDistExtraFactor=1.0;
    $nr = intmg( $nDistExtraFactor*($toc*$chord)/$ds + 1.5 );
    $nTheta $nr
  boundary conditions
    -1 -1 5 0 
  share
     0  0 5 0 
  mappingName
    $coreBName
exit
#
#  --- solid backGround Grid ---
#
 rectangle
  set corners
    $nx=int( $solidFactor*$toc*$chord/$ds+1.5 );
    $ny=int( $solidFactor*$chord/$ds+1.5 );
    $nn=max($nx,$ny);
    $xas=$cx-.2; $xbs=$chord+$cx+.2; $yas=-$chord*.5+$cy; $ybs=$chord*.5+$cy; 
    $xas $xbs $yas $ybs
    $refineFactor=1.;
  lines
    $nn $nn
  boundary conditions
    0 0 0 0 
  share 
    0 0 0 0 
  mappingName
    $solidBackGroundName
  exit
