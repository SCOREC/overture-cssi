#
#   Create grids for a deforming disk
# Input:
#      $count, $thickness, $chord, $theta, $cx, $cy, $r1, $r2, $x1, $x2
#
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
#
# --- make a nurbs for a circle ----
#   **** Later we could make different curves for the interface ****
  smoothedPolygon  
    $halfT=$thickness*0.5;
    $beamXa=-$chord*.5;
    $beamXb=$chord*.5;
    $beamYa=-$thickness*.5;
    $beamYb=+$thickness*.5;
    # 
    $spacingFactor=1.2;  # add extra points around beam
    $nl=int($spacingFactor*(($thickness+0.01)+($chord*2+0.02))/$ds);
    $nx=$nl;
#    $ny=int(7);
#    $ndist=($ny-1)*$ds;
    # wdh:
    $ny=int(7);
    $ndist=($ny-3)*$ds;
    periodicity
      2
    vertices 
      7
      $beamXb,0
      $beamXb,$beamYa 
      $beamXa,$beamYa
      $beamXa,0
      $beamXa,$beamYb
      $beamXb,$beamYb
      $beamXb,0
    curve or area (toggle)
    lines 
    $nx
    sharpness
      $sharpness
      $sharpness
      $sharpness
      $sharpness
      $sharpness
      $sharpness
      $sharpness
    t-stretch
      0.0, 10.
      0.1, $tStretch
      0.1, $tStretch
      0.1, $tStretch
      0.1, $tStretch
      0.1, $tStretch
      0.0, 10.
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
    $theta
    $cx $cy
    exit
#
#  -- fluid beam grid ----
# 
  hyperbolic
    backward
    distance to march $ndist
    lines to march $ny
    points on initial curve $nx
    # wdh: 
    # geometric stretch factor 1.1
    volume smooths 50
    # This next line is important to keep the ghost lines
    # near the corner of good quality:
    apply boundary conditions to start curve 1
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
    $ndist=min( $thickness*.15, ($nr0-3)*$ds);
    forward
    distance to march $ndist
    lines to march $nr
    points on initial curve $nx
    volume smooths 200
    # This next line is important to keep the ghost lines
    # near the corner of good quality:
    apply boundary conditions to start curve 1
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
  Annulus
    $xcore=$cx+cos($theta*$pi/180)*$chord*($x1-0.5);
    $ycore=$cy+sin($theta*$pi/180)*$chord*($x1-0.5);
  center: $xcore $ycore
  inner and outer radii
    #$ndist=min( $thickness*.15, ($nr0-3)*$ds);
    $innerRad=$r1;
    $outerRad=$r1+$ndist;
    $innerRad $outerRad
  lines
    $nTheta = intmg( 2.*$pi*($thickness)/$ds + 3.5 );
    $nDistExtraFactor=1.0;
    $nr = intmg( $nDistExtraFactor*($thickness)/$ds + 1.5 );
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
    $xcore=$cx+cos($theta*$pi/180)*$chord*($x2-0.5);
    $ycore=$cy+sin($theta*$pi/180)*$chord*($x2-0.5);
  center: $xcore $ycore
  inner and outer radii
    #$ndist=min( $thickness*.15, ($nr0-3)*$ds);
    $innerRad=$r2;
    $outerRad=$r2+$ndist;
    $innerRad $outerRad
  lines
    $nTheta = intmg( 2.*$pi*($thickness)/$ds + 3.5 );
    $nDistExtraFactor=1.0;
    $nr = intmg( $nDistExtraFactor*($thickness)/$ds + 1.5 );
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
    $nx=int( $solidFactor*$thickness/$ds+1.5 );
    $ny=int( $solidFactor*$chord/$ds+1.5 );
    $nn=max($nx,$ny);
    $xas=-$chord*.5+$cx; $xbs=$chord*.5+$cx; $yas=-$chord*.5+$cy; $ybs=$chord*.5+$cy; 
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
