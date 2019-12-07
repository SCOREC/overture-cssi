#
#  ----- Given arrays of points xc,yc, build curves and grids -----
#
# Update lists 
$count +=1; $shareInterface += 1; $bcInterface+=1;
$fluidGridName="fluidInterface$count";
$solidGridName="solidInterface$count";
$solidBackGroundName="solidBackGround$count";
#
$solidCoreName="solidCore$count";
## $solidGridNames .= "\n $solidBackGroundName\n $solidCoreName\n $solidGridName"; 
# -- remove the core grid 
$solidGridNames .= "\n $solidBackGroundName\n$solidGridName"; 
#
$fluidGridNames .= "\n $fluidGridName";
## $solidDomains .= "\n specify a domain\n solidDomain$count\n $solidBackGroundName\n $solidCoreName\n $solidGridName\n done";
$solidDomains .= "\n specify a domain\n solidDomain$count\n $solidBackGroundName\n $solidGridName\n done";
#
if( $degree eq "" ) { $degree=3; } # degree of Nurbs 
if( $numberOfVolumeSmooths eq "" ){ $numberOfVolumeSmooths=20; }  #
# rotation angle (degress)
if( $angle eq "" ){ $angle=0.; }  
if( $pi eq "" ){ $pi=4.*atan2(1.,1.); }  
# Core offset from center
if( $xCore eq "" ){ $xCore=0.; }
if( $yCore eq "" ){ $yCore=0.; }
if( $coreRadius eq "" ){ $coreRadius=.05; }
#
$xcMin=1e10; $xcMax=-1e10; $ycMin=1e10; $ycMax=-1e10;  # for inner background grid 
$cmd="#";  $ns=0;  $arcLength=0.; $x0=$xc[0]; $y0=$yc[0];
for( $ic=0; $ic<$nc-1; $ic++ ){ $numpt=$numPerSegment; if( $ic == $nc-2 ){ $numpt=$numPerSegment+1;} \
for( $i=0; $i<$numpt; $i++ ){ $s=($i)/($numPerSegment); $ns=$ns+1;  \
   $xs=(1.-$s)*$xc[$ic]+$s*$xc[$ic+1]; \
   $ys=(1.-$s)*$yc[$ic]+$s*$yc[$ic+1];   \
   $ct=cos($angle*$pi/180.); $st=sin($angle*$pi/180.);  \
   $x=$cx + $ct*($xs-$cx)-$st*($ys-$cy);   \
   $y=$cy + $st*($xs-$cx)+$ct*($ys-$cy);   \
   $xcMin = min($xcMin,$x); $xcMax = max($xcMax,$x);  \
   $ycMin = min($ycMin,$y); $ycMax = max($ycMax,$y);  \
   $arcLength=$arcLength + sqrt( ($x-$x0)**2 + ($y-$y0)**2 ); printf(" ns=$ns: x0=$x0 x=$x y0=$y0 y=$y arcLength=$arcLength\n");  $x0=$x; $y0=$y; \
   $cmd .= "\n $x $y 1."; } }\
  $knots="#"; for( $i=$degree-1; $i<$ns-($degree-1); $i++ ){ $s=$i/($ns-2); $knots .= "\n $s"; } \
#
nurbs (curve)
  periodicity
   2
  enter control points
    $degree
    $ns
    $knots
    $cmd 
 #
 parameterize by chord length
 #
 lines
  $lines=intmg($arcLength/$ds + 1.5 );
  $lines
 mappingName
  $curveName="boundaryInitial$count";
  $curveName
exit
# -- interpolate the initial NURBS so that we have an arclength parameterization --
nurbs (curve)
  interpolate from a mapping
    $curveName
  mappingName
   $curveName2="boundaryCurve$count";
   $curveName2
exit 
# 
# -- Make EXTERIOR hyperbolic grid --
#
  $nr = intmg( 6 + 3*($order-2)/2 );
  hyperbolic
    forward
    $nDist=($nr-3.25)*$ds;
    distance to march $nDist
    $nrm=$nr-1; 
    lines to march $nrm
    $nTheta = int($arcLength/$ds+1.5);
    points on initial curve $nTheta
    uniform dissipation 0.05
    volume smooths $numberOfVolumeSmooths
    equidistribution 0 (in [0,1])
    #
    spacing: geometric
    geometric stretch factor 1.05 
    #
    generate
    # open graphics
    # pause
    boundary conditions
      -1 -1 $bcInterface 0 0 0
    share 
       0 0 $shareInterface 0 0 0
    name $fluidGridName
  exit
# 
# -- Make INTERIOR hyperbolic grid --
#
   $dsi=$ds*.75;  # target grid spacing for interior
# 
  $nr = intmg( 6 + 3*($order-2)/2 );
  hyperbolic
    backward
    $nDist=($nr-3.25)*$dsi;
    distance to march $nDist
    $nrm=$nr-1; 
    lines to march $nrm
    $nTheta = int($arcLength/$ds+1.5);
    points on initial curve $nTheta
    uniform dissipation 0.05
    volume smooths 20 
    equidistribution 0 (in [0,1])
    #
    spacing: geometric
    geometric stretch factor 1.05 
    #
    generate
    # open graphics
    # pause
    boundary conditions
      -1 -1 $bcInterface 0 0 0
    share 
       0 0 $shareInterface 0 0 0
    name $solidGridName
  exit
#
#   Inner background -- make a bit finer than target
#
   rectangle
   set corners
    $xai=$xcMin; $xbi=$xcMax;
    $yai=$ycMin; $ybi=$ycMax;
    $xai $xbi $yai $ybi
   lines
    $nxi = int( ($xbi-$xai)/$dsi +1.5 ); 
    $nyi = int( ($ybi-$yai)/$dsi +1.5 ); 
    $nxi $nyi
    boundary conditions
      0 0 0 0 
    mappingName
      $solidBackGroundName
  exit 
#
#  Inner hole (core)
#
#
#   Solid annular core *** NOT CURRENTLY USED *****
#
#-   $coreName="solidCore";
#-   Annulus
#-   #
#-   $x=$cx + $ct*($xCore)-$st*($yCore); 
#-   $y=$cy + $st*($xCore)+$ct*($yCore);  
#-   center: $x $y
#-   inner and outer radii
#-     $innerRad=$coreRadius;
#-     $nr =  intmg( 6 );
#-     $nDist = ($nr-2)*$dsi;
#-     $outerRad=$coreRadius+$nDist;
#-     $innerRad $outerRad
#-   lines
#-     $nTheta = intmg( $pi*($innerRad+$outerRad)/$dsi + .5 );
#-     $nDistExtraFactor=1.0;
#-     # $nr = intmg( $nDistExtraFactor*($outerRad-$innerRad)/$dsi + .5 );
#-     $nTheta $nr
#-   boundary conditions
#-     -1 -1 5 0 
#-   share
#-      0  0 5 0 
#-   mappingName
#-     $solidCoreName
#- exit
#
