*
* Deforming sphere
*
* usage: ogen [noplot] sphereDeform -factor=<num> -order=[2/4/6/8] -interp=[e/i] -nrExtra=<>
*
*  nrExtra: extra lines to add in the radial direction on the sphere grids 
* 
* examples:
*     ogen noplot sphereDeform -factor=1 -order=2
*     ogen noplot sphereDeform -factor=1 -order=4
*     ogen noplot sphereDeform -factor=2 -order=4 -interp=e -nrExtra=8    (for cgmx : add extra grid lines in the normal direction)
* 
*     ogen noplot sphereDeform -factor=1 -order=2 -interp=e
*     ogen noplot sphereDeform -factor=2 -order=2 -interp=e
*     ogen noplot sphereDeform -factor=4 -order=2 -interp=e
*     ogen noplot sphereDeform -factor=8 -order=2 -interp=e
*     ogen noplot sphereDeform -factor=16 -order=2 -interp=e
*
$xa=-1.5; $xb=1.5; $ya=-1.5; $yb=1.5; $za=-1.5; $zb=1.5; $nrExtra=4; 
$order=2; $factor=1; $interp="i"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids"; $dse=0.; 
* 
* get command line arguments
GetOptions( "order=i"=>\$order,"factor=i"=> \$factor,"nrExtra=i"=> \$nrExtra,"interp=s"=> \$interp);
* 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; $dse=1.; }
* 
$suffix = ".order$order"; 
$name = "sphereDeform" . "$interp$factor" . $suffix . ".hdf";
* 
$ds=.2/$factor;
$rad = 1;  # radius of the sphere 
$sa0 = .5+$ds;  # cap extent parameter, 2=cap-covers half the sphere 
$deltaPhi=.05+$ds/4.;  # bounds on phi for the central sphere grid 
$pi = 4.*atan2(1.,1.);
* nr = number of lines in the radial direction
$nr=2+$order; if( $interp eq "e" ){ $nr=$nr+$order; } 
$innerRad=$rad; $outerRad=$innerRad+$distToMarch; 
$nr=$nr + $nrExtra; # make the grid finer in the normal direction (for stretching)
$linesToMarch=5; 
$distToMarch = $linesToMarch*$ds*0.5; 
*
$bcInterface0=100;  # bc for interfaces
$bcInterface1=100;  
$bcInterface2=100;  
$shareInterface=100;        # share value for interfaces
*
* 
* ---------------------------------------
** turn off graphics
* ---------------------------------------
*
create mappings
* first make grid for the surface of a sphere
Sphere
  inner and outer radii
    $innerRad $outerRad
  surface or volume (toggle)
  mappingName
    sphereSurface
exit
* 
* now make a mapping for the north pole
*
reparameterize
  transform which mapping?
    sphere
  orthographic
    * sa=2 --> patches just match (not including ghost points)
    $sa = $sa0; $sb=$sa; 
    specify sa,sb
      $sa $sb
  exit
  lines
    $nTheta=int( 2.5*$sa0*($innerRad+$outerRad)*.5/$ds +1.5 );    
    $nTheta $nTheta $nr
  boundary conditions
    0 0 0 0 1 0
  share
    0 0 0 0 1 0
  mappingName
    north-pole-surface-un-rotated
exit
* 
  rotate/scale/shift
    rotate
      -90 1
      0 0 0
    mappingName
      north-pole-surface
  exit
* 
  hyperbolic
    lines to march $linesToMarch
    distance to march $distToMarch
    BC: left free floating
    BC: right free floating
    BC: bottom free floating
    BC: top free floating
    march along normals 1
    apply boundary conditions to start curve 1
    backward
    boundary offset 1, 1, 1, 1, 0, 0 (l r b t b f)
    generate
    fourth order
    name north-pole
    boundary conditions
      0 0 0 0 $bcInterface0 0
    share
      0 0 0 0 $shareInterface 0
    * save the start curve in the data-base file: 
    save reference surface when put
    exit
*
* now make a mapping for the south pole
*
reparameterize
  transform which mapping?
    sphere
  orthographic
    choose north or south pole
      -1
    specify sa,sb
      $sa $sb
  exit
  lines
    $nTheta $nTheta $nr
  boundary conditions
    0 0 0 0 $bcInterface0 0
  share
    0 0 0 0 $shareInterface 0
  mappingName
    south-pole-surface-unrotated
exit
* 
  rotate/scale/shift
    rotate
      -90 1
      0 0 0
    mappingName
      sorth-pole-surface
  exit
* 
  hyperbolic
    lines to march $linesToMarch
    distance to march $distToMarch
    BC: left free floating
    BC: right free floating
    BC: bottom free floating
    BC: top free floating
    march along normals 1
    apply boundary conditions to start curve 1
    backward
    boundary offset 1, 1, 1, 1, 0, 0 (l r b t b f)
    generate
    fourth order
    name south-pole
    boundary conditions
      0 0 0 0 $bcInterface1 0
    share
      0 0 0 0 $shareInterface 0
    * save the start curve in the data-base file: 
    save reference surface when put
    exit
*
* Here is the box
*
Box
  set corners
    $xa $xb $ya $yb $za $zb
  lines
    $nx = int( ($xb-$xa)/$ds +1.5);
    $ny = int( ($yb-$ya)/$ds +1.5);
    $nz = int( ($zb-$za)/$ds +1.5);
    $nx $ny $nz
  boundary conditions
    0 0 0 0 0 0
  mappingName
    box
  exit
*
* Here is the main sphere grid 
*
 Sphere
  $nr=3+$order; if( $interp eq "e" ){ $nr=$nr+$order; } 
  $innerRad=$rad; $outerRad=$innerRad+($nr-1)*$ds;
  $nr=$nr + $nrExtra; 
  inner and outer radii
    $innerRad $outerRad
  bounds on phi (latitude)
    $phi0=$deltaPhi; $phi1=1.-$deltaPhi; 
    $phi0 $phi1 
  boundary conditions
    0 0 -1 -1 $bcInterface2 0
  share
    0 0  0  0 $shareInterface 0
  $nTheta = int( 2.*$pi*($innerRad+$outerRad)*.5/$ds +1.5 );
  $ns     = int( ($phi1-$phi0)*$pi*($innerRad+$outerRad)*.5/$ds +1.5 );
  lines 
    $ns $nTheta $nr 
  surface or volume (toggle)
  mappingName
    sphereNoEnds-unrotated 
 exit
* 
  rotate/scale/shift
    rotate
      -90 1
      0 0 0
    mappingName
      sphereNoEnds
  exit
* -- build a hyperbolic grid on the main sphere ---
  hyperbolic
    lines to march $linesToMarch
    distance to march $distToMarch
    BC: left free floating
    BC: right free floating
    * BC: bottom free floating
    * BC: top free floating
    march along normals 1
    apply boundary conditions to start curve 1
    backward
    boundary conditions
      0 0 -1 -1 $bcInterface2 0
    share
      0 0 0 0 $shareInterface 0
    generate
    name sphere
    * save the start curve in the data-base file: 
    save reference surface when put
    exit
*
exit
*
generate an overlapping grid
  box
  sphere
  north-pole
  south-pole
  done
  change parameters
    * improve quality of interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ng $ng $ng $ng $ng $ng 
  exit
* 
  compute overlap
* 
exit
* save an overlapping grid
save a grid (compressed)
$name
sphere
exit
