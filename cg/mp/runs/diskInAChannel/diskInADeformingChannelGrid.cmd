#
# Grid for an rigid disk in a deforming channel
#
#
# usage: ogen [noplot] diskInADeformingChannelGrid -factor=<num> -order=[2/4/6/8] -interp=[e/i] 
# Options:
#    -ra,rb
# 
# examples:
#   ogen -noplot diskInADeformingChannelGrid -interp=e -factor=2
#
#
$prefix="diskInADeformingChannelGrid";
$order=2; $factor=1; $interp="e"; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids"; $ml=0; 
$name=""; 
$xa=-2; $xb=2; $ya=-1.; $yb=1.; 
$ra=.5; $rb=.75;  # rigid disk, inner and outer radii
$solidThickness=.5;  # thickness of the solid wall
$numberOfVolumeSmooths=0; 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"ra=f"=> \$ra,"rb=f"=> \$rb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"per=i"=>\$per,\
            "xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "dsx=f"=>\$dsx,"prefix=s"=>\$prefix );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $per eq 1 ){ $suffix .= "p"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# 
$ds=.1/$factor;
$pi = 4.*atan2(1.,1.);
# 
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
sub max{ local($n,$m)=@_; if( $n>$m ){ return $n; }else{ return $m; } }
#
$bcInterface=100;  # bc for interfaces
$shareInterface=100;        # share value for interfaces
#
create mappings
#
# --- make a line for the top interface from a nurbs ---
#
$cmd1=""; $degree=3; 
$nx = intmg( ($xb-$xa)/$ds + .5 );
$dx = ($xb-$xa)/($nx-1); 
for($i=0; $i<$nx; $i++){$x= $xa + $i*$dx; $y=$yb; $cmd1=$cmd1 . "$x $y\n";}
$arcLength=$xb-$xa; 
#
#  --- upper boundary 
  nurbs (curve)
    parameterize by index (uniform)
    enter points
    $nx $degree
    $cmd1
    lines
      $nx
    mappingName
      upperBoundary
    exit
# 
# -- Make a hyperbolic grid --
#
  $nr = intmg( 5 + $order/2 );
  hyperbolic
    forward
    # Fixed normal distance
    $nDistInterface=.2;
    # Make sure the number of grid cells doubles when ds is halved.
    $nDistExtraFactor=1.2; # make grid a bit finer in the normal direction
    $nrm = intmg(($nDistInterface*$nDistExtraFactor)/$ds+.5);
    # $nDistInterface=($nr-3)*$ds;
    # $nrm=$nr-1; 
    distance to march $nDistInterface
    lines to march $nrm
    points on initial curve $nx
    uniform dissipation 0.05
    volume smooths $numberOfVolumeSmooths
    equidistribution 0 (in [0,1])
    #
    BC: left fix x, float y and z
    BC: right fix x, float y and z
    #
    # spacing: geometric
    # geometric stretch factor 1.05 
    #
    generate
    # use fourth order interpolant to define the mapping: *wdh* May 8, 2018
    fourth order
    # evaluate as nurbs 1
    boundary conditions
      1 2 100 0 
    share 
      1 2 100 0 
    name fluidTopInterface
  exit
#
#  --- Solid wall on top ---
#
 rectangle
  set corners
    $yc = $yb + $solidThickness; 
    $xa $xb $yb $yc
  lines
    $ny = int( $nDistExtraFactor*($yc-$yb)/$ds + .5 ); 
    $nx $ny
  boundary conditions
    1 2 100 4 
  share 
    0 0 100 0 
  mappingName
    solidTopWall
  exit
#
#
# -- Rigid body in the fluid channel ----
#
Annulus
  inner and outer radii
    $innerRad=$ra; $outerRad=$rb; 
    $innerRad $outerRad
  lines
    $nTheta = intmg( 2.*$pi*($innerRad+$outerRad)*.5/$ds + 3.5 );
    $nr = intmg( $nDistExtraFactor*($outerRad-$innerRad)/$ds + 1.5 );
    $nTheta $nr
  boundary conditions
    -1 -1 5 0 
  mappingName
    fluidDisk
exit
#
#  --- fluid Channel Grid ---
#
 rectangle
  set corners
    # $ybb=$yb - $nDistInterface;
    $ybb=$yb+ .1; # add extra at top to allow for a displacement
    $xa $xb $ya $ybb
  lines
    $nx = int( ($xb -$xa)/$ds + .5 ); 
    $ny = int( ($ybb-$ya)/$ds + .5 ); 
    $nx $ny
  boundary conditions
    1 2 3 0
  share 
    1 2 0 0
  mappingName
    fluidChannel
  exit
#
exit this menu
#
generate an overlapping grid
  solidTopWall
  fluidChannel
  fluidTopInterface
  fluidDisk
  done choosing mappings
  # 
  change parameters
    specify a domain
     solidDomain
       solidTopWall
    done
    specify a domain
     fluidDomain
       fluidChannel
       fluidTopInterface
       fluidDisk
    done
    # choose implicit or explicit interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      2 2 2 2 2 2
    exit
  # open graphics
  compute overlap
exit
#
# save an overlapping grid
save a grid (compressed)
$name
elasticBalloonGrid
exit

