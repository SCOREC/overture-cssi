#
# Grid for a falling "stick" a channel
# 
# Examples:
#     ogen -noplot dropStickGrid -factor=1
#     ogen -noplot dropStickGrid -factor=1 -angle=10.
#
# -- bigger domain:
#      ogen -noplot dropStickGrid -factor=1 -xa=-5. -xb=5. -ya=-1.5 -yb=1. -name="dropStickGridBigDomain1.hdf" 
#      ogen -noplot dropStickGrid -factor=2 -xa=-4. -xb=4. -ya=-1.5 -yb=1. -name="dropStickGridBigDomain2.hdf" 
# 
#      ogen -noplot dropStickGrid -factor=2 -xa=-2.5 -xb=2.5 -ya=-1.5 -yb=2.0 -name="dropStickGrid2a.hdf" 
#      ogen -noplot dropStickGrid -factor=2 -xa=-2.5 -xb=2.5 -ya=-2.5 -yb=2.0 -name="dropStickGrid2b.hdf" 
#      ogen -noplot dropStickGrid -factor=2 -xa=-3.0 -xb=3.0 -ya=-2.5 -yb=2.0 -name="dropStickGrid2c.hdf" 
#
$prefix="dropStickGrid";  
$order=2; $factor=1; $interp="i"; $ml=0; # default values
$orderOfAccuracy = "second order"; $ng=2; $interpType = "implicit for all grids";
$xa =-1.; $xb=1.; $ya=-1.; $yb=1.; $angle=60.; 
$name=""; 
$cx=0.; $cy=0.;  # center for the stick 
# 
# get command line arguments
GetOptions( "order=i"=>\$order,"factor=f"=> \$factor,"xa=f"=>\$xa,"xb=f"=>\$xb,"ya=f"=>\$ya,"yb=f"=>\$yb,\
            "interp=s"=> \$interp,"name=s"=> \$name,"ml=i"=>\$ml,"blf=f"=> \$blf, "prefix=s"=> \$prefix,\
            "cx=f"=>\$cx,"cy=f"=>\$cy,"rgd=s"=> \$rgd,"angle=f"=>\$angle );
# 
if( $order eq 4 ){ $orderOfAccuracy="fourth order"; $ng=2; }\
elsif( $order eq 6 ){ $orderOfAccuracy="sixth order"; $ng=4; }\
elsif( $order eq 8 ){ $orderOfAccuracy="eighth order"; $ng=6; }
if( $interp eq "e" ){ $interpType = "explicit for all grids"; }
# 
$suffix = ".order$order"; 
if( $ml ne 0 ){ $suffix .= ".ml$ml"; }
if( $name eq "" ){$name = $prefix . "$interp$factor" . $suffix . ".hdf";}
# -- convert a number so that it is a power of 2 plus 1 --
#    ml = number of multigrid levels 
$ml2 = 2**$ml; 
sub intmg{ local($n)=@_; $n = int(int($n+$ml2-2)/$ml2)*$ml2+1; return $n; }
# 
$ds=.025/$factor;
$pi = 4.*atan2(1.,1.);
#
create mappings
#
rectangle
  set corners
    $xa $xb $ya $yb 
  lines
    $nx = intmg( ($xb-$xa)/$ds+1.5);  
    $ny = intmg( ($yb-$ya)/$ds+1.5);
    $nx $ny 
  boundary conditions
     1 2 3 4 
  mappingName
   channel
exit
#
SmoothedPolygon
# start on a side so that the polygon is symmetric
  vertices 
    6
    -.05  .00
    -.05  .25
     .05  .25
     .05 -.25
    -.05 -.25
    -.05  .00
  n-stretch
   1. 3.0 0.
  n-dist
    $nr = intmg( 6 );
    $nDist= ($nr-1)*$ds;
    fixed normal distance
    $nDist
  periodicity
    2
  lines
   $ns = intmg( 1.2*1.2/$ds + 1.5 );
   $nr = intmg( $nDist/$ds + 2.5 );
   $ns $nr 
  t-stretch
    0. 1.
    1. 5.
    1. 5.
    1. 5.
    1. 5.
    0. 1.
  boundary conditions
    -1 -1 5 0
  mappingName
    unrotated-stick
  exit
# 
  rotate/scale/shift
    transform which mapping?
     unrotated-stick
    rotate
     $angle
     0 0 0
    mappingName
     drop
  exit
#
exit
generate an overlapping grid
    channel
    drop
  done
  change parameters
 # choose implicit or explicit interpolation
    interpolation type
      $interpType
    order of accuracy 
      $orderOfAccuracy
    ghost points
      all
      $ngp = $ng+1;
      $ng $ng $ng $ngp $ng $ng
  exit
#  display intermediate results
  compute overlap
#  pause
  exit
#
save an overlapping grid
$name
dropStickGrid
exit

